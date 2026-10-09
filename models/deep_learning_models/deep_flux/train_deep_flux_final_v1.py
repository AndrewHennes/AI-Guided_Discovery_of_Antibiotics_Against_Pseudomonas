#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Train GI60 and GI80 ensembles and score molecular representations.

Set the data paths and run settings in the configuration block below, or pass
them as command-line arguments. Use absolute paths for inputs and outputs.
The training CSV must contain Minimol_Representation and the two configured
binary target columns. Representations must be serialized lists of numbers;
negative target values denote missing labels. The optional screening CSV must
contain Minimol_Representation. Set SCREENING_DATA_PATH to None to train only.

RUN_HYPERPARAMETER_SEARCH enables ten-fold cross-validation for GI60, GI80,
and their joint model before ensemble training. Outputs include tuning results,
the best ensemble checkpoints, and scored screening CSVs under OUTPUT_ROOT.
"""

from __future__ import annotations

import argparse
import ast
import json
import os
from pathlib import Path
import tempfile
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import roc_auc_score, average_precision_score

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

import pytorch_lightning as pl
from pytorch_lightning.callbacks import EarlyStopping, ModelCheckpoint

import ray
from ray import tune, train
from ray.tune.search.optuna import OptunaSearch
from ray.tune.callback import Callback as RayTuneCallback

##############################################################################
# 1) Read Data
##############################################################################

def _strip(s: str) -> str:
    return s.strip() if isinstance(s, str) else s

def read_data(file_path: str, rep_col_name: str, target_col_names: Sequence[str]) -> pd.DataFrame:
    df = pd.read_csv(file_path, sep=",")
    # Robust to stray tabs/spaces in header names
    df.columns = [c.strip() for c in df.columns]

    rep_col_name = _strip(rep_col_name)
    target_col_names = [_strip(c) for c in target_col_names]

    print("Columns in the DataFrame are:", df.columns.tolist())

    required = [rep_col_name] + list(target_col_names)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(
            "Missing expected columns in CSV after stripping whitespace:\n"
            f"  Missing: {missing}\n"
            f"  Available: {df.columns.tolist()}"
        )

    df = df[required].copy()
    return df

##############################################################################
# 2) Split Data
##############################################################################

def train_val_test_split(
    df: pd.DataFrame,
    split_fracs: Tuple[float, float, float] = (0.85, 0.1, 0.05),
    random_state: int = 42
):
    assert len(split_fracs) == 3, "split_fracs must have length 3."
    assert abs(sum(split_fracs) - 1.0) < 1e-7, "split_fracs must sum to 1."
    train_frac, val_frac, test_frac = split_fracs

    df_train, df_temp = train_test_split(
        df, test_size=(1.0 - train_frac), random_state=random_state
    )
    relative_val_frac = val_frac / (val_frac + test_frac)
    df_val, df_test = train_test_split(
        df_temp, test_size=(1.0 - relative_val_frac), random_state=random_state
    )
    return df_train, df_val, df_test

##############################################################################
# 3) DataModule
##############################################################################

def parse_feature_tensor(rep_values: Sequence[object]) -> torch.Tensor:
    parsed_features = []
    for row_index, item in enumerate(rep_values):
        try:
            array = np.asarray(
                ast.literal_eval(item) if isinstance(item, str) else item, dtype=np.float32
            )
        except (ValueError, SyntaxError, TypeError, OverflowError) as error:
            raise ValueError(f"Invalid feature representation at row {row_index}.") from error
        if array.ndim != 1 or array.size == 0 or not np.isfinite(array).all():
            raise ValueError(f"Features at row {row_index} must be a nonempty, finite vector.")
        if parsed_features and array.shape != parsed_features[0].shape:
            raise ValueError(f"Feature vector length at row {row_index} differs from the first row.")
        parsed_features.append(array)
    if not parsed_features:
        raise ValueError("The dataset must contain at least one feature vector.")
    return torch.from_numpy(np.stack(parsed_features, axis=0))

class MultiTaskDataset(Dataset):
    def __init__(self, df: pd.DataFrame, rep_col_name: str, target_col_names: Sequence[str]):
        super().__init__()
        if len(target_col_names) < 1:
            raise ValueError("target_col_names must contain at least one target column name.")

        self.target_col_names = list(target_col_names)
        self.features = parse_feature_tensor(df[rep_col_name])
        # shape: (N, num_tasks)
        self.targets = torch.from_numpy(df[self.target_col_names].values.astype(np.float32))

    def __len__(self):
        return len(self.features)

    def __getitem__(self, index: int):
        return {"x": self.features[index], "y": self.targets[index]}

class MultiTaskDataModule(pl.LightningDataModule):
    def __init__(
        self,
        df_train: pd.DataFrame,
        df_val: pd.DataFrame,
        df_test: pd.DataFrame,
        rep_col_name: str,
        target_col_names: Sequence[str],
        batch_size: int = 32,
        df_extra_test: Optional[pd.DataFrame] = None,
    ):
        super().__init__()
        self.df_train = df_train
        self.df_val = df_val
        self.df_test = df_test
        self.rep_col_name = rep_col_name
        self.target_col_names = list(target_col_names)
        self.batch_size = batch_size
        self.df_extra_test = df_extra_test
        self.extra_test_dataset = None

    def setup(self, stage: Optional[str] = None):
        self.train_dataset = MultiTaskDataset(
            self.df_train, self.rep_col_name, self.target_col_names
        )
        self.val_dataset = MultiTaskDataset(
            self.df_val, self.rep_col_name, self.target_col_names
        )
        self.test_dataset = MultiTaskDataset(
            self.df_test, self.rep_col_name, self.target_col_names
        )
        if self.df_extra_test is not None and len(self.df_extra_test) > 0:
            self.extra_test_dataset = MultiTaskDataset(
                self.df_extra_test, self.rep_col_name, self.target_col_names
            )

    def train_dataloader(self):
        batch_size = min(self.batch_size, len(self.train_dataset))
        return DataLoader(
            self.train_dataset,
            batch_size=batch_size,
            shuffle=True,
            drop_last=True,
        )

    def val_dataloader(self):
        return DataLoader(self.val_dataset, batch_size=self.batch_size)

    def test_dataloader(self):
        return DataLoader(self.test_dataset, batch_size=self.batch_size)

    def extra_test_dataloader(self):
        if self.extra_test_dataset is not None:
            return DataLoader(self.extra_test_dataset, batch_size=self.batch_size)
        return None

##############################################################################
# 4) Model
##############################################################################

class MLPLayer(nn.Module):
    """Linear -> (BatchNorm) -> Activation -> (Dropout)."""

    def __init__(self, in_dim, out_dim, activation, use_batch_norm=False, dropout_rate=0.0):
        super().__init__()
        self.use_batch_norm = use_batch_norm
        self.activation = activation

        self.linear = nn.Linear(in_dim, out_dim)
        self.bn = nn.BatchNorm1d(out_dim) if self.use_batch_norm else None
        self.dropout = nn.Dropout(p=dropout_rate) if dropout_rate > 0 else nn.Identity()

    def forward(self, x):
        x = self.linear(x)
        if self.bn is not None:
            x = self.bn(x)
        x = self.activation(x)
        x = self.dropout(x)
        return x

class MultiTaskFeedForward(pl.LightningModule):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        dim_size: int = 128,
        dim_shrinking_scale: float = 0.5,
        num_layers: int = 3,
        learning_rate: float = 1e-3,
        batch_size: int = 32,
        L2_weight_norm: float = 0.0,
        L1_weight_norm: float = 0.0,
        activation_function: str = "relu",
        dropout_rate: float = 0.0,
        use_residual: bool = False,
        use_batch_norm: bool = False,
        scheduler_step_size: int = 5,
        scheduler_gamma: float = 0.5,
        task1_loss_scale: float = 1.0,
    ):
        super().__init__()
        self.save_hyperparameters()

        if output_dim < 1:
            raise ValueError("output_dim must be >= 1.")

        self.input_dim = input_dim
        self.output_dim = output_dim

        self.dim_size = dim_size
        self.dim_shrinking_scale = dim_shrinking_scale
        self.num_layers = num_layers
        self.learning_rate = learning_rate
        self.batch_size = batch_size
        self.L2_weight_norm = L2_weight_norm
        self.L1_weight_norm = L1_weight_norm
        self.activation_function = activation_function
        self.dropout_rate = dropout_rate
        self.use_residual = use_residual
        self.use_batch_norm = use_batch_norm

        self.scheduler_step_size = scheduler_step_size
        self.scheduler_gamma = scheduler_gamma

        self.task1_loss_scale = task1_loss_scale

        activations = {
            "relu": nn.ReLU(),
            "tanh": nn.Tanh(),
            "sigmoid": nn.Sigmoid(),
            "leaky_relu": nn.LeakyReLU(),
        }
        self.act_fn = activations.get(self.activation_function, nn.ReLU())

        self.shared_layers = self._build_shared_layers()
        final_dim = self._get_final_dim()
        self.head = nn.Linear(final_dim, self.output_dim)

    def _build_shared_layers(self):
        layers = nn.ModuleList()
        in_dim = self.input_dim
        current_dim = self.dim_size

        for _ in range(self.num_layers):
            layers.append(
                MLPLayer(
                    in_dim=in_dim,
                    out_dim=current_dim,
                    activation=self.act_fn,
                    use_batch_norm=self.use_batch_norm,
                    dropout_rate=self.dropout_rate,
                )
            )
            in_dim = current_dim
            current_dim = max(1, int(current_dim * self.dim_shrinking_scale))

        return layers

    def _get_final_dim(self):
        dim = self.dim_size
        for _ in range(1, self.num_layers):
            dim = max(1, int(dim * self.dim_shrinking_scale))
        return dim

    def forward(self, x):
        out = self._forward_shared_layers(x, self.shared_layers)
        return self.head(out)

    def _forward_shared_layers(self, x, layers):
        out = x
        for layer in layers:
            new_out = layer(out)
            if self.use_residual and (new_out.shape == out.shape):
                out = out + new_out
            else:
                out = new_out
        return out

    def _task_loss_scales(self) -> List[float]:
        return [self.task1_loss_scale] + [1.0] * (self.output_dim - 1)

    def _compute_loss(self, logits, y):
        """
        Multi-task BCE-with-logits with partial labels:
          - For each task t, samples with y[:, t] < 0 are ignored for that task.
          - Task losses are scaled by task{t}_loss_scale hyperparameters.
        """
        if logits.ndim != 2:
            raise ValueError(f"Expected logits shape (N, T), got {tuple(logits.shape)}")
        if y.ndim != 2:
            raise ValueError(f"Expected y shape (N, T), got {tuple(y.shape)}")
        if logits.shape[1] != y.shape[1]:
            raise ValueError(
                f"Logits/task dim mismatch: logits has {logits.shape[1]} tasks, y has {y.shape[1]} tasks."
            )

        scales = self._task_loss_scales()
        total = None
        any_valid = False

        for t in range(y.shape[1]):
            mask = y[:, t] >= 0.0
            if not mask.any():
                continue
            any_valid = True
            task_logits = logits[mask, t]
            task_y = y[mask, t]
            task_loss = F.binary_cross_entropy_with_logits(task_logits, task_y, reduction="mean")
            task_loss = task_loss * float(scales[t])
            total = task_loss if total is None else (total + task_loss)

        if not any_valid:
            return torch.tensor(0.0, device=logits.device, requires_grad=True)
        return total

    def training_step(self, batch, batch_idx):
        x, y = batch["x"], batch["y"]
        logits = self(x)
        loss = self._compute_loss(logits, y)

        if self.L1_weight_norm > 0:
            l1_penalty = sum(param.abs().sum() for param in self.parameters())
            loss = loss + self.L1_weight_norm * l1_penalty

        self.log("train_loss", loss)
        return loss

    def validation_step(self, batch, batch_idx):
        x, y = batch["x"], batch["y"]
        logits = self(x)
        val_loss = self._compute_loss(logits, y)
        self.log("val_loss", val_loss, prog_bar=True)
        return val_loss

    def test_step(self, batch, batch_idx):
        x, y = batch["x"], batch["y"]
        logits = self(x)
        test_loss = self._compute_loss(logits, y)
        self.log("test_loss", test_loss, prog_bar=True)
        return test_loss

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(
            self.parameters(),
            lr=self.learning_rate,
            weight_decay=self.L2_weight_norm,
        )
        scheduler_config = {
            "scheduler": torch.optim.lr_scheduler.StepLR(
                optimizer,
                step_size=self.scheduler_step_size,
                gamma=self.scheduler_gamma,
            ),
            "interval": "epoch",
            "frequency": 1,
        }
        return [optimizer], [scheduler_config]

##############################################################################
# 5) Metrics
##############################################################################

def compute_auroc_auprc(model: nn.Module, dataloader: DataLoader, device: torch.device):
    """
    Computes AUROC/AUPRC per task (column) and returns:
      - auroc_per_task, auprc_per_task (lists)
      - auroc, auprc (nan-mean across tasks)
    """
    model.eval()
    all_logits = []
    all_targets = []

    with torch.no_grad():
        for batch in dataloader:
            x = batch["x"].to(device)
            y = batch["y"].to(device)
            logits = model(x)
            all_logits.append(logits.detach().cpu())
            all_targets.append(y.detach().cpu())

    logits = torch.cat(all_logits, dim=0)  # (N, T)
    targets = torch.cat(all_targets, dim=0)  # (N, T)

    probs = torch.sigmoid(logits).numpy()
    targets_np = targets.numpy()

    num_tasks = probs.shape[1]
    auroc_per_task: List[float] = []
    auprc_per_task: List[float] = []

    for t in range(num_tasks):
        y_true = targets_np[:, t]
        y_pred = probs[:, t]

        valid_mask = y_true >= 0.0
        y_true = y_true[valid_mask]
        y_pred = y_pred[valid_mask]

        if len(y_true) == 0:
            auroc_val, auprc_val = float("nan"), float("nan")
        else:
            try:
                auroc_val = roc_auc_score(y_true, y_pred)
            except ValueError:
                auroc_val = float("nan")
            try:
                auprc_val = average_precision_score(y_true, y_pred)
            except ValueError:
                auprc_val = float("nan")

        auroc_per_task.append(float(auroc_val))
        auprc_per_task.append(float(auprc_val))

    auroc_avg = float(np.nanmean(auroc_per_task)) if np.isfinite(auroc_per_task).any() else float("nan")
    auprc_avg = float(np.nanmean(auprc_per_task)) if np.isfinite(auprc_per_task).any() else float("nan")

    print(f"auROCs per task: {auroc_per_task} | avg: {auroc_avg}")
    print(f"auPRCs per task: {auprc_per_task} | avg: {auprc_avg}")

    return {
        "auroc": auroc_avg,
        "auprc": auprc_avg,
        "auroc_per_task": auroc_per_task,
        "auprc_per_task": auprc_per_task,
    }

##############################################################################
# 6) train_with_tune
##############################################################################

def fit_best_model(
    model: MultiTaskFeedForward,
    data_module: MultiTaskDataModule,
    max_epochs: int,
    patience: int,
    monitor: str,
    checkpoint_directory: Optional[str],
    checkpoint_path: Optional[str] = None,
):
    if model.use_batch_norm and min(data_module.batch_size, len(data_module.df_train)) < 2:
        raise ValueError("Batch normalization requires at least two training examples per batch.")
    if checkpoint_directory is not None:
        Path(checkpoint_directory).mkdir(parents=True, exist_ok=True)
    mode = "min" if monitor == "val_loss" else "max"

    with tempfile.TemporaryDirectory(dir=checkpoint_directory, prefix=".training-") as temporary:
        early_stopping = EarlyStopping(monitor=monitor, patience=patience, mode=mode, verbose=False)
        checkpoint_callback = ModelCheckpoint(
            dirpath=temporary,
            filename="best",
            monitor=monitor,
            mode=mode,
            save_top_k=1,
        )
        trainer = pl.Trainer(
            max_epochs=max_epochs,
            callbacks=[early_stopping, checkpoint_callback],
            enable_checkpointing=True,
            enable_progress_bar=False,
            logger=False,
            default_root_dir=temporary,
        )
        trainer.fit(model, data_module)

        best_path = checkpoint_callback.best_model_path
        if not best_path:
            raise RuntimeError(f"Training did not produce a checkpoint for {monitor!r}.")
        checkpoint = torch.load(best_path, map_location=model.device, weights_only=False)
        model.load_state_dict(checkpoint["state_dict"])
        if checkpoint_path is not None:
            os.replace(best_path, checkpoint_path)

    return model


def train_with_tune(
    config,
    df_train: pd.DataFrame,
    df_val: pd.DataFrame,
    df_test: pd.DataFrame,
    rep_col: str,
    target_cols: Sequence[str],
    max_epochs: int = 10,
    df_extra_test: Optional[pd.DataFrame] = None,
):
    df_cross_validation = pd.concat([df_train, df_val], ignore_index=True)

    fold_auroc_avg_values: List[float] = []
    fold_auprc_avg_values: List[float] = []
    fold_auroc_per_task_values: List[List[float]] = []
    fold_auprc_per_task_values: List[List[float]] = []

    fold_splitter = KFold(n_splits=10, shuffle=True, random_state=42)

    for fold_train_indices, fold_val_indices in fold_splitter.split(df_cross_validation):
        fold_train_dataframe = df_cross_validation.iloc[fold_train_indices]
        fold_val_dataframe = df_cross_validation.iloc[fold_val_indices]

        data_module = MultiTaskDataModule(
            fold_train_dataframe,
            fold_val_dataframe,
            fold_val_dataframe,
            rep_col_name=rep_col,
            target_col_names=target_cols,
            batch_size=int(config["batch_size"]),
            df_extra_test=df_extra_test,
        )
        data_module.setup()

        input_dim = data_module.train_dataset.features.shape[1]
        output_dim = len(list(target_cols))

        model = MultiTaskFeedForward(
            input_dim=input_dim,
            output_dim=output_dim,
            dim_size=int(config["dim_size"]),
            dim_shrinking_scale=float(config["dim_shrinking_scale"]),
            num_layers=int(config["num_layers"]),
            learning_rate=float(config["learning_rate"]),
            batch_size=int(config["batch_size"]),
            L2_weight_norm=float(config["L2_weight_norm"]),
            L1_weight_norm=float(config["L1_weight_norm"]),
            activation_function=str(config["activation_function"]),
            dropout_rate=float(config["dropout_rate"]),
            use_residual=bool(config["use_residual"]),
            use_batch_norm=bool(config["use_batch_norm"]),
            scheduler_step_size=int(config["scheduler_step_size"]),
            scheduler_gamma=float(config["scheduler_gamma"]),
            task1_loss_scale=float(config.get("task1_loss_scale", 1.0)),
        )

        model = fit_best_model(
            model=model,
            data_module=data_module,
            max_epochs=max_epochs,
            patience=int(config["patience"]),
            monitor=config.get("early_stop_metric", "val_loss"),
            checkpoint_directory=config.get("model_file_path"),
        )

        fold_stats = compute_auroc_auprc(model, data_module.val_dataloader(), model.device)

        fold_auroc_avg_values.append(float(fold_stats["auroc"]))
        fold_auprc_avg_values.append(float(fold_stats["auprc"]))
        fold_auroc_per_task_values.append(list(fold_stats["auroc_per_task"]))
        fold_auprc_per_task_values.append(list(fold_stats["auprc_per_task"]))

    auroc_avg = float(np.nanmean(fold_auroc_avg_values))
    auroc_std = float(np.nanstd(fold_auroc_avg_values))
    auprc_avg = float(np.nanmean(fold_auprc_avg_values))
    auprc_std = float(np.nanstd(fold_auprc_avg_values))

    auroc_matrix = np.array(fold_auroc_per_task_values, dtype=float)  # (folds, tasks)
    auprc_matrix = np.array(fold_auprc_per_task_values, dtype=float)

    val_auroc_per_task_mean = np.nanmean(auroc_matrix, axis=0).astype(float).tolist()
    val_auroc_per_task_std = np.nanstd(auroc_matrix, axis=0).astype(float).tolist()
    val_auprc_per_task_mean = np.nanmean(auprc_matrix, axis=0).astype(float).tolist()
    val_auprc_per_task_std = np.nanstd(auprc_matrix, axis=0).astype(float).tolist()

    report = tune.report if hasattr(tune, "get_context") else train.report
    report(
        {
            "val_auprc_avg": auprc_avg,  # objective metric for OptunaSearch
            "val_auroc_avg": auroc_avg,
            "val_auroc_std": auroc_std,
            "val_auprc_std": auprc_std,
            "val_auroc_per_task_mean": val_auroc_per_task_mean,
            "val_auroc_per_task_std": val_auroc_per_task_std,
            "val_auprc_per_task_mean": val_auprc_per_task_mean,
            "val_auprc_per_task_std": val_auprc_per_task_std,
        }
    )

##############################################################################
# 6b) Ensemble training + deployment
##############################################################################

def build_model_from_config(config: Dict[str, object], input_dim: int, output_dim: int):
    model_kwargs = {
        key: config[key]
        for key in [
            "dim_size",
            "dim_shrinking_scale",
            "num_layers",
            "learning_rate",
            "batch_size",
            "L2_weight_norm",
            "L1_weight_norm",
            "activation_function",
            "dropout_rate",
            "use_residual",
            "use_batch_norm",
            "scheduler_step_size",
            "scheduler_gamma",
            "task1_loss_scale",
        ]
        if key in config
    }
    return MultiTaskFeedForward(input_dim=input_dim, output_dim=output_dim, **model_kwargs)

def train_ensemble_models(
    config: Dict[str, object],
    df_train: pd.DataFrame,
    df_val: pd.DataFrame,
    df_test: pd.DataFrame,
    rep_col_name: str,
    target_col_names: Sequence[str],
    number_ensemble_members: int,
    ensemble_members_path: str,
    max_epochs: int = 25,
    df_extra_test: Optional[pd.DataFrame] = None,
    seed: int = 42,
) -> List[str]:
    if number_ensemble_members < 1:
        raise ValueError("number_ensemble_members must be at least 1.")
    ensemble_members_path = Path(ensemble_members_path)
    data_module = MultiTaskDataModule(
        df_train,
        df_val,
        df_test,
        rep_col_name=rep_col_name,
        target_col_names=target_col_names,
        batch_size=int(config["batch_size"]),
        df_extra_test=df_extra_test,
    )
    data_module.setup()

    input_dim = data_module.train_dataset.features.shape[1]
    output_dim = len(list(target_col_names))

    os.makedirs(ensemble_members_path, exist_ok=True)
    checkpoint_paths: List[str] = []

    for member_index in range(number_ensemble_members):
        pl.seed_everything(seed + member_index, workers=True)

        model = build_model_from_config(config, input_dim=input_dim, output_dim=output_dim)
        checkpoint_path = os.path.join(
            ensemble_members_path, f"ensemble_member_{member_index}.ckpt"
        )
        fit_best_model(
            model=model,
            data_module=data_module,
            max_epochs=max_epochs,
            patience=int(config["patience"]),
            monitor=config.get("early_stop_metric", "val_loss"),
            checkpoint_directory=str(ensemble_members_path),
            checkpoint_path=checkpoint_path,
        )
        checkpoint_paths.append(checkpoint_path)

    with tempfile.NamedTemporaryFile(
        mode="w", dir=ensemble_members_path, suffix=".json", delete=False
    ) as handle:
        manifest_path = Path(handle.name)
        json.dump([Path(path).name for path in checkpoint_paths], handle, indent=2)
    try:
        os.replace(manifest_path, ensemble_members_path / "ensemble_members.json")
    finally:
        manifest_path.unlink(missing_ok=True)

    return checkpoint_paths

def deploy_models(
    ensemble_members_path: str,
    file_to_score: str,
    model_column_prefix: str,
    rep_col_name: str = "Minimol_Representation",
    output_path: Optional[str] = None,
    batch_size: int = 2048,
    chunk_rows: int = 50000,
    checkpoint_paths: Optional[Sequence[str]] = None,
) -> str:
    rep_col_name = _strip(rep_col_name)
    ensemble_members_path = Path(ensemble_members_path)
    file_to_score = Path(file_to_score)
    if checkpoint_paths is None:
        manifest_path = ensemble_members_path / "ensemble_members.json"
        if manifest_path.is_file():
            with manifest_path.open() as handle:
                checkpoint_paths = [ensemble_members_path / name for name in json.load(handle)]
        else:
            checkpoint_paths = sorted(ensemble_members_path.glob("*.ckpt"))
    checkpoint_paths = [Path(path) for path in checkpoint_paths]
    if not checkpoint_paths:
        raise FileNotFoundError(f"No .ckpt files found in: {ensemble_members_path}")

    if output_path is None:
        base, ext = os.path.splitext(file_to_score)
        output_path = f"{base}_scored{ext or '.csv'}"
    output_path = Path(output_path)
    if output_path == file_to_score or (
        output_path.exists() and file_to_score.exists()
        and os.path.samefile(output_path, file_to_score)
    ):
        raise ValueError("The input CSV and output CSV must be different files.")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    models = [
        MultiTaskFeedForward.load_from_checkpoint(path, map_location=device).to(device).eval()
        for path in checkpoint_paths
    ]
    expected_feature_size = int(models[0].input_dim)
    if any(
        model.input_dim != expected_feature_size or model.output_dim != models[0].output_dim
        for model in models
    ):
        raise ValueError("Ensemble checkpoints must have matching input and output dimensions.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=output_path.parent, suffix=".csv", delete=False) as handle:
        temporary_path = Path(handle.name)
    try:
        header_written = False
        with pd.read_csv(file_to_score, sep=",", chunksize=chunk_rows) as reader:
            for chunk_dataframe in reader:
                chunk_dataframe.columns = [c.strip() for c in chunk_dataframe.columns]
                if rep_col_name not in chunk_dataframe.columns:
                    raise KeyError(
                        f"Missing representation column '{rep_col_name}'. Available: {chunk_dataframe.columns.tolist()}"
                    )

                keep_mask = []
                parsed_features = []

                for item in chunk_dataframe[rep_col_name]:
                    try:
                        array = np.asarray(
                            ast.literal_eval(item) if isinstance(item, str) else item,
                            dtype=np.float32,
                        ).ravel()
                    except (ValueError, SyntaxError, TypeError, OverflowError):
                        keep_mask.append(False)
                        continue

                    keep = array.size == expected_feature_size and np.isfinite(array).all()
                    keep_mask.append(keep)
                    if keep:
                        parsed_features.append(array)

                chunk_dataframe = chunk_dataframe.loc[keep_mask].copy()
                if chunk_dataframe.empty:
                    continue

                features = torch.from_numpy(np.stack(parsed_features, axis=0))
                dataloader = DataLoader(
                    torch.utils.data.TensorDataset(features),
                    batch_size=batch_size,
                    shuffle=False,
                    pin_memory=(device.type == "cuda"),
                )

                ensemble_count = 0
                mean_predictions = None
                mean_square_error = None

                for model_index, model in enumerate(models):
                    batch_predictions: List[np.ndarray] = []
                    with torch.inference_mode():
                        for (x_batch,) in dataloader:
                            logits = model(x_batch.to(device))
                            batch_predictions.append(torch.sigmoid(logits).cpu().numpy())

                    predictions = np.concatenate(batch_predictions, axis=0)

                    output_dim = predictions.shape[1]
                    column_names = [f"{model_column_prefix}{model_index}"] if output_dim == 1 else [
                        f"{model_column_prefix}{model_index}_task{task_index}" for task_index in range(output_dim)
                    ]
                    for task_index, column_name in enumerate(column_names):
                        chunk_dataframe[column_name] = predictions[:, task_index]

                    ensemble_count += 1
                    if mean_predictions is None:
                        mean_predictions = predictions.astype(np.float64, copy=True)
                        mean_square_error = np.zeros_like(mean_predictions)
                    else:
                        delta = predictions - mean_predictions
                        mean_predictions += delta / ensemble_count
                        mean_square_error += delta * (predictions - mean_predictions)

                variance_predictions = mean_square_error / max(ensemble_count, 1)
                std_predictions = np.sqrt(variance_predictions)

                output_dim = mean_predictions.shape[1]
                for task_index in range(output_dim):
                    suffix = "" if output_dim == 1 else f"_task{task_index}"
                    chunk_dataframe[f"{model_column_prefix}mean{suffix}"] = mean_predictions[:, task_index]
                    chunk_dataframe[f"{model_column_prefix}std{suffix}"] = std_predictions[:, task_index]

                chunk_dataframe.to_csv(
                    temporary_path,
                    sep=",",
                    index=False,
                    mode="a",
                    header=not header_written,
                )
                header_written = True
        if not header_written:
            raise ValueError("The input CSV contains no valid feature representations to score.")
        os.replace(temporary_path, output_path)
    finally:
        temporary_path.unlink(missing_ok=True)

    return str(output_path)


##############################################################################
# 7) Append Metrics Callback (hyperparameter spreadsheet / CSV)
##############################################################################

def format_pythonic_list(values: List[float]) -> str:
    return repr([float(value) for value in values])

class AppendMetricsCallback(RayTuneCallback):
    """
    Appends each trial's final result to an ongoing summary_metrics.csv file.
    Lists are written as pythonic list strings (e.g., "[0.12, 0.09, ...]").
    """

    def __init__(self, csv_path: str):
        self.spreadsheet_path = csv_path
        self.columns: Optional[List[str]] = None

    def on_trial_complete(self, iteration, trials, trial, **info):
        last_result = dict(trial.last_result) if trial.last_result else {}
        configuration = dict(trial.config) if trial.config else {}

        row_payload: Dict[str, object] = {}
        for key, value in configuration.items():
            row_payload[f"hyperparam_{key}"] = value
        for key, value in last_result.items():
            if isinstance(value, list):
                row_payload[key] = format_pythonic_list(value)
            else:
                row_payload[key] = value

        spreadsheet_path = Path(self.spreadsheet_path)
        spreadsheet_path.parent.mkdir(parents=True, exist_ok=True)
        has_rows = spreadsheet_path.is_file() and spreadsheet_path.stat().st_size > 0
        if self.columns is None:
            self.columns = (
                pd.read_csv(spreadsheet_path, nrows=0).columns.tolist()
                if has_rows else list(row_payload)
            )

        new_columns = [column for column in row_payload if column not in self.columns]
        if new_columns:
            columns = self.columns + new_columns
            if has_rows:
                existing = pd.read_csv(
                    spreadsheet_path, dtype=str, keep_default_na=False
                ).reindex(columns=columns)
                with tempfile.NamedTemporaryFile(
                    dir=spreadsheet_path.parent, suffix=".csv", delete=False
                ) as handle:
                    temporary_path = Path(handle.name)
                try:
                    existing.to_csv(temporary_path, index=False)
                    os.replace(temporary_path, spreadsheet_path)
                finally:
                    temporary_path.unlink(missing_ok=True)
            self.columns = columns

        row_df = pd.DataFrame([row_payload]).reindex(columns=self.columns)
        row_df.to_csv(spreadsheet_path, mode="a", index=False, header=not has_rows)

##############################################################################
# 8) Search space builder
##############################################################################

def build_search_space(model_dir: str) -> Dict[str, object]:
    """Build the hyperparameter search space."""
    space: Dict[str, object] = {
        "dim_size": tune.lograndint(128, 4096),
        "dim_shrinking_scale": tune.uniform(0.5, 1.0),
        "num_layers": tune.randint(1, 3),
        "learning_rate": tune.loguniform(1e-6, 1e-1),
        "batch_size": tune.choice([32, 64, 128]),
        "L2_weight_norm": tune.loguniform(1e-8, 1e-2),
        "L1_weight_norm": tune.loguniform(1e-12, 1e-2),
        "activation_function": tune.choice(["relu", "tanh", "leaky_relu"]),
        "patience": tune.choice([2]),
        "dropout_rate": tune.uniform(0.0, 0.9),
        "use_residual": tune.choice([True, False]),
        "use_batch_norm": tune.choice([True, False]),
        "early_stop_metric": tune.choice(["val_loss"]),
        "scheduler_step_size": tune.lograndint(2, 20),
        "scheduler_gamma": tune.uniform(0.1, 0.9),
        "model_file_path": tune.choice([model_dir]),
    }

    space["task1_loss_scale"] = tune.loguniform(0.2, 2.0)

    return space

##############################################################################
# 9) main()
##############################################################################

def run_hyperparameter_search(
    df_train: pd.DataFrame,
    df_val: pd.DataFrame,
    df_test: pd.DataFrame,
    rep_col_name: str,
    output_root: Path,
    num_trials: int,
    max_epochs: int,
):
    if not ray.is_initialized():
        ray.init()

    output_root.mkdir(parents=True, exist_ok=True)
    df_extra_test = None
    experiments = [
        ("GI_60", [GI_60_COL]),
        ("GI_80", [GI_80_COL]),
        ("GI_60_GI_80", [GI_60_COL, GI_80_COL]),
    ]
    best_configs = {}

    for experiment_name, target_cols in experiments:
        experiment_dir = os.path.join(output_root, experiment_name)
        os.makedirs(experiment_dir, exist_ok=True)

        # Ensure a clean metrics file per experiment run
        summary_metrics_path = os.path.join(experiment_dir, "summary_metrics.csv")
        if os.path.exists(summary_metrics_path):
            os.remove(summary_metrics_path)

        append_metrics_callback = AppendMetricsCallback(csv_path=summary_metrics_path)

        model_dir = os.path.join(experiment_dir, "models")
        os.makedirs(model_dir, exist_ok=True)

        search_space = build_search_space(model_dir=model_dir)
        algo = OptunaSearch(metric="val_auprc_avg", mode="max")

        analysis = tune.run(
            tune.with_parameters(
                train_with_tune,
                df_train=df_train,
                df_val=df_val,
                df_test=df_test,
                rep_col=rep_col_name,
                target_cols=target_cols,
                max_epochs=max_epochs,
                df_extra_test=df_extra_test,
            ),
            config=search_space,
            search_alg=algo,
            num_samples=num_trials,
            resources_per_trial={"cpu": 1, "gpu": 0},
            metric="val_auprc_avg",
            mode="max",
            storage_path=str(output_root),
            name=experiment_name,
            callbacks=[append_metrics_callback],
        )

        best_config = analysis.get_best_config(metric="val_auprc_avg", mode="max")
        print(f"[{experiment_name}] Best hyperparameters found were:", best_config)

        best_configs[experiment_name] = best_config
        best_config_payload = {
            "experiment_name": experiment_name,
            "target_cols": list(target_cols),
            "best_config": best_config,
        }
        with open(os.path.join(experiment_dir, "best_config.json"), "w") as f:
            json.dump(best_config_payload, f, indent=2)

    return best_configs


TRAINING_DATA_PATH = None
SCREENING_DATA_PATH = None
OUTPUT_ROOT = None
RUN_HYPERPARAMETER_SEARCH = False
NUM_TRIALS = 600
NUMBER_ENSEMBLE_MEMBERS = 10
MAX_EPOCHS = 25

REP_COL_NAME = "Minimol_Representation"
GI_60_COL = "Del6_Binarized_Growth_Inhibition_at_60_Percent"
GI_80_COL = "Del6_Binarized_Growth_Inhibition_at_80_Percent"

CONFIG_FOR_GI_60 = {
    "dim_size": 1490,
    "dim_shrinking_scale": 0.9588910743889216,
    "num_layers": 1,
    "learning_rate": 0.0006374894557852803,
    "batch_size": 64,
    "L2_weight_norm": 6.0810680232663006e-05,
    "L1_weight_norm": 3.456979785934819e-08,
    "activation_function": "relu",
    "patience": 2,
    "dropout_rate": 0.44778427790914066,
    "use_residual": True,
    "use_batch_norm": True,
    "early_stop_metric": "val_loss",
    "scheduler_step_size": 14,
    "scheduler_gamma": 0.8846764337685703,
    "task1_loss_scale": 0.2027645831728871,
}

CONFIG_FOR_GI_80 = {
    "dim_size": 827,
    "dim_shrinking_scale": 0.529706502311749,
    "num_layers": 1,
    "learning_rate": 0.0013255364740731724,
    "batch_size": 128,
    "L2_weight_norm": 1.7131658909380193e-07,
    "L1_weight_norm": 9.792956870780199e-09,
    "activation_function": "leaky_relu",
    "patience": 2,
    "dropout_rate": 0.5476514432568991,
    "use_residual": True,
    "use_batch_norm": True,
    "early_stop_metric": "val_loss",
    "scheduler_step_size": 12,
    "scheduler_gamma": 0.8033576541765975,
    "task1_loss_scale": 0.2895772459745544,
}

def main():
    parser = argparse.ArgumentParser(
        description="Train DeepFlux ensembles with optional hyperparameter search and screening."
    )
    parser.add_argument(
        "--training-data", type=Path, default=TRAINING_DATA_PATH,
        help="Training CSV containing representations and both binary target columns.",
    )
    parser.add_argument(
        "--screening-data", type=Path, default=SCREENING_DATA_PATH,
        help="Optional CSV to score with the trained GI60 and GI80 ensembles.",
    )
    parser.add_argument(
        "--output-root", type=Path, default=OUTPUT_ROOT,
        help="Directory for checkpoints, tuning results, and scored CSVs.",
    )
    parser.add_argument(
        "--tune", action="store_true", default=RUN_HYPERPARAMETER_SEARCH,
        help="Run hyperparameter search before training the ensembles.",
    )
    parser.add_argument("--num-trials", type=int, default=NUM_TRIALS)
    parser.add_argument("--ensemble-members", type=int, default=NUMBER_ENSEMBLE_MEMBERS)
    parser.add_argument("--max-epochs", type=int, default=MAX_EPOCHS)
    args = parser.parse_args()

    if args.training_data is None or args.output_root is None:
        parser.error("Set TRAINING_DATA_PATH and OUTPUT_ROOT, or pass --training-data and --output-root.")
    if min(args.num_trials, args.ensemble_members, args.max_epochs) < 1:
        parser.error("Trial, ensemble-member, and epoch counts must be positive.")

    training_data_path = Path(args.training_data).expanduser().resolve()
    output_root = Path(args.output_root).expanduser().resolve()
    screening_data_path = (
        Path(args.screening_data).expanduser().resolve() if args.screening_data is not None else None
    )
    if not training_data_path.is_file():
        parser.error(f"Training CSV not found: {training_data_path}")
    if screening_data_path is not None and not screening_data_path.is_file():
        parser.error(f"Screening CSV not found: {screening_data_path}")

    df = read_data(str(training_data_path), REP_COL_NAME, [GI_80_COL, GI_60_COL])
    df_train, df_val, df_test = train_val_test_split(df, (0.80, 0.10, 0.10))

    configurations = {"GI_60": CONFIG_FOR_GI_60, "GI_80": CONFIG_FOR_GI_80}
    if args.tune:
        configurations = run_hyperparameter_search(
            df_train=df_train,
            df_val=df_val,
            df_test=df_test,
            rep_col_name=REP_COL_NAME,
            output_root=output_root / "Run_Results",
            num_trials=args.num_trials,
            max_epochs=args.max_epochs,
        )

    input_file = screening_data_path
    experiments = [("GI_60", GI_60_COL, "GI60"), ("GI_80", GI_80_COL, "GI80")]
    scored_suffix = ""
    for experiment_name, target_col, ensemble_name in experiments:
        checkpoint_directory = output_root / "final_ensembles" / ensemble_name.lower()
        checkpoint_paths = train_ensemble_models(
            config=configurations[experiment_name],
            df_train=df_train,
            df_val=df_val,
            df_test=df_test,
            rep_col_name=REP_COL_NAME,
            target_col_names=[target_col],
            number_ensemble_members=args.ensemble_members,
            ensemble_members_path=str(checkpoint_directory),
            max_epochs=args.max_epochs,
        )

        if input_file is not None:
            scored_suffix += f"_{ensemble_name}"
            output_file = output_root / "Scores" / f"{screening_data_path.stem}{scored_suffix}.csv"
            input_file = Path(deploy_models(
                ensemble_members_path=str(checkpoint_directory),
                file_to_score=str(input_file),
                model_column_prefix=f"{experiment_name}_model_",
                rep_col_name=REP_COL_NAME,
                output_path=str(output_file),
                checkpoint_paths=checkpoint_paths,
            ))


if __name__ == "__main__":
    main()
