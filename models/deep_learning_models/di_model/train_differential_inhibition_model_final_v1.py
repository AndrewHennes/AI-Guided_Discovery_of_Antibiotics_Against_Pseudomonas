#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Train, tune, and deploy DeltaFlux ensembles using configured paths or CLI arguments."""

import argparse
import os
import ast
import json
import glob
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy
import pandas

from sklearn.model_selection import KFold
from sklearn.metrics import average_precision_score, roc_auc_score, r2_score

import torch
import torch.nn as nn
import torch.nn.functional as torch_functional
from torch.utils.data import Dataset, DataLoader

import pytorch_lightning as pytorch_lightning
from pytorch_lightning.callbacks import EarlyStopping, ModelCheckpoint
from pytorch_lightning import seed_everything

import ray
from ray import tune
from ray.tune.search.optuna import OptunaSearch
from ray.tune.callback import Callback as RayTuneCallback

def _write_json_atomic(payload, path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".json", delete=False) as handle:
        temporary_path = Path(handle.name)
    try:
        with temporary_path.open("w") as handle:
            json.dump(payload, handle, indent=4)
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)


def _write_csv_atomic(dataframe: pandas.DataFrame, path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".csv", delete=False) as handle:
        temporary_path = Path(handle.name)
    try:
        dataframe.to_csv(temporary_path, sep=",", index=False)
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)

##############################################################################
# 1) Read Data
##############################################################################

def read_data(
    file_path: str,
    representation_column_name: str,
    target_column_names: List[str],
    keep_all_columns: bool = False,
) -> pandas.DataFrame:
    """
    Reads a CSV and returns a dataframe containing:
      - representation_column_name
      - target_column_names (length >= 1)

    If keep_all_columns=True, returns all columns from the CSV, but still
    enforces that representation_column_name and all target columns exist and
    are non-null/numeric (targets coerced to numeric).
    """
    if not isinstance(target_column_names, list) or len(target_column_names) < 1:
        raise ValueError(
            "target_column_names must be a non-empty list of column names, "
            f"but got: {target_column_names}"
        )

    dataframe = pandas.read_csv(file_path, sep=",")
    print("Columns in the DataFrame are:", dataframe.columns.tolist())

    required_columns = list(dict.fromkeys([representation_column_name] + target_column_names))
    missing_columns = [column for column in required_columns if column not in dataframe.columns]
    if missing_columns:
        raise KeyError(f"Missing columns in input file: {missing_columns}")

    if not keep_all_columns:
        dataframe = dataframe[required_columns].copy()
    else:
        dataframe = dataframe.copy()

    # Coerce all targets to numeric
    for target_column_name in dict.fromkeys(target_column_names):
        dataframe[target_column_name] = pandas.to_numeric(dataframe[target_column_name], errors="coerce")

    # Drop rows missing representation or ANY target
    dataframe = dataframe.dropna(subset=[representation_column_name] + target_column_names).reset_index(drop=True)
    if dataframe.empty:
        raise ValueError("No rows remain after removing missing representations and targets.")
    if not numpy.isfinite(dataframe[target_column_names].to_numpy(dtype=numpy.float32)).all():
        raise ValueError("Target values must be finite numbers.")

    return dataframe

##############################################################################
# 2) Write + Load 10-Fold Cross-Validation Splits (for reproducibility)
##############################################################################

def write_kfold_splits(
    splits_output_directory_path: str,
    number_of_samples: int,
    number_of_splits: int = 10,
    random_state: int = 42,
    dataframe: Optional[pandas.DataFrame] = None,
    write_fold_csvs: bool = True,
) -> str:
    """
    Writes KFold train/test indices to disk as JSON files:

      splits_output_directory_path/
        splits_meta.json
        fold_00.json
        fold_01.json
        ...
        fold_09.json

    Additionally, if dataframe is provided and write_fold_csvs=True, writes CSVs
    containing the rows used in each fold:

      splits_output_directory_path/
        fold_00_train.csv
        fold_00_test.csv
        ...
        fold_09_train.csv
        fold_09_test.csv

    The CSVs include a 'cv_index' column that corresponds to the row index used
    by the splits (0..N-1 after cleaning/reset_index in read_data()).

    Returns the splits_output_directory_path.
    """
    kfold = KFold(n_splits=number_of_splits, shuffle=True, random_state=random_state)
    fold_indices = list(kfold.split(numpy.arange(number_of_samples)))
    os.makedirs(splits_output_directory_path, exist_ok=True)

    if dataframe is not None:
        if len(dataframe) != int(number_of_samples):
            raise ValueError(
                "Provided dataframe length does not match number_of_samples. "
                f"len(dataframe)={len(dataframe)} vs number_of_samples={number_of_samples}"
            )

    meta_path = os.path.join(splits_output_directory_path, "splits_meta.json")
    meta = {
        "number_of_samples": int(number_of_samples),
        "number_of_splits": int(number_of_splits),
        "random_state": int(random_state),
        "splitter": "KFold",
        "shuffle": True,
        "fold_csvs_written": bool(write_fold_csvs and (dataframe is not None)),
    }
    with open(meta_path, "w") as meta_file:
        json.dump(meta, meta_file, indent=4)

    dataframe_with_cv_index: Optional[pandas.DataFrame] = None
    if dataframe is not None and bool(write_fold_csvs):
        dataframe_with_cv_index = dataframe.copy()
        if "cv_index" not in dataframe_with_cv_index.columns:
            dataframe_with_cv_index.insert(0, "cv_index", numpy.arange(len(dataframe_with_cv_index), dtype=int))
        else:
            dataframe_with_cv_index["cv_index"] = numpy.arange(len(dataframe_with_cv_index), dtype=int)

    for fold_index, (training_indices, testing_indices) in enumerate(fold_indices):
        fold_path = os.path.join(splits_output_directory_path, f"fold_{fold_index:02d}.json")
        fold_payload = {
            "fold_index": int(fold_index),
            "train_indices": [int(index) for index in training_indices],
            "test_indices": [int(index) for index in testing_indices],
        }
        with open(fold_path, "w") as fold_file:
            json.dump(fold_payload, fold_file, indent=4)

        if dataframe_with_cv_index is not None:
            fold_train_csv_path = os.path.join(splits_output_directory_path, f"fold_{fold_index:02d}_train.csv")
            fold_test_csv_path = os.path.join(splits_output_directory_path, f"fold_{fold_index:02d}_test.csv")

            train_dataframe = dataframe_with_cv_index.iloc[training_indices].copy()
            test_dataframe = dataframe_with_cv_index.iloc[testing_indices].copy()

            train_dataframe.to_csv(fold_train_csv_path, sep=",", index=False)
            test_dataframe.to_csv(fold_test_csv_path, sep=",", index=False)

    print(f"Wrote {number_of_splits} cross-validation splits to: {splits_output_directory_path}")
    if dataframe_with_cv_index is not None:
        print("Also wrote per-fold train/test CSV files (rows used in each fold).")
    return splits_output_directory_path


def load_kfold_splits(
    splits_output_directory_path: str,
    number_of_splits: int = 10,
    number_of_samples: Optional[int] = None,
) -> List[Dict[str, List[int]]]:
    """
    Loads the fold JSON files written by write_kfold_splits().
    """
    with open(os.path.join(splits_output_directory_path, "splits_meta.json")) as file_handle:
        metadata = json.load(file_handle)
    if int(metadata["number_of_splits"]) != int(number_of_splits):
        raise ValueError("The saved fold count does not match number_of_splits.")
    saved_sample_count = int(metadata["number_of_samples"])
    if number_of_samples is not None and saved_sample_count != int(number_of_samples):
        raise ValueError("The saved splits do not match the current dataset size.")
    expected_indices = set(range(saved_sample_count))
    all_test_indices = []
    splits: List[Dict[str, List[int]]] = []
    for fold_index in range(number_of_splits):
        fold_path = os.path.join(splits_output_directory_path, f"fold_{fold_index:02d}.json")
        with open(fold_path, "r") as fold_file:
            fold_payload = json.load(fold_file)
        train_indices = fold_payload["train_indices"]
        test_indices = fold_payload["test_indices"]
        if (
            not train_indices or not test_indices
            or any(type(index) is not int for index in train_indices + test_indices)
            or len(set(train_indices)) != len(train_indices)
            or len(set(test_indices)) != len(test_indices)
            or set(train_indices).intersection(test_indices)
            or set(train_indices).union(test_indices) != expected_indices
        ):
            raise ValueError(f"Invalid or overlapping indices in fold {fold_index}.")
        all_test_indices.extend(test_indices)
        splits.append(
            {
                "train_indices": fold_payload["train_indices"],
                "test_indices": fold_payload["test_indices"],
            }
        )
    if sorted(all_test_indices) != list(range(saved_sample_count)):
        raise ValueError("Each sample must occur in exactly one test fold.")
    return splits

##############################################################################
# 3) Dataset + DataModule (array-backed to avoid re-parsing each fold)
##############################################################################

def _parse_representation(item) -> numpy.ndarray:
    array = numpy.asarray(ast.literal_eval(item) if isinstance(item, str) else item, dtype=numpy.float32)
    if array.ndim != 1 or array.size == 0 or not numpy.isfinite(array).all():
        raise ValueError("Representations must be nonempty, finite vectors.")
    return array


def parse_representation_column(representation_series: pandas.Series) -> numpy.ndarray:
    """
    Parses a column of list-like representations into a 2D numpy array (N, D).
    Accepts either:
      - string representations of Python lists (e.g., "[0, 1, 0, ...]")
      - actual list/array objects
    """
    parsed_rows: List[numpy.ndarray] = []
    for row_index, item in enumerate(representation_series.values):
        try:
            parsed_array = _parse_representation(item)
        except (ValueError, SyntaxError, TypeError, OverflowError) as error:
            raise ValueError(f"Invalid representation at row {row_index}.") from error
        parsed_rows.append(parsed_array)

    feature_lengths = {row.shape[0] for row in parsed_rows}
    if len(feature_lengths) != 1:
        raise ValueError(f"Inconsistent feature vector sizes found: {sorted(feature_lengths)}")

    return numpy.stack(parsed_rows, axis=0).astype(numpy.float32)


class ArrayRegressionDataset(Dataset):
    def __init__(self, features: numpy.ndarray, targets: numpy.ndarray, indices: List[int]):
        super().__init__()
        self.features_tensor = torch.from_numpy(features[indices])
        self.targets_tensor = torch.from_numpy(targets[indices])

    def __len__(self) -> int:
        return int(self.features_tensor.shape[0])

    def __getitem__(self, index: int) -> Dict[str, torch.Tensor]:
        return {"features": self.features_tensor[index], "targets": self.targets_tensor[index]}


class RegressionDataModule(pytorch_lightning.LightningDataModule):
    def __init__(
        self,
        features: numpy.ndarray,
        targets: numpy.ndarray,
        training_indices: List[int],
        validation_indices: List[int],
        testing_indices: List[int],
        batch_size: int,
    ):
        super().__init__()
        self.features = features
        self.targets = targets
        self.training_indices = training_indices
        self.validation_indices = validation_indices
        self.testing_indices = testing_indices
        self.batch_size = int(batch_size)
        if self.batch_size < 1 or len(self.training_indices) == 0:
            raise ValueError("Training requires a positive batch size and at least one training sample.")

    def setup(self, stage: Optional[str] = None) -> None:
        self.training_dataset = ArrayRegressionDataset(self.features, self.targets, self.training_indices)
        self.validation_dataset = ArrayRegressionDataset(self.features, self.targets, self.validation_indices)
        self.testing_dataset = ArrayRegressionDataset(self.features, self.targets, self.testing_indices)

    def train_dataloader(self) -> DataLoader:
        batch_size = min(self.batch_size, len(self.training_dataset))
        return DataLoader(self.training_dataset, batch_size=batch_size, shuffle=True, num_workers=0, drop_last=True)

    def val_dataloader(self) -> DataLoader:
        return DataLoader(self.validation_dataset, batch_size=self.batch_size, shuffle=False, num_workers=0)

    def test_dataloader(self) -> DataLoader:
        return DataLoader(self.testing_dataset, batch_size=self.batch_size, shuffle=False, num_workers=0)

##############################################################################
# 4) Model
##############################################################################

class MLPLayer(nn.Module):
    """
    One block of: Linear -> (BatchNorm) -> Activation -> (Dropout).
    """

    def __init__(
        self,
        input_dimension: int,
        output_dimension: int,
        activation: nn.Module,
        use_batch_norm: bool = False,
        dropout_rate: float = 0.0,
    ):
        super().__init__()
        self.linear = nn.Linear(input_dimension, output_dimension)
        self.batch_norm = nn.BatchNorm1d(output_dimension) if use_batch_norm else nn.Identity()
        self.activation = activation
        self.dropout = nn.Dropout(p=float(dropout_rate)) if float(dropout_rate) > 0 else nn.Identity()

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        output = self.linear(features)
        output = self.batch_norm(output)
        output = self.activation(output)
        output = self.dropout(output)
        return output


class SigmoidMSERegressor(pytorch_lightning.LightningModule):
    """
    Constrained regression:
      - Network produces logits
      - If output_dimension == 1: apply sigmoid(logits) => predictions in (0, 1)
      - If output_dimension  > 1: apply softmax(logits, dim=1) => predictions sum to 1 across outputs
      - Training loss is MSE(predictions, targets)
    """

    def __init__(
        self,
        input_dimension: int,
        output_dimension: int = 1,
        softmax_group_sizes: Optional[List[int]] = None,
        loss_weights: Optional[List[float]] = None,
        hidden_dimension: int = 256,
        shrinking_scale: float = 0.75,
        number_of_layers: int = 3,
        learning_rate: float = 1e-3,
        batch_size: int = 32,
        l2_weight_norm: float = 0.0,
        l1_weight_norm: float = 0.0,
        activation_function: str = "relu",
        dropout_rate: float = 0.0,
        use_residual: bool = False,
        use_batch_norm: bool = False,
        scheduler_step_size: int = 5,
        scheduler_gamma: float = 0.5,
    ):
        super().__init__()
        self.save_hyperparameters()

        self.output_dimension = int(output_dimension)
        if self.output_dimension < 1:
            raise ValueError(f"output_dimension must be >= 1, got: {self.output_dimension}")

        self.softmax_group_sizes = [int(size) for size in softmax_group_sizes] if softmax_group_sizes else None
        if self.softmax_group_sizes is not None:
            if self.output_dimension == 1:
                raise ValueError("softmax_group_sizes cannot be used when output_dimension == 1.")
            if sum(self.softmax_group_sizes) != self.output_dimension or min(self.softmax_group_sizes) < 1:
                raise ValueError(
                    "softmax_group_sizes must be positive and sum to output_dimension. "
                    f"softmax_group_sizes={self.softmax_group_sizes}, output_dimension={self.output_dimension}"
                )

        if loss_weights is not None and len(loss_weights) != self.output_dimension:
            raise ValueError(
                "loss_weights must have length output_dimension. "
                f"len(loss_weights)={len(loss_weights)}, output_dimension={self.output_dimension}"
            )
        self.register_buffer(
            "loss_weights",
            torch.tensor(loss_weights, dtype=torch.float32) if loss_weights is not None else None,
        )

        activations = {
            "relu": nn.ReLU(),
            "tanh": nn.Tanh(),
            "sigmoid": nn.Sigmoid(),
            "leaky_relu": nn.LeakyReLU(),
        }
        activation = activations.get(str(activation_function), nn.ReLU())

        layers = nn.ModuleList()
        current_input_dimension = int(input_dimension)
        current_hidden_dimension = int(hidden_dimension)
        for _ in range(int(number_of_layers)):
            layers.append(
                MLPLayer(
                    input_dimension=current_input_dimension,
                    output_dimension=current_hidden_dimension,
                    activation=activation,
                    use_batch_norm=bool(use_batch_norm),
                    dropout_rate=float(dropout_rate),
                )
            )
            current_input_dimension = current_hidden_dimension
            current_hidden_dimension = max(1, int(current_hidden_dimension * float(shrinking_scale)))

        self.layers = layers
        self.logits_head = nn.Linear(current_input_dimension, int(self.output_dimension))

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        output = features
        for layer in self.layers:
            new_output = layer(output)
            if bool(self.hparams.use_residual) and new_output.shape == output.shape:
                output = output + new_output
            else:
                output = new_output

        logits = self.logits_head(output)

        if int(self.output_dimension) == 1:
            return torch.sigmoid(logits)
        if self.softmax_group_sizes is None:
            return torch.softmax(logits, dim=1)

        logits_groups = torch.split(logits, self.softmax_group_sizes, dim=1)
        return torch.cat([torch.softmax(group, dim=1) for group in logits_groups], dim=1)

    def _mse_loss(self, predictions: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        if predictions.shape != targets.shape:
            raise ValueError(f"Prediction and target shapes must match: {predictions.shape} vs {targets.shape}.")
        squared_error = torch_functional.mse_loss(predictions, targets, reduction="none")
        if self.loss_weights is not None:
            squared_error = squared_error * self.loss_weights
        return squared_error.mean()

    def training_step(self, batch: Dict[str, torch.Tensor], batch_index: int) -> torch.Tensor:
        features_batch = batch["features"]
        targets_batch = batch["targets"]
        predictions = self(features_batch)
        loss = self._mse_loss(predictions, targets_batch)

        if float(self.hparams.l1_weight_norm) > 0:
            l1_penalty = torch.tensor(0.0, device=self.device)
            for parameter in self.parameters():
                l1_penalty = l1_penalty + parameter.abs().sum()
            loss = loss + float(self.hparams.l1_weight_norm) * l1_penalty

        if float(self.hparams.l2_weight_norm) > 0:
            l2_penalty = torch.tensor(0.0, device=self.device)
            for parameter in self.parameters():
                l2_penalty = l2_penalty + (parameter ** 2).sum()
            loss = loss + float(self.hparams.l2_weight_norm) * l2_penalty

        self.log("train_loss", loss, on_step=False, on_epoch=True, prog_bar=False)
        return loss

    def validation_step(self, batch: Dict[str, torch.Tensor], batch_index: int) -> torch.Tensor:
        features_batch = batch["features"]
        targets_batch = batch["targets"]
        predictions = self(features_batch)
        loss = self._mse_loss(predictions, targets_batch)
        self.log("val_loss", loss, on_step=False, on_epoch=True, prog_bar=True)
        return loss

    def test_step(self, batch: Dict[str, torch.Tensor], batch_index: int) -> torch.Tensor:
        features_batch = batch["features"]
        targets_batch = batch["targets"]
        predictions = self(features_batch)
        loss = self._mse_loss(predictions, targets_batch)
        self.log("test_loss", loss, on_step=False, on_epoch=True, prog_bar=True)
        return loss

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=float(self.hparams.learning_rate), weight_decay=0.0)
        scheduler = torch.optim.lr_scheduler.StepLR(
            optimizer,
            step_size=int(self.hparams.scheduler_step_size),
            gamma=float(self.hparams.scheduler_gamma),
        )
        return {
            "optimizer": optimizer,
            "lr_scheduler": {
                "scheduler": scheduler,
                "interval": "epoch",
                "frequency": 1,
            },
        }

##############################################################################
# 5) Regression Metrics
##############################################################################

def compute_regression_metrics(
    model: nn.Module,
    dataloader: DataLoader,
    device: torch.device,
) -> Dict[str, object]:
    """
    Computes regression-style metrics by flattening predictions/targets across
    all outputs (works for both single-output and multi-output).
    """
    model.eval()
    predictions_list: List[torch.Tensor] = []
    targets_list: List[torch.Tensor] = []

    with torch.no_grad():
        for batch in dataloader:
            features_batch = batch["features"].to(device)
            targets_batch = batch["targets"].to(device)
            predictions_list.append(model(features_batch).detach().cpu())
            targets_list.append(targets_batch.detach().cpu())

    if not predictions_list:
        raise ValueError("Cannot compute regression metrics for an empty dataset.")
    predictions_all = torch.cat(predictions_list, dim=0)
    targets_all = torch.cat(targets_list, dim=0)

    residuals = targets_all - predictions_all
    squared_residuals = residuals ** 2
    mean_squared_error_value = squared_residuals.mean().item()
    root_mean_squared_error_value = float(numpy.sqrt(mean_squared_error_value))
    mean_absolute_error_value = residuals.abs().mean().item()

    per_column_mean_squared_error_values = squared_residuals.mean(dim=0).tolist()
    per_column_root_mean_squared_error_values = torch.sqrt(squared_residuals.mean(dim=0)).tolist()
    per_column_mean_absolute_error_values = residuals.abs().mean(dim=0).tolist()

    r_squared_value = r2_score(targets_all.numpy().ravel(), predictions_all.numpy().ravel())
    if len(targets_all) < 2:
        per_column_r_squared_values = [float("nan")] * int(targets_all.shape[1])
    else:
        per_column_r_squared_values = numpy.atleast_1d(r2_score(
            targets_all.numpy(), predictions_all.numpy(), multioutput="raw_values"
        )).tolist()

    return {
        "mse": float(mean_squared_error_value),
        "rmse": float(root_mean_squared_error_value),
        "mae": float(mean_absolute_error_value),
        "r2": float(r_squared_value),
        "mse_per_column": [float(value) for value in per_column_mean_squared_error_values],
        "rmse_per_column": [float(value) for value in per_column_root_mean_squared_error_values],
        "mae_per_column": [float(value) for value in per_column_mean_absolute_error_values],
        "r2_per_column": [float(value) for value in per_column_r_squared_values],
    }


def format_pythonic_list(values: List[float]) -> str:
    return repr([float(value) for value in values])

##############################################################################
# 6) 10-Fold CV Training per Ray Tune Trial
##############################################################################

def train_with_tune(
    configuration: Dict,
    dataset_file_path: str,
    representation_column_name: str,
    target_column_names: List[str],
    splits_output_directory_path: str,
    number_of_splits: int,
    max_epochs: int,
    random_state: int,
    softmax_group_sizes: Optional[List[int]] = None,
    model_class=SigmoidMSERegressor,
    binary_dataset_file_path: Optional[str] = None,
    binary_target_column_name: Optional[str] = None,
    validation_fraction: float = 0.1,
):
    seed_everything(int(random_state), workers=True)

    dataframe = read_data(dataset_file_path, representation_column_name, target_column_names, keep_all_columns=False)
    features = parse_representation_column(dataframe[representation_column_name])

    targets = dataframe[target_column_names].values.astype(numpy.float32)
    if targets.ndim != 2:
        raise ValueError(f"Targets must be 2D (N, C). Got shape: {targets.shape}")

    output_dimension = int(targets.shape[1])

    loss_weights = [
        float(configuration.get(f"loss_weight_{column_index}", 1.0))
        for column_index in range(output_dimension)
    ]

    binary_features = None
    binary_targets = None
    if binary_dataset_file_path is not None:
        if binary_target_column_name is None:
            raise ValueError("binary_target_column_name is required when binary_dataset_file_path is set.")
        binary_dataframe = read_data(
            binary_dataset_file_path,
            representation_column_name,
            [binary_target_column_name],
            keep_all_columns=False,
        )
        binary_features = parse_representation_column(binary_dataframe[representation_column_name])
        binary_targets = binary_dataframe[binary_target_column_name].to_numpy()
        if not numpy.isin(binary_targets, [0, 1]).all():
            raise ValueError("Binary evaluation targets must contain only 0 and 1.")
        binary_targets = binary_targets.astype(numpy.int64)
        if int(binary_features.shape[1]) != int(features.shape[1]):
            raise ValueError(
                "Binary dataset feature dimension does not match training dataset. "
                f"{binary_features.shape[1]} vs {features.shape[1]}"
            )

    binary_targets_have_both_classes = binary_targets is not None and len(numpy.unique(binary_targets)) >= 2

    splits = load_kfold_splits(
        splits_output_directory_path, number_of_splits=number_of_splits, number_of_samples=len(dataframe)
    )

    fold_test_mse_values: List[float] = []
    fold_test_rmse_values: List[float] = []
    fold_test_mae_values: List[float] = []
    fold_test_r2_values: List[float] = []

    fold_test_mse_per_column_values = []
    fold_test_rmse_per_column_values = []
    fold_test_mae_per_column_values = []
    fold_binary_auroc_per_column_values = []
    fold_binary_auprc_per_column_values = []

    for fold_index, fold in enumerate(splits):
        outer_training_indices = fold["train_indices"]
        testing_indices = fold["test_indices"]
        inner_training_indices, inner_validation_indices = _make_train_val_split_indices(
            len(outer_training_indices), validation_fraction, int(random_state) + fold_index
        )
        training_indices = [outer_training_indices[index] for index in inner_training_indices]
        validation_indices = [outer_training_indices[index] for index in inner_validation_indices]

        data_module = RegressionDataModule(
            features=features,
            targets=targets,
            training_indices=training_indices,
            validation_indices=validation_indices,
            testing_indices=testing_indices,
            batch_size=int(configuration["batch_size"]),
        )
        data_module.setup()

        model = model_class(
            input_dimension=int(features.shape[1]),
            output_dimension=int(output_dimension),
            softmax_group_sizes=softmax_group_sizes,
            loss_weights=loss_weights,
            hidden_dimension=int(configuration["hidden_dimension"]),
            shrinking_scale=float(configuration["shrinking_scale"]),
            number_of_layers=int(configuration["number_of_layers"]),
            learning_rate=float(configuration["learning_rate"]),
            batch_size=int(configuration["batch_size"]),
            l2_weight_norm=float(configuration["l2_weight_norm"]),
            l1_weight_norm=float(configuration["l1_weight_norm"]),
            activation_function=str(configuration["activation_function"]),
            dropout_rate=float(configuration["dropout_rate"]),
            use_residual=bool(configuration["use_residual"]),
            use_batch_norm=bool(configuration["use_batch_norm"]),
            scheduler_step_size=int(configuration["scheduler_step_size"]),
            scheduler_gamma=float(configuration["scheduler_gamma"]),
        )

        from ray import train
        context = tune.get_context() if hasattr(tune, "get_context") else train.get_context()
        trial_directory_path = context.get_trial_dir()
        if not trial_directory_path:
            raise RuntimeError("train_with_tune must run inside a Ray Tune trial.")

        fold_output_directory_path = os.path.join(trial_directory_path, f"fold_{fold_index:02d}")
        os.makedirs(fold_output_directory_path, exist_ok=True)
        _write_json_atomic({
            "train_indices": training_indices,
            "validation_indices": validation_indices,
            "test_indices": testing_indices,
        }, os.path.join(fold_output_directory_path, "split_indices.json"))
        if model.hparams.use_batch_norm and min(data_module.batch_size, len(training_indices)) < 2:
            raise ValueError("Batch normalization requires at least two training samples per batch.")

        early_stopping_callback = EarlyStopping(
            monitor="val_loss",
            patience=int(configuration["patience"]),
            mode="min",
            verbose=False,
        )

        checkpoint_callback = ModelCheckpoint(
            dirpath=fold_output_directory_path,
            monitor="val_loss",
            mode="min",
            save_top_k=1,
            filename=f"best-fold-{fold_index:02d}" + "-{epoch:02d}-{val_loss:.6f}",
        )

        trainer = pytorch_lightning.Trainer(
            max_epochs=int(max_epochs),
            callbacks=[early_stopping_callback, checkpoint_callback],
            enable_checkpointing=True,
            enable_progress_bar=False,
            logger=False,
            default_root_dir=fold_output_directory_path,
            accelerator="cpu",
            devices=1,
            deterministic=True,
        )

        trainer.fit(model, data_module)

        best_checkpoint_path = checkpoint_callback.best_model_path
        if best_checkpoint_path:
            best_model = model_class.load_from_checkpoint(best_checkpoint_path, map_location="cpu")
        else:
            raise RuntimeError(f"Fold {fold_index} did not produce a checkpoint.")

        device = torch.device("cpu")
        best_model.to(device)

        fold_metrics = compute_regression_metrics(best_model, data_module.test_dataloader(), device)
        fold_test_mse_values.append(fold_metrics["mse"])
        fold_test_rmse_values.append(fold_metrics["rmse"])
        fold_test_mae_values.append(fold_metrics["mae"])
        fold_test_r2_values.append(fold_metrics["r2"])

        fold_test_mse_per_column_values.append(fold_metrics["mse_per_column"])
        fold_test_rmse_per_column_values.append(fold_metrics["rmse_per_column"])
        fold_test_mae_per_column_values.append(fold_metrics["mae_per_column"])

        if binary_features is not None:
            binary_predictions = _predict_numpy_in_batches(
                best_model, binary_features, int(configuration["batch_size"]), device
            )
            if bool(binary_targets_have_both_classes):
                binary_auroc_per_column = [
                    float(roc_auc_score(binary_targets, binary_predictions[:, column_index]))
                    for column_index in range(output_dimension)
                ]
            else:
                binary_auroc_per_column = [float("nan")] * output_dimension
            binary_auprc_per_column = [
                float(average_precision_score(binary_targets, binary_predictions[:, column_index]))
                for column_index in range(output_dimension)
            ]

            fold_binary_auroc_per_column_values.append(binary_auroc_per_column)
            fold_binary_auprc_per_column_values.append(binary_auprc_per_column)

    cv_test_mse_mean = float(numpy.mean(fold_test_mse_values))
    cv_test_mse_std = float(numpy.std(fold_test_mse_values))

    cv_test_rmse_mean = float(numpy.mean(fold_test_rmse_values))
    cv_test_mae_mean = float(numpy.mean(fold_test_mae_values))
    cv_test_r2_mean = float(numpy.mean(fold_test_r2_values))

    cv_test_mse_mean_per_column = numpy.mean(fold_test_mse_per_column_values, axis=0).tolist()
    cv_test_rmse_mean_per_column = numpy.mean(fold_test_rmse_per_column_values, axis=0).tolist()
    cv_test_mae_mean_per_column = numpy.mean(fold_test_mae_per_column_values, axis=0).tolist()

    cv_binary_auroc_mean_per_column = []
    cv_binary_auprc_mean_per_column = []
    if fold_binary_auroc_per_column_values:
        cv_binary_auroc_mean_per_column = numpy.mean(fold_binary_auroc_per_column_values, axis=0).tolist()
        cv_binary_auprc_mean_per_column = numpy.mean(fold_binary_auprc_per_column_values, axis=0).tolist()

    if "1_Minus_Delta_6_Minus_WT" in target_column_names:
        optimization_column_index = target_column_names.index("1_Minus_Delta_6_Minus_WT")
    elif "1_Minus_Delta_6_Minus_Delta_40" in target_column_names:
        optimization_column_index = target_column_names.index("1_Minus_Delta_6_Minus_Delta_40")
    else:
        optimization_column_index = 0

    cv_test_mse_mean_optimized_column = float(cv_test_mse_mean_per_column[optimization_column_index])

    from ray import train

    report = tune.report if hasattr(tune, "get_context") else train.report
    report({
        "optimization_target_column_name": str(target_column_names[optimization_column_index]),
        "optimization_target_column_index": int(optimization_column_index),
        "cv_test_mse_mean_optimized_column": cv_test_mse_mean_optimized_column,

        "cv_test_mse_mean": cv_test_mse_mean,
        "cv_test_mse_std": cv_test_mse_std,
        "cv_test_rmse_mean": cv_test_rmse_mean,
        "cv_test_mae_mean": cv_test_mae_mean,
        "cv_test_r2_mean": cv_test_r2_mean,

        "cv_test_mse_per_fold": format_pythonic_list(fold_test_mse_values),
        "cv_test_rmse_per_fold": format_pythonic_list(fold_test_rmse_values),
        "cv_test_mae_per_fold": format_pythonic_list(fold_test_mae_values),
        "cv_test_r2_per_fold": format_pythonic_list(fold_test_r2_values),

        "cv_test_mse_mean_per_column": format_pythonic_list(cv_test_mse_mean_per_column),
        "cv_test_rmse_mean_per_column": format_pythonic_list(cv_test_rmse_mean_per_column),
        "cv_test_mae_mean_per_column": format_pythonic_list(cv_test_mae_mean_per_column),
        "binary_dataset_auroc_mean_per_column": format_pythonic_list(cv_binary_auroc_mean_per_column),
        "binary_dataset_auprc_mean_per_column": format_pythonic_list(cv_binary_auprc_mean_per_column),
    })

##############################################################################
# 7) Append Metrics Callback (hyperparameter spreadsheet / TSV)
##############################################################################

class AppendMetricsCallback(RayTuneCallback):
    """
    A custom Ray Tune callback that appends each trial's final result
    to an ongoing summary_metrics.tsv (or .csv) file.
    """

    def __init__(self, spreadsheet_path: str):
        self.spreadsheet_path = spreadsheet_path
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

        path = Path(self.spreadsheet_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        has_rows = path.is_file() and path.stat().st_size > 0
        if self.columns is None:
            self.columns = pandas.read_csv(path, nrows=0).columns.tolist() if has_rows else list(row_payload)
        new_columns = [column for column in row_payload if column not in self.columns]
        if new_columns:
            columns = self.columns + new_columns
            if has_rows:
                existing = pandas.read_csv(path, dtype=str, keep_default_na=False).reindex(columns=columns)
                _write_csv_atomic(existing, path)
            self.columns = columns
        row_dataframe = pandas.DataFrame([row_payload]).reindex(columns=self.columns)
        row_dataframe.to_csv(path, sep=",", mode="a", index=False, header=not has_rows)

##############################################################################
# 8) NEW: Ensemble train() + deploy()
##############################################################################

def _ensure_directory(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def _make_train_val_split_indices(
    number_of_samples: int,
    validation_fraction: float,
    seed: int,
) -> Tuple[List[int], List[int]]:
    """
    Deterministic train/validation split using a per-member seed.
    """
    if int(number_of_samples) < 2:
        raise ValueError("At least two samples are required for a training/validation split.")
    if not (0.0 < float(validation_fraction) < 1.0):
        raise ValueError(f"validation_fraction must be in (0,1), got: {validation_fraction}")

    rng = numpy.random.RandomState(int(seed))
    permuted = rng.permutation(int(number_of_samples))

    validation_size = int(round(float(validation_fraction) * float(number_of_samples)))
    validation_size = max(1, min(int(number_of_samples) - 1, validation_size))

    validation_indices = permuted[:validation_size].tolist()
    training_indices = permuted[validation_size:].tolist()
    return training_indices, validation_indices


def _resolve_device(accelerator: str) -> torch.device:
    """
    Maps an accelerator string to a torch.device.
    """
    accel = str(accelerator).lower().strip()
    if accel == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if torch.backends.mps.is_available():
            return torch.device("mps")
    if accel in {"gpu", "cuda"}:
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if accel == "mps":
        return torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    return torch.device("cpu")


def _predict_numpy_in_batches(
    model: nn.Module,
    features_numpy: numpy.ndarray,
    batch_size: int,
    device: torch.device,
) -> numpy.ndarray:
    """
    Returns predictions as a numpy array of shape (N, C).
    """
    model.eval()
    tensor_dataset = torch.utils.data.TensorDataset(torch.from_numpy(features_numpy))
    dataloader = DataLoader(tensor_dataset, batch_size=int(batch_size), shuffle=False, num_workers=0)

    predictions: List[torch.Tensor] = []
    with torch.no_grad():
        for (features_batch,) in dataloader:
            features_batch = features_batch.to(device)
            batch_predictions = model(features_batch).detach().cpu()
            predictions.append(batch_predictions)

    if not predictions:
        return numpy.zeros((0, int(getattr(model, "output_dimension", 1))), dtype=numpy.float32)

    predictions_tensor = torch.cat(predictions, dim=0)
    if predictions_tensor.ndim == 1:
        predictions_tensor = predictions_tensor.unsqueeze(1)
    return predictions_tensor.numpy()


class EnsembleWorkflow:
    """
    Train and deploy an ensemble of Lightning checkpoints using a dictionary config.

      - train(config, dataset_path, ensemble_out_path, number_of_ensemble_members):
          trains number_of_ensemble_members models using config and writes all checkpoints to ensemble_out_path

      - deploy(checkpoints_path, dataset_path_to_score, scored_dataset_output_path):
          loads checkpoints and writes scored dataset to scored_dataset_output_path
    """

    MODEL_REGISTRY = {
        "SigmoidMSERegressor": SigmoidMSERegressor,
    }

    @staticmethod
    def train(
        config: Dict,
        dataset_path: str,
        ensemble_out_path: str,
        number_of_ensemble_members: int,
        *,
        representation_column_name: Optional[str] = None,
        target_column_names: Optional[List[str]] = None,
        softmax_group_sizes: Optional[List[int]] = None,
        max_epochs: int = 25,
        random_state: int = 42,
        validation_fraction: float = 0.1,
        accelerator: str = "cpu",
        devices: int = 1,
        model_architecture_name: str = "SigmoidMSERegressor",
    ) -> Dict[str, object]:
        """
        Trains an ensemble of models using a fixed hyperparameter configuration.

        Writes to ensemble_out_path:
          - ensemble_metadata.json
          - checkpoints_manifest.json
          - member_XX-best-...ckpt (one per member)
        """
        if not isinstance(config, dict):
            raise TypeError(f"config must be a dict, got: {type(config)}")

        if int(number_of_ensemble_members) < 1:
            raise ValueError(f"number_of_ensemble_members must be >= 1, got: {number_of_ensemble_members}")

        ensemble_out_path = _ensure_directory(str(ensemble_out_path))

        resolved_representation_column_name = (
            str(representation_column_name)
            if representation_column_name is not None
            else str(config.get("representation_column_name", "Minimol_Representation"))
        )

        resolved_target_column_names = (
            list(target_column_names)
            if target_column_names is not None
            else config.get("target_column_names", None)
        )
        if resolved_target_column_names is None or not isinstance(resolved_target_column_names, list) or len(resolved_target_column_names) < 1:
            raise ValueError(
                "target_column_names must be provided either as a keyword argument "
                "or inside config['target_column_names'] as a non-empty list."
            )

        model_class = EnsembleWorkflow.MODEL_REGISTRY.get(str(model_architecture_name))
        if model_class is None:
            raise KeyError(
                f"Unknown model_architecture_name={model_architecture_name}. "
                f"Known: {sorted(EnsembleWorkflow.MODEL_REGISTRY.keys())}"
            )

        dataframe = read_data(
            file_path=str(dataset_path),
            representation_column_name=resolved_representation_column_name,
            target_column_names=resolved_target_column_names,
            keep_all_columns=False,
        )
        features = parse_representation_column(dataframe[resolved_representation_column_name])
        targets = dataframe[resolved_target_column_names].values.astype(numpy.float32)
        if targets.ndim != 2:
            raise ValueError(f"Targets must be 2D (N, C). Got shape: {targets.shape}")

        number_of_samples = int(features.shape[0])
        input_dimension = int(features.shape[1])
        output_dimension = int(targets.shape[1])

        loss_weights = [
            float(config.get(f"loss_weight_{column_index}", 1.0))
            for column_index in range(output_dimension)
        ]

        ensemble_metadata = {
            "dataset_path": str(dataset_path),
            "ensemble_out_path": str(ensemble_out_path),
            "number_of_ensemble_members": int(number_of_ensemble_members),
            "random_state": int(random_state),
            "validation_fraction": float(validation_fraction),
            "representation_column_name": str(resolved_representation_column_name),
            "target_column_names": [str(c) for c in resolved_target_column_names],
            "input_dimension": int(input_dimension),
            "output_dimension": int(output_dimension),
            "softmax_group_sizes": [int(x) for x in softmax_group_sizes] if softmax_group_sizes is not None else None,
            "model_architecture_name": str(model_architecture_name),
            "config": dict(config),
        }
        metadata_path = os.path.join(ensemble_out_path, "ensemble_metadata.json")

        checkpoint_paths: List[str] = []
        member_records: List[Dict[str, object]] = []

        for member_index in range(int(number_of_ensemble_members)):
            member_seed = int(random_state) + int(member_index)
            seed_everything(int(member_seed), workers=True)

            training_indices, validation_indices = _make_train_val_split_indices(
                number_of_samples=number_of_samples,
                validation_fraction=float(validation_fraction),
                seed=int(member_seed),
            )

            data_module = RegressionDataModule(
                features=features,
                targets=targets,
                training_indices=training_indices,
                validation_indices=validation_indices,
                testing_indices=validation_indices,
                batch_size=int(config["batch_size"]),
            )
            data_module.setup()
            if bool(config["use_batch_norm"]) and min(int(config["batch_size"]), len(training_indices)) < 2:
                raise ValueError("Batch normalization requires at least two training samples per batch.")

            model = model_class(
                input_dimension=int(input_dimension),
                output_dimension=int(output_dimension),
                softmax_group_sizes=softmax_group_sizes,
                loss_weights=loss_weights,
                hidden_dimension=int(config["hidden_dimension"]),
                shrinking_scale=float(config["shrinking_scale"]),
                number_of_layers=int(config["number_of_layers"]),
                learning_rate=float(config["learning_rate"]),
                batch_size=int(config["batch_size"]),
                l2_weight_norm=float(config["l2_weight_norm"]),
                l1_weight_norm=float(config["l1_weight_norm"]),
                activation_function=str(config["activation_function"]),
                dropout_rate=float(config["dropout_rate"]),
                use_residual=bool(config["use_residual"]),
                use_batch_norm=bool(config["use_batch_norm"]),
                scheduler_step_size=int(config["scheduler_step_size"]),
                scheduler_gamma=float(config["scheduler_gamma"]),
            )

            early_stopping_callback = EarlyStopping(
                monitor="val_loss",
                patience=int(config.get("patience", 2)),
                mode="min",
                verbose=False,
            )

            checkpoint_callback = ModelCheckpoint(
                dirpath=str(ensemble_out_path),
                monitor="val_loss",
                mode="min",
                save_top_k=1,
                filename=f"member_{member_index:02d}-best" + "-{epoch:02d}-{val_loss:.6f}",
            )

            trainer = pytorch_lightning.Trainer(
                max_epochs=int(max_epochs),
                callbacks=[early_stopping_callback, checkpoint_callback],
                enable_checkpointing=True,
                enable_progress_bar=False,
                logger=False,
                default_root_dir=str(ensemble_out_path),
                accelerator=str(accelerator),
                devices=int(devices),
                deterministic=True,
            )

            trainer.fit(model, data_module)

            best_checkpoint_path = checkpoint_callback.best_model_path
            if not best_checkpoint_path or not os.path.exists(best_checkpoint_path):
                raise RuntimeError(
                    f"Member {member_index:02d} did not produce a checkpoint. "
                    f"best_model_path={best_checkpoint_path}"
                )

            checkpoint_paths.append(str(best_checkpoint_path))
            member_records.append(
                {
                    "member_index": int(member_index),
                    "seed": int(member_seed),
                    "best_checkpoint_path": Path(best_checkpoint_path).name,
                    "num_train": int(len(training_indices)),
                    "num_val": int(len(validation_indices)),
                }
            )

        manifest = {
            "checkpoint_paths": [Path(path).name for path in checkpoint_paths],
            "members": member_records,
        }
        manifest_path = os.path.join(ensemble_out_path, "checkpoints_manifest.json")
        _write_json_atomic(ensemble_metadata, metadata_path)
        _write_json_atomic(manifest, manifest_path)

        return {
            "ensemble_out_path": str(ensemble_out_path),
            "metadata_path": str(metadata_path),
            "manifest_path": str(manifest_path),
            "checkpoint_paths": checkpoint_paths,
        }

    @staticmethod
    def deploy(
        checkpoints_path: str,
        dataset_path_to_score: str,
        scored_dataset_output_path: str,
        *,
        representation_column_name: Optional[str] = None,
        batch_size: int = 1024,
        accelerator: str = "cpu",
    ) -> str:
        """
        Loads the manifested checkpoints, scores the input CSV, and writes the output CSV.

        Output columns added:
          - prediction_mean_<target> (or prediction_mean_<i>)
          - prediction_std_<target> (or prediction_std_<i>)

        Notes:
          - Preserves original row count. Invalid or missing representations receive NaN predictions.
        """
        checkpoints_path = str(checkpoints_path)

        metadata_path = os.path.join(checkpoints_path, "ensemble_metadata.json")
        metadata: Dict[str, object] = {}
        if os.path.exists(metadata_path):
            with open(metadata_path, "r") as file_handle:
                metadata = json.load(file_handle)

        resolved_representation_column_name = (
            str(representation_column_name)
            if representation_column_name is not None
            else str(metadata.get("representation_column_name", "Minimol_Representation"))
        )
        target_column_names = metadata.get("target_column_names", None)
        model_architecture_name = str(metadata.get("model_architecture_name", "SigmoidMSERegressor"))

        model_class = EnsembleWorkflow.MODEL_REGISTRY.get(model_architecture_name)
        if model_class is None:
            raise ValueError(f"Unknown model architecture in ensemble metadata: {model_architecture_name}")

        manifest_path = os.path.join(checkpoints_path, "checkpoints_manifest.json")
        checkpoint_files: List[str] = []
        if os.path.exists(manifest_path):
            with open(manifest_path, "r") as file_handle:
                manifest = json.load(file_handle)
            for entry in manifest.get("checkpoint_paths", []):
                path = Path(entry)
                checkpoint_files.append(str(Path(checkpoints_path) / (path.name if path.is_absolute() else path)))
        else:
            checkpoint_files = sorted(glob.glob(os.path.join(checkpoints_path, "**", "*.ckpt"), recursive=True))

        if not checkpoint_files:
            raise FileNotFoundError(f"No .ckpt files found under: {checkpoints_path}")
        for checkpoint_file in checkpoint_files:
            if not Path(checkpoint_file).is_file():
                raise FileNotFoundError(f"Manifested checkpoint not found: {checkpoint_file}")

        dataframe = pandas.read_csv(str(dataset_path_to_score), sep=",")
        if resolved_representation_column_name not in dataframe.columns:
            raise KeyError(
                f"Representation column '{resolved_representation_column_name}' not found in scoring dataset. "
                f"Columns={dataframe.columns.tolist()}"
            )

        device = _resolve_device(accelerator=str(accelerator))
        first_model = model_class.load_from_checkpoint(checkpoint_files[0], map_location=device)
        expected_input_dimension = int(first_model.hparams.input_dimension)
        valid_mask = numpy.zeros(len(dataframe), dtype=bool)
        parsed_rows = []
        for row_index, item in enumerate(dataframe[resolved_representation_column_name]):
            try:
                array = _parse_representation(item)
            except (ValueError, SyntaxError, TypeError, OverflowError):
                continue
            if array.size == expected_input_dimension:
                valid_mask[row_index] = True
                parsed_rows.append(array)
        features = (
            numpy.stack(parsed_rows) if parsed_rows
            else numpy.empty((0, expected_input_dimension), dtype=numpy.float32)
        )

        member_predictions: List[numpy.ndarray] = []
        output_dimension: Optional[int] = None

        for member_index, checkpoint_path in enumerate(checkpoint_files):
            model = first_model if member_index == 0 else model_class.load_from_checkpoint(checkpoint_path, map_location=device)
            if int(model.hparams.input_dimension) != expected_input_dimension:
                raise ValueError("Ensemble checkpoints have inconsistent input dimensions.")
            model.to(device)

            predictions = _predict_numpy_in_batches(
                model=model,
                features_numpy=features,
                batch_size=int(batch_size),
                device=device,
            )

            if predictions.ndim == 1:
                predictions = predictions.reshape(-1, 1)

            if output_dimension is None:
                output_dimension = int(predictions.shape[1])
            elif int(predictions.shape[1]) != int(output_dimension):
                raise ValueError(
                    "Inconsistent output dimensions across checkpoints. "
                    f"Expected {output_dimension}, got {predictions.shape[1]} from {checkpoint_path}"
                )

            member_predictions.append(predictions.astype(numpy.float32))

        stacked = numpy.stack(member_predictions, axis=0)  # (M, N_valid, C)
        predictions_mean = stacked.mean(axis=0)            # (N_valid, C)
        predictions_std = stacked.std(axis=0)              # (N_valid, C)

        if isinstance(target_column_names, list) and output_dimension is not None and len(target_column_names) == int(output_dimension):
            if len(set(target_column_names)) != len(target_column_names):
                target_column_names = [f"{index}_{name}" for index, name in enumerate(target_column_names)]
            mean_column_names = [f"prediction_mean_{name}" for name in target_column_names]
            std_column_names = [f"prediction_std_{name}" for name in target_column_names]
        else:
            mean_column_names = [f"DeltaFlux_Mean_{i}" for i in range(int(output_dimension or 1))]
            std_column_names = [f"DeltaFlux_Std_{i}" for i in range(int(output_dimension or 1))]

        # Initialize columns with NaN for all rows, then fill only valid rows
        for column_name in mean_column_names + std_column_names:
            dataframe[column_name] = numpy.nan

        for column_index, column_name in enumerate(mean_column_names):
            dataframe.loc[valid_mask, column_name] = predictions_mean[:, column_index]

        for column_index, column_name in enumerate(std_column_names):
            dataframe.loc[valid_mask, column_name] = predictions_std[:, column_index]

        _write_csv_atomic(dataframe, scored_dataset_output_path)
        return str(scored_dataset_output_path)

##############################################################################
# 9) main()
##############################################################################

def run_hyperparameter_search(
    dataset_file_path: str,
    representation_column_name: str,
    experiments_output_directory_path: str,
    number_of_trials: int = 1000,
    number_of_splits: int = 10,
    max_epochs: int = 25,
    random_state: int = 42,
    validation_fraction: float = 0.1,
    binary_dataset_file_path: Optional[str] = None,
    binary_target_column_name: str = "1_Minus_PA_Hit",
    target_configurations=None,
):
    best_configurations = {}
    target_column_name_candidates = [
        "1_Minus_Delta_6_Minus_WT",
        "Delta_6_Minus_Delta_3",
        "Delta_3_Minus_WT",
        "1_Minus_Delta_6_Minus_Delta_40",
        "Delta_6_Minus_Delta_3",
        "Delta_3_Minus_Delta_40",
    ]

    if target_configurations is None:
        target_configurations = [
            ("column_1", [target_column_name_candidates[0]], None),
            ("column_4", [target_column_name_candidates[3]], None),
        ]

    model_architectures = {
        "sigmoid_mse_regressor": SigmoidMSERegressor,
    }

    os.makedirs(experiments_output_directory_path, exist_ok=True)

    if not ray.is_initialized():
        ray.init(ignore_reinit_error=True)


    search_space = {
        "hidden_dimension": tune.lograndint(128, 4096),
        "shrinking_scale": tune.uniform(0.5, 1.0),
        "number_of_layers": tune.randint(1, 4),
        "learning_rate": tune.loguniform(1e-6, 1e-1),
        "batch_size": tune.choice([32, 64, 128]),
        "l2_weight_norm": tune.loguniform(1e-8, 1e-2),
        "l1_weight_norm": tune.loguniform(1e-12, 1e-2),
        "activation_function": tune.choice(["relu", "tanh", "leaky_relu"]),
        "patience": tune.choice([2]),
        "dropout_rate": tune.uniform(0.0, 0.9),
        "use_residual": tune.choice([True, False]),
        "use_batch_norm": tune.choice([True, False]),
        "scheduler_step_size": tune.lograndint(2, 20),
        "scheduler_gamma": tune.uniform(0.1, 0.9),
    }

    maximum_number_of_targets = max(
        len(target_column_names) for _, target_column_names, _ in target_configurations
    )
    for target_index in range(maximum_number_of_targets):
        search_space[f"loss_weight_{target_index}"] = tune.uniform(0.2, 2.0)

    for configuration_name, target_column_names, softmax_group_sizes in target_configurations:
        configuration_output_directory_path = os.path.join(experiments_output_directory_path, configuration_name)
        performance_metrics_directory_path = os.path.join(configuration_output_directory_path, "performance_metrics")
        checkpoints_output_directory_path = os.path.join(configuration_output_directory_path, "checkpoints")
        splits_output_directory_path = os.path.join(configuration_output_directory_path, "dataset_splits")

        for directory_path in [
            performance_metrics_directory_path,
            checkpoints_output_directory_path,
            splits_output_directory_path,
        ]:
            os.makedirs(directory_path, exist_ok=True)

        dataframe = read_data(
            dataset_file_path,
            representation_column_name,
            target_column_names,
            keep_all_columns=True,
        )

        write_kfold_splits(
            splits_output_directory_path=splits_output_directory_path,
            number_of_samples=len(dataframe),
            number_of_splits=number_of_splits,
            random_state=random_state,
            dataframe=dataframe,
            write_fold_csvs=True,
        )

        for model_architecture_name, model_class in model_architectures.items():
            algorithm = OptunaSearch(metric="cv_test_mse_mean_optimized_column", mode="min")

            summary_metrics_path = os.path.join(
                performance_metrics_directory_path,
                f"metrics_{model_architecture_name}.csv",
            )
            append_metrics_callback = AppendMetricsCallback(spreadsheet_path=summary_metrics_path)

            tuner = tune.run(
                tune.with_parameters(
                    train_with_tune,
                    dataset_file_path=dataset_file_path,
                    representation_column_name=representation_column_name,
                    target_column_names=target_column_names,
                    splits_output_directory_path=splits_output_directory_path,
                    number_of_splits=number_of_splits,
                    max_epochs=max_epochs,
                    random_state=random_state,
                    softmax_group_sizes=softmax_group_sizes,
                    model_class=model_class,
                    binary_dataset_file_path=binary_dataset_file_path,
                    binary_target_column_name=binary_target_column_name,
                    validation_fraction=validation_fraction,
                ),
                config=search_space,
                search_alg=algorithm,
                num_samples=number_of_trials,
                resources_per_trial={"cpu": 1, "gpu": 0},
                metric="cv_test_mse_mean_optimized_column",
                mode="min",
                callbacks=[append_metrics_callback],
                storage_path=os.path.join(checkpoints_output_directory_path, model_architecture_name),
            )

            best_configuration = tuner.get_best_config(metric="cv_test_mse_mean_optimized_column", mode="min")
            if best_configuration is None:
                raise RuntimeError(f"No valid tuning result for {configuration_name}.")
            best_configurations[tuple(target_column_names)] = best_configuration
            _write_json_atomic({
                "target_column_names": target_column_names,
                "softmax_group_sizes": softmax_group_sizes,
                "config": best_configuration,
            }, os.path.join(configuration_output_directory_path, "best_config.json"))
            print(
                "Best hyperparameters found were:",
                configuration_name,
                model_architecture_name,
                best_configuration,
            )

    return best_configurations


TRAINING_DATA_PATH = None
SCREENING_DATA_PATH = None
BINARY_EVALUATION_DATA_PATH = None
OUTPUT_ROOT = None
CHECKPOINTS_PATH = None
SCORED_DATASET_OUTPUT_PATH = None

REPRESENTATION_COLUMN_NAME = "Minimol_Representation"
TARGET_COLUMN_NAMES = ["1_Minus_Delta_6_Minus_WT"]
SOFTMAX_GROUP_SIZES = None
BINARY_TARGET_COLUMN_NAME = "1_Minus_PA_Hit"

ENSEMBLE_CONFIG = {
    "hidden_dimension": 170,
    "shrinking_scale": 0.5654,
    "number_of_layers": 3,
    "learning_rate": 0.00118,
    "batch_size": 64,
    "l2_weight_norm": 1.25e-07,
    "l1_weight_norm": 2.23e-07,
    "activation_function": "leaky_relu",
    "patience": 2,
    "dropout_rate": 0.24,
    "use_residual": True,
    "use_batch_norm": True,
    "scheduler_step_size": 4,
    "scheduler_gamma": 0.8168,
    "loss_weight_0": 0.95,
}


def main():
    parser = argparse.ArgumentParser(description="Train, tune, and score DeltaFlux ensembles.")
    parser.add_argument("--training-data", type=Path, default=TRAINING_DATA_PATH)
    parser.add_argument("--screening-data", type=Path, default=SCREENING_DATA_PATH)
    parser.add_argument("--binary-evaluation-data", type=Path, default=BINARY_EVALUATION_DATA_PATH)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    parser.add_argument("--checkpoints", type=Path, default=CHECKPOINTS_PATH,
                        help="Checkpoint directory; defaults to OUTPUT_ROOT/Final_Model_Checkpoints.")
    parser.add_argument("--scored-output", type=Path, default=SCORED_DATASET_OUTPUT_PATH)
    parser.add_argument("--tune", action="store_true", help="Tune hyperparameters before ensemble training.")
    parser.add_argument("--skip-training", action="store_true", help="Use existing checkpoints or run tuning only.")
    parser.add_argument("--representation-column", default=None,
                        help="Training representation column, or the saved metadata column when scoring.")
    parser.add_argument("--target-columns", nargs="+", default=TARGET_COLUMN_NAMES)
    parser.add_argument("--softmax-group-sizes", nargs="+", type=int, default=SOFTMAX_GROUP_SIZES)
    parser.add_argument("--binary-target-column", default=BINARY_TARGET_COLUMN_NAME)
    parser.add_argument("--ensemble-members", type=int, default=10)
    parser.add_argument("--num-trials", type=int, default=1000)
    parser.add_argument("--num-folds", type=int, default=10)
    parser.add_argument("--max-epochs", type=int, default=25)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--validation-fraction", type=float, default=0.1)
    parser.add_argument("--accelerator", choices=["cpu", "cuda", "gpu", "mps", "auto"], default="cpu")
    args = parser.parse_args()
    representation_column_name = args.representation_column or REPRESENTATION_COLUMN_NAME

    if args.output_root is None:
        parser.error("Set OUTPUT_ROOT or pass --output-root.")
    if (not args.skip_training or args.tune) and args.training_data is None:
        parser.error("Set TRAINING_DATA_PATH or pass --training-data.")
    if args.skip_training and not args.tune and args.screening_data is None:
        parser.error("Provide --screening-data or --tune when using --skip-training.")
    if min(args.ensemble_members, args.num_trials, args.max_epochs) < 1 or args.num_folds < 2:
        parser.error("Run counts must be positive and --num-folds must be at least 2.")
    if not 0.0 < args.validation_fraction < 1.0:
        parser.error("--validation-fraction must be between 0 and 1.")
    if args.softmax_group_sizes is not None and (
        len(args.target_columns) == 1 or min(args.softmax_group_sizes) < 1
        or sum(args.softmax_group_sizes) != len(args.target_columns)
    ):
        parser.error("Softmax groups must be positive and sum to the number of output targets, which must exceed one.")

    output_root = Path(args.output_root).expanduser().resolve()
    for name in ("training_data", "screening_data", "binary_evaluation_data"):
        value = getattr(args, name)
        if value is not None:
            path = Path(value).expanduser().resolve()
            if not path.is_file():
                parser.error(f"Input CSV not found: {path}")
            setattr(args, name, str(path))
    checkpoints_path = (
        Path(args.checkpoints).expanduser().resolve() if args.checkpoints is not None
        else output_root / "final_model_checkpoints"
    )

    seed_everything(args.seed, workers=True)
    config = dict(ENSEMBLE_CONFIG)
    if args.tune:
        target_configurations = [
            ("column_1", ["1_Minus_Delta_6_Minus_WT"], None),
            ("column_4", ["1_Minus_Delta_6_Minus_Delta_40"], None),
        ]
        if args.target_columns not in [columns for _, columns, _ in target_configurations]:
            target_configurations = [("selected_targets", args.target_columns, args.softmax_group_sizes)]
        best_configurations = run_hyperparameter_search(
            dataset_file_path=args.training_data,
            representation_column_name=representation_column_name,
            experiments_output_directory_path=str(output_root / "Run_Results"),
            number_of_trials=args.num_trials,
            number_of_splits=args.num_folds,
            max_epochs=args.max_epochs,
            random_state=args.seed,
            validation_fraction=args.validation_fraction,
            binary_dataset_file_path=args.binary_evaluation_data,
            binary_target_column_name=args.binary_target_column,
            target_configurations=target_configurations,
        )
        config = best_configurations[tuple(args.target_columns)]

    if not args.skip_training:
        EnsembleWorkflow.train(
            config=config,
            dataset_path=args.training_data,
            ensemble_out_path=str(checkpoints_path),
            number_of_ensemble_members=args.ensemble_members,
            target_column_names=args.target_columns,
            representation_column_name=representation_column_name,
            softmax_group_sizes=args.softmax_group_sizes,
            max_epochs=args.max_epochs,
            random_state=args.seed,
            validation_fraction=args.validation_fraction,
            accelerator=args.accelerator,
        )

    if args.screening_data is not None:
        scored_output = (
            Path(args.scored_output).expanduser().resolve() if args.scored_output is not None
            else output_root / "Scores" / f"{Path(args.screening_data).stem}_DeltaFlux.csv"
        )
        EnsembleWorkflow.deploy(
            checkpoints_path=str(checkpoints_path),
            dataset_path_to_score=args.screening_data,
            scored_dataset_output_path=str(scored_output),
            representation_column_name=args.representation_column,
            accelerator=args.accelerator,
        )


if __name__ == "__main__":
    main()
