"""DeepFlux ensemble inference, hit-cohort ROC analysis, and molecular properties."""

from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

import figure2_config as config

ROOT = config.METHODS_DIRECTORY
WT_COADD = "PA_Binarized_at_80_Percent_Growth_Inhibition"
D5_COADD = "PA_MexB_Binarized_at_80_Percent_Growth_Inhibition"
WT_8K = "WT_Binarized_at_80_Percent_Growth_Inhibition"
D6_8K = "Delta6_Binarized_at_80_Percent_Growth_Inhibition"
PROPERTIES = ("log10IPC", "logP", "Number_of_Aromatic_Rings")
CATEGORIES = ("Delta6_only_active", "WT_and_Delta6_active")
ROC_DEFINITION = {
    "positive_class": "Delta5_hit",
    "negative_class": "inactive_in_both_strains",
    "positive_rule": "PA_MexB_Binarized_at_80_Percent_Growth_Inhibition == 1",
    "negative_rule": "PA_MexB_Binarized_at_80_Percent_Growth_Inhibition == 0 and PA_Binarized_at_80_Percent_Growth_Inhibition == 0",
    "label_column": "ROC_Delta5_Hit_Positive_Label",
    "score_transformation": "original ensemble_mean",
    "growth_inhibition_hit_threshold_percent": 80,
    "analysis_unit": "one_per_COADD_compound_ID",
    "deduplication_key": "COADD_Compounds_ID",
    "dataset_variant": "filtered",
    "representative_selection": "first eligible source record; labels, structure, and scores must agree within ID",
    "x_axis": "False positive rate",
    "y_axis": "True positive rate",
}
ROC_COHORT_NAME = "COADD_Delta5_hits_and_inactive_ROC_cohort.csv"
ROC_INPUT_NAME = "COADD_roc_input_scores.csv.gz"


def input_path(dataset):
    if config.DATA_MODE == "saved":
        return ROOT / ("2a/saved_inputs/COADD_input.csv.gz" if dataset == "COADD" else "2d/saved_inputs/8K_input.csv.gz")
    value = config.DATASET_PATHS[dataset]
    if not value:
        raise ValueError(f"Set the absolute {dataset} dataset path in figure2_config.py.")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise ValueError(f"The {dataset} dataset path must be absolute.")
    return path


def hit_categories(frame):
    wt = pd.to_numeric(frame[WT_COADD], errors="coerce")
    d5 = pd.to_numeric(frame[D5_COADD], errors="coerce")
    return np.select(
        [wt.eq(1), wt.eq(0) & d5.eq(1), wt.eq(0) & d5.eq(0)],
        ["WT_hit", "Delta5_only_hit", "non_hit_in_both_strains"],
        default="incomplete_paired_measurement",
    )


def add_roc_columns(frame, one_per_id=False):
    frame = frame.drop(columns=[column for column in frame if column.startswith("ROC_")]).copy()
    frame["COADD_Hit_Category"] = hit_categories(frame)
    wt = pd.to_numeric(frame[WT_COADD], errors="coerce")
    d5 = pd.to_numeric(frame[D5_COADD], errors="coerce")
    positive = d5.eq(1)
    negative = d5.eq(0) & wt.eq(0)
    frame["ROC_Cohort_Category"] = np.select(
        [positive, negative],
        [ROC_DEFINITION["positive_class"], ROC_DEFINITION["negative_class"]],
        default="excluded",
    )
    frame["ROC_Eligible"] = positive | negative
    frame["ROC_Included"] = frame.ROC_Eligible
    frame[ROC_DEFINITION["label_column"]] = np.where(
        frame.ROC_Included, positive.astype(int), np.nan
    )
    for model in ("GI60", "GI80"):
        frame[f"ROC_{model}_Delta5_Hit_Score"] = pd.to_numeric(frame[f"DeepFlux_{model}_mean"], errors="raise")
    if one_per_id:
        eligible = frame.loc[frame.ROC_Eligible].copy()
        if eligible.COADD_Compounds_ID.isna().any() or eligible.COADD_Compounds_ID.astype(str).str.strip().eq("").any():
            raise ValueError("Every eligible record requires a CO-ADD compound ID.")
        groups = eligible.groupby("COADD_Compounds_ID", sort=False)
        for column in (ROC_DEFINITION["label_column"], "Standardized_SMILES", "DeepFlux_GI60_mean", "DeepFlux_GI80_mean"):
            if groups[column].nunique(dropna=False).gt(1).any():
                raise ValueError(f"Conflicting {column} values within a CO-ADD compound ID.")
        first = eligible.drop_duplicates("COADD_Compounds_ID", keep="first")
        frame["ROC_Included"] = False
        frame.loc[first.index, "ROC_Included"] = True
        frame["ROC_Eligible_Records_For_ID"] = frame.COADD_Compounds_ID.map(groups.size()).fillna(0).astype(int)
    return frame


def score_coadd():
    import torch
    from shared.figure2_models import DeepFlux

    torch.set_num_threads(config.WORKER_THREADS)
    source = input_path("COADD")
    frame = pd.read_csv(source, dtype=str, keep_default_na=False)
    cache = np.load(ROOT / "2a/saved_inputs/COADD_minimol_features.npz", allow_pickle=False)
    cached_index = {smiles: index for index, smiles in enumerate(cache["smiles"])}
    cached_features = cache["features"]
    smiles = sorted(set(frame.Standardized_SMILES) - {""})
    index = {value: position for position, value in enumerate(smiles)}
    features = np.full((len(smiles), 512), np.nan, dtype=np.float32)
    for s, position in index.items():
        if s in cached_index:
            features[position] = cached_features[cached_index[s]]
    if "Minimol_Representation" in frame:
        for s, value in frame[["Standardized_SMILES", "Minimol_Representation"]].itertuples(index=False, name=None):
            if s not in index or not value:
                continue
            vector = np.fromstring(value.strip().strip("[]"), sep=",", dtype=np.float32)
            if len(vector) == 512 and np.isfinite(vector).all():
                features[index[s]] = vector
    valid = np.isfinite(features).all(axis=1)
    positions = np.asarray([index.get(s, -1) for s in frame.Standardized_SMILES])
    row_valid = np.asarray([i >= 0 and valid[i] for i in positions])
    frame["Source_Row_Index"] = np.arange(len(frame))
    frame["Minimol_Feature_Index"] = positions
    frame["DeepFlux_Scoring_Status"] = np.where(row_valid, "scored", "missing_representation")
    x = torch.from_numpy(features[valid])
    checkpoint_report = []
    for threshold in ("GI60", "GI80"):
        member_scores = []
        for number in range(10):
            path = ROOT / "2a/checkpoints" / threshold.lower() / f"ensemble_member_{number}.ckpt"
            checkpoint = torch.load(path, map_location="cpu", weights_only=True)
            model = DeepFlux(checkpoint["hyper_parameters"])
            model.load_state_dict(checkpoint["state_dict"], strict=True)
            model.eval()
            scores = np.full(len(features), np.nan)
            with torch.inference_mode():
                predictions = [torch.sigmoid(model(batch)).numpy().ravel() for batch in x.split(config.INFERENCE_BATCH_SIZE)]
            scores[valid] = np.concatenate(predictions) if predictions else np.array([])
            member_scores.append(scores)
            values = np.full(len(frame), np.nan)
            values[row_valid] = scores[positions[row_valid]]
            frame[f"DeepFlux_{threshold}_member_{number}"] = values
            checkpoint_report.append({"model": threshold, "member": number, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "parameters": checkpoint["hyper_parameters"]})
        scores = np.stack(member_scores)
        for suffix, aggregate in (("mean", scores[:, valid].mean(axis=0)), ("std", scores[:, valid].std(axis=0, ddof=0))):
            lookup = np.full(len(features), np.nan)
            lookup[valid] = aggregate
            values = np.full(len(frame), np.nan)
            values[row_valid] = lookup[positions[row_valid]]
            frame[f"DeepFlux_{threshold}_{suffix}"] = values
    frame = add_roc_columns(frame, one_per_id=True)
    output = ROOT / "2a"
    frame.to_csv(output / "COADD_DeepFlux_GI60_GI80_scored.csv", index=False)
    frame.loc[~row_valid].to_csv(output / "COADD_records_without_minimol_representation.csv", index=False)
    feature_archive = output / "saved_inputs/COADD_minimol_features.npz"
    if not (np.array_equal(cache["smiles"], np.asarray(smiles)) and np.array_equal(cached_features, features, equal_nan=True)):
        feature_archive = output / "COADD_scoring_features.npz"
        np.savez_compressed(feature_archive, smiles=np.asarray(smiles), features=features)
    report = {
        "rows": len(frame), "unique_standardized_smiles": len(smiles), "scored_rows": int(row_valid.sum()),
        "scored_unique_structures": int(valid.sum()), "missing_representations": int((~valid).sum()),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "checkpoints": checkpoint_report,
        "feature_archive": str(feature_archive.relative_to(output)),
        "roc_definition": ROC_DEFINITION,
    }
    (output / "scoring_report.json").write_text(json.dumps(report, indent=2))
    cohort = frame.loc[frame.ROC_Included].copy()
    cohort.to_csv(output / ROC_COHORT_NAME, index=False)
    frame.loc[frame.ROC_Eligible].to_csv(output / "saved_inputs" / ROC_INPUT_NAME, index=False)
    print(f"Scored {row_valid.sum():,} of {len(frame):,} CO-ADD records.", flush=True)
    return frame


@lru_cache(maxsize=1)
def roc_tables():
    if config.DATA_MODE == "saved":
        frame = pd.read_csv(ROOT / "2a/saved_inputs" / ROC_INPUT_NAME, dtype=str, keep_default_na=False)
    else:
        frame = score_coadd()
    frame = add_roc_columns(frame, one_per_id=True)
    cohort = frame.loc[frame.ROC_Included].copy()
    y = cohort[ROC_DEFINITION["label_column"]].astype(int).to_numpy()
    cohort[ROC_DEFINITION["label_column"]] = y
    if set(y) != {0, 1}:
        raise ValueError("ROC analysis requires both Delta5 hits and inactive records.")
    curves = {}
    rows = []
    metrics = []
    for model in ("GI60", "GI80"):
        score_column = f"ROC_{model}_Delta5_Hit_Score"
        scores = cohort[score_column].to_numpy()
        if not np.isfinite(scores).all():
            raise ValueError("Each ROC cohort record must have a finite ensemble score.")
        fpr, tpr, thresholds = roc_curve(y, scores, pos_label=1, drop_intermediate=False)
        plot_fpr, plot_tpr, plot_thresholds = roc_curve(y, scores, pos_label=1, drop_intermediate=True)
        auc = float(roc_auc_score(y, scores))
        curves[model] = {"fpr": plot_fpr, "tpr": plot_tpr, "threshold": plot_thresholds, "auc": auc}
        rows.extend({"Model": model, "False_Positive_Rate": f, "True_Positive_Rate": t, "Score_Threshold": s, "Score_Column": score_column} for f, t, s in zip(fpr, tpr, thresholds))
        metrics.append({"Model": model, "AUC": auc, "Positive_Class": ROC_DEFINITION["positive_class"], "Negative_Class": ROC_DEFINITION["negative_class"], "Delta5_hit_compound_IDs": int(y.sum()), "Inactive_compound_IDs": int((1-y).sum()), "Analysis_Unit": ROC_DEFINITION["analysis_unit"], "Hit_Cutoff_Percent": 80, "Score_Column": score_column, "Score_Transformation": ROC_DEFINITION["score_transformation"]})
    cohort.to_csv(ROOT / "2a" / ROC_COHORT_NAME, index=False)
    pd.DataFrame(rows).to_csv(ROOT / "2a/ROC_curve_coordinates.csv", index=False)
    pd.DataFrame(metrics).to_csv(ROOT / "2a/ROC_metrics.csv", index=False)
    (ROOT / "2a/ROC_definition.json").write_text(json.dumps(ROC_DEFINITION, indent=2) + "\n")
    counts = frame.loc[frame.ROC_Eligible].groupby("COADD_Compounds_ID", sort=False).agg(
        Eligible_Source_Records=("COADD_Compounds_ID", "size"),
        Label=(ROC_DEFINITION["label_column"], "first"),
        Representative_Source_Row=("Source_Row_Index", "first"),
    )
    counts.to_csv(ROOT / "2a/ROC_compound_ID_record_counts.csv")
    return curves


@lru_cache(maxsize=1)
def property_tables():
    frame = pd.read_csv(input_path("8K"), dtype=str, keep_default_na=False)
    wt = pd.to_numeric(frame[WT_8K], errors="coerce")
    delta6 = pd.to_numeric(frame[D6_8K], errors="coerce")
    selected = delta6.eq(1) & wt.isin([0, 1])
    hits = frame.loc[selected].drop(columns=["Minimol_Representation"], errors="ignore").copy()
    hits["Source_Row_Index"] = np.flatnonzero(selected)
    hits["Hit_Category"] = np.where(wt[selected].eq(1), "WT_and_Delta6_active", "Delta6_only_active")
    values = []
    for smiles in hits.Standardized_SMILES:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"Unable to parse hit SMILES: {smiles}")
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            ipc = float(Descriptors.Ipc(molecule))
        if not np.isfinite(ipc) or ipc <= 0:
            raise ValueError(f"IPC must be positive and finite: {smiles}")
        values.append((ipc, np.log10(ipc), Crippen.MolLogP(molecule), rdMolDescriptors.CalcNumAromaticRings(molecule)))
    hits[["IPC", *PROPERTIES]] = np.asarray(values)
    hits["Number_of_Aromatic_Rings"] = hits.Number_of_Aromatic_Rings.astype(int)
    fits = []
    tests = []
    for prop in PROPERTIES:
        for category in CATEGORIES:
            data = hits.loc[hits.Hit_Category.eq(category), prop].to_numpy(dtype=float)
            mean, sigma = stats.norm.fit(data)
            percentiles = stats.norm.ppf([0.1, 0.25, 0.5, 0.75, 0.9], loc=mean, scale=sigma)
            fits.append({"Property": prop, "Hit_Category": category, "N": len(data), "Gaussian_Mean": mean, "Gaussian_Sigma_MLE": sigma, **dict(zip(["P10", "P25", "P50", "P75", "P90"], percentiles))})
        first = hits.loc[hits.Hit_Category.eq(CATEGORIES[0]), prop].to_numpy(dtype=float)
        second = hits.loc[hits.Hit_Category.eq(CATEGORIES[1]), prop].to_numpy(dtype=float)
        result = stats.ttest_ind(first, second, equal_var=False, alternative="two-sided")
        tests.append({"Property": prop, "Test": "Welch_t_test_two_sided", "t_statistic": result.statistic, "degrees_of_freedom": result.df, "p_value": result.pvalue})
    fit_frame = pd.DataFrame(fits)
    hits.to_csv(ROOT / "2d/8K_WT_and_Delta6_only_hit_properties.csv", index=False)
    fit_frame.to_csv(ROOT / "2d/Gaussian_fit_parameters_and_percentiles.csv", index=False)
    pd.DataFrame(tests).to_csv(ROOT / "2d/Property_comparison_statistics.csv", index=False)
    return hits, fit_frame, pd.DataFrame(tests)
