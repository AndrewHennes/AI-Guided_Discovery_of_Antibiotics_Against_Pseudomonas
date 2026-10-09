"""Standardize molecular structures and fit a joint chemical-space embedding."""

from functools import lru_cache
import gzip
import hashlib
import importlib.metadata
import json
import os
import platform
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.decomposition import PCA
from openTSNE import TSNE
from threadpoolctl import threadpool_limits

import config

DIRECTORY = config.DIRECTORY
RESULTS = DIRECTORY / "analysis"
os.environ.setdefault("NUMBA_CACHE_DIR", str(RESULTS / "numba_cache"))
STANDARDIZATION = {
    "input_column": "SMILES in both source datasets",
    "steps": ["RDKit sanitized SMILES parsing", "rdMolStandardize.Cleanup", "rdMolStandardize.FragmentParent with skipStandardize=True", "rdMolStandardize.Uncharger", "RemoveHs", "canonical isomeric SMILES serialization"],
    "tautomer_canonicalization": False,
    "stereochemistry": "Retained in the standardized isomeric SMILES",
    "convergence": "Repeat cleanup and parent preparation to a fixed point; choose the lexicographically first canonical SMILES if charge redistribution cycles",
    "deduplication": "One embedding coordinate per standardized structure in the joint dataset",
}


@lru_cache(maxsize=None)
def standardize(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Missing SMILES")
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Cannot parse SMILES: {smiles}")
    original = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    current, states = original, []
    for _ in range(12):
        if current in states:
            standardized = min(states[states.index(current):])
            return original, standardized
        states.append(current)
        molecule = rdMolStandardize.Cleanup(Chem.MolFromSmiles(current))
        molecule = rdMolStandardize.FragmentParent(molecule, skipStandardize=True)
        molecule = rdMolStandardize.Uncharger().uncharge(molecule)
        molecule = Chem.RemoveHs(molecule)
        Chem.SanitizeMol(molecule)
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
        current = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    raise ValueError("Molecular standardization did not converge")


def binary_column(values):
    values = pd.to_numeric(values, errors="raise")
    if not values.dropna().isin([0, 1]).all():
        raise ValueError("Binary hit columns must contain 0, 1, or missing values")
    return values.astype(float)


def load_source(dataset, source_manifest):
    override = config.DATASET_PATHS[dataset]
    if override:
        path = Path(override).expanduser()
        if not path.is_absolute():
            raise ValueError("External dataset paths must be absolute")
        raw = path.read_bytes()
    else:
        path = DIRECTORY / "inputs" / source_manifest[dataset]["saved_filename"]
        raw = gzip.decompress(path.read_bytes())
        if hashlib.sha256(raw).hexdigest() != source_manifest[dataset]["source_sha256"]:
            raise ValueError(f"Saved source checksum mismatch: {dataset}")
    from io import BytesIO
    frame = pd.read_csv(BytesIO(raw), float_precision="round_trip", low_memory=False)
    source_hash = hashlib.sha256(raw).hexdigest()
    frame["Source_Row_Index"] = np.arange(len(frame))
    frame["Source_Dataset"] = dataset
    original, standardized = zip(*(standardize(value) for value in frame.SMILES))
    frame["Input_Canonical_Isomeric_SMILES"] = original
    if "Standardized_SMILES" in frame:
        frame = frame.rename(columns={"Standardized_SMILES": "Source_Provided_Standardized_SMILES"})
    frame["Standardized_SMILES"] = standardized
    frame["Structure_Changed_By_Standardization"] = frame.Input_Canonical_Isomeric_SMILES.ne(frame.Standardized_SMILES)
    if dataset == "8K":
        frame["WT_Hit_80"] = binary_column(frame.WT_Binarized_at_80_Percent_Growth_Inhibition)
        frame["Delta6_Hit_80"] = binary_column(frame.Delta6_Binarized_at_80_Percent_Growth_Inhibition)
    else:
        for source, target in (("WT PGI", "WT_Hit_80"), ("∆6 PGI", "Delta6_Hit_80")):
            values = pd.to_numeric(frame[source], errors="raise")
            frame[target] = values.ge(config.GROWTH_INHIBITION_CUTOFF).astype(float).where(values.notna())
    frame.to_csv(RESULTS / f"{dataset}_standardized_records.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    summary = {
        "source_filename": Path(override).name if override else source_manifest[dataset]["filename"],
        "source_sha256": source_hash,
        "records": len(frame),
        "original_unique_canonical_structures": frame.Input_Canonical_Isomeric_SMILES.nunique(),
        "standardized_unique_structures": frame.Standardized_SMILES.nunique(),
        "records_changed_by_standardization": int(frame.Structure_Changed_By_Standardization.sum()),
        "wt_hit_records": int(frame.WT_Hit_80.eq(1).sum()),
        "delta6_hit_records": int(frame.Delta6_Hit_80.eq(1).sum()),
        "missing_wt_calls": int(frame.WT_Hit_80.isna().sum()),
        "missing_delta6_calls": int(frame.Delta6_Hit_80.isna().sum()),
    }
    print(f"Standardized {dataset}: {len(frame)} records, {summary['standardized_unique_structures']} unique structures", flush=True)
    return frame, summary


def prepare_records():
    RESULTS.mkdir(exist_ok=True)
    RDLogger.DisableLog("rdApp.*")
    source_manifest = json.loads((DIRECTORY / "inputs" / "source_manifest.json").read_text())
    frames, summaries = {}, {}
    for dataset in ("8K", "New102"):
        frames[dataset], summaries[dataset] = load_source(dataset, source_manifest)
    all_keys = sorted(set(frames["8K"].Standardized_SMILES) | set(frames["New102"].Standardized_SMILES))
    table = pd.DataFrame(index=pd.Index(all_keys, name="Standardized_SMILES"))
    for dataset, frame in frames.items():
        grouped = frame.groupby("Standardized_SMILES", sort=True)
        table[f"in_{dataset}"] = table.index.isin(grouped.groups)
        table[f"Records_{dataset}"] = grouped.size().reindex(table.index, fill_value=0)
        for strain in ("WT", "Delta6"):
            table[f"{strain}_{dataset}"] = grouped[f"{strain}_Hit_80"].max().reindex(table.index)
    for key in all_keys:
        if standardize(key)[1] != key:
            raise ValueError("SMILES standardization is not idempotent")
    new_frame = frames["New102"]
    grouped = new_frame.groupby("Standardized_SMILES", sort=True)
    table["New102_Compound_IDs_JSON"] = grouped["Batch Id"].agg(lambda values: json.dumps(values.tolist())).reindex(table.index).fillna("[]")
    table["New102_Source_Rows_JSON"] = grouped.Source_Row_Index.agg(lambda values: json.dumps(list(map(int, values)))).reindex(table.index).fillna("[]")
    table = table.reset_index()
    table.to_csv(RESULTS / "unique_standardized_molecules.csv", index=False)
    changed = pd.concat([frame.loc[frame.Structure_Changed_By_Standardization, ["Source_Dataset", "Source_Row_Index", "SMILES", "Input_Canonical_Isomeric_SMILES", "Standardized_SMILES"]] for frame in frames.values()], ignore_index=True)
    changed.to_csv(RESULTS / "standardization_changes.csv", index=False)
    summary = {
        "datasets": summaries, "standardization": STANDARDIZATION,
        "joint_unique_structures": len(table),
        "structures_shared_between_datasets": int((table.in_8K & table.in_New102).sum()),
        "standardization_idempotence_verified": True,
        "new_compounds_all_present": bool(len(new_frame) == 102 and new_frame.SMILES.notna().all()),
    }
    (RESULTS / "standardization_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return table, frames, summary


def fit_coordinates(table, frames, summary):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=config.FINGERPRINT_RADIUS, fpSize=config.FINGERPRINT_BITS, includeChirality=False)
    fingerprints = np.empty((len(table), config.FINGERPRINT_BITS), dtype=np.float32)
    for i, smiles in enumerate(table.Standardized_SMILES):
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(Chem.MolFromSmiles(smiles)), fingerprints[i])
    np.savez_compressed(RESULTS / "morgan_fingerprints.npz", fingerprints=fingerprints.astype(np.uint8), standardized_smiles=table.Standardized_SMILES.to_numpy(dtype=str))
    with threadpool_limits(limits=config.WORKER_THREADS):
        pca = PCA(n_components=config.PCA_COMPONENTS, svd_solver="randomized", random_state=config.RANDOM_SEED)
        reduced = pca.fit_transform(fingerprints).astype(np.float32)
        if not np.isfinite(reduced).all() or not np.isfinite(pca.explained_variance_ratio_).all():
            raise ValueError("Nonfinite PCA projection")
        embedding = TSNE(n_components=2, perplexity=config.TSNE_PERPLEXITY, initialization="pca", metric="euclidean", neighbors="annoy", negative_gradient_method="fft", random_state=config.RANDOM_SEED, n_jobs=config.WORKER_THREADS, early_exaggeration_iter=config.EARLY_EXAGGERATION_ITERATIONS, n_iter=config.TSNE_ITERATIONS, verbose=True).fit(reduced)
    xy = np.asarray(embedding)
    if not np.isfinite(xy).all():
        raise ValueError("Nonfinite t-SNE coordinates")
    table[["tSNE_1", "tSNE_2"]] = xy
    table.to_csv(RESULTS / "joint_tsne_coordinates.csv", index=False)
    np.savez_compressed(RESULTS / "pca_projection.npz", reduced_coordinates=reduced, components=pca.components_, mean=pca.mean_, explained_variance_ratio=pca.explained_variance_ratio_)
    coordinate_lookup = table.set_index("Standardized_SMILES")[["tSNE_1", "tSNE_2"]]
    for dataset, frame in frames.items():
        plotted = frame.join(coordinate_lookup, on="Standardized_SMILES", validate="many_to_one")
        plotted.to_csv(RESULTS / f"{dataset}_records_with_tsne_coordinates.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
        if dataset == "New102":
            plotted.to_csv(DIRECTORY / "102_compounds_standardized_with_tsne_coordinates.csv", index=False)
    settings = {
        "standardization": STANDARDIZATION,
        "fingerprints": {"type": "Morgan binary", "radius": config.FINGERPRINT_RADIUS, "bits": config.FINGERPRINT_BITS, "includeChirality": False},
        "pca": {"components": config.PCA_COMPONENTS, "svd_solver": "randomized", "random_state": config.RANDOM_SEED, "explained_variance_ratio_sum": float(pca.explained_variance_ratio_.sum())},
        "tsne": {"dimensions": 2, "perplexity": config.TSNE_PERPLEXITY, "initialization": "pca", "metric": "euclidean", "neighbors": "annoy", "negative_gradient_method": "fft", "random_state": config.RANDOM_SEED, "n_jobs": config.WORKER_THREADS, "early_exaggeration_iter": config.EARLY_EXAGGERATION_ITERATIONS, "n_iter": config.TSNE_ITERATIONS, "kl_divergence": float(embedding.kl_divergence)},
        "python_version": platform.python_version(),
        "package_versions": {package: importlib.metadata.version(package) for package in ["numpy", "pandas", "rdkit", "scikit-learn", "openTSNE", "annoy", "threadpoolctl", "matplotlib", "PyMuPDF"]},
        "input_and_standardization_summary": summary,
    }
    (RESULTS / "analysis_settings.json").write_text(json.dumps(settings, indent=2) + "\n")
    print(f"Saved {len(table)} joint coordinates and all {len(frames['New102'])} new compound records", flush=True)
    return table


def main():
    table, frames, summary = prepare_records()
    fit_coordinates(table, frames, summary)


if __name__ == "__main__":
    main()
