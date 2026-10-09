"""Prepare and jointly embed all standardized source structures."""

from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import importlib.metadata
from io import BytesIO
import json
import os
from pathlib import Path
import platform

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.decomposition import PCA
from openTSNE import TSNE
from threadpoolctl import threadpool_limits

import config
from standardization import STANDARDIZATION, initialize_worker, standardize_record

DIRECTORY = config.DIRECTORY
RESULTS = DIRECTORY / "analysis"
os.environ.setdefault("NUMBA_CACHE_DIR", str(RESULTS / "numba_cache"))
SOURCES = ("8K", "COADD", "ECBD", "New102")
ID_COLUMNS = {"8K": "Catalog ID", "COADD": "COADD_Compounds_ID", "ECBD": "Inch_Key", "New102": "Batch Id"}


def load_sources():
    source_manifest = json.loads((DIRECTORY / "inputs" / "source_manifest.json").read_text())
    frames, metadata = {}, {}
    for dataset in SOURCES:
        override = config.DATASET_PATHS[dataset]
        if override:
            path = Path(override).expanduser()
            if not path.is_absolute():
                raise ValueError("External dataset paths must be absolute")
            raw = path.read_bytes()
            filename = path.name
        else:
            record = source_manifest[dataset]
            raw = gzip.decompress((DIRECTORY / "inputs" / record["saved_filename"]).read_bytes())
            if hashlib.sha256(raw).hexdigest() != record["source_sha256"]:
                raise ValueError(f"Saved input checksum differs: {dataset}")
            filename = record["filename"]
        columns = {"SMILES", ID_COLUMNS[dataset], "Common_Compound_Name", "Short Name"}
        if dataset == "New102":
            frame = pd.read_csv(BytesIO(raw), float_precision="round_trip", low_memory=False)
        else:
            frame = pd.read_csv(BytesIO(raw), usecols=lambda column: column in columns, float_precision="round_trip", low_memory=False)
        frame["Source_Row_Index"] = np.arange(len(frame))
        frame["Source_Dataset"] = dataset
        frames[dataset] = frame
        metadata[dataset] = {"source_filename": filename, "source_sha256": hashlib.sha256(raw).hexdigest(), "records": len(frame)}
        print(f"Loaded {dataset}: {len(frame):,} records", flush=True)
    return frames, metadata


def prepare_records():
    RESULTS.mkdir(exist_ok=True)
    RDLogger.DisableLog("rdApp.*")
    frames, summaries = load_sources()
    input_smiles = sorted({value for frame in frames.values() for value in frame.SMILES if isinstance(value, str) and value.strip()})
    cache_path = RESULTS / "standardization_cache.csv.gz"
    procedure_hash = hashlib.sha256((DIRECTORY / "standardization.py").read_bytes()).hexdigest()
    source_hashes = {key: value["source_sha256"] for key, value in summaries.items()}
    cache_key = {"procedure_sha256": procedure_hash, "sources": source_hashes, "rdkit_version": importlib.metadata.version("rdkit")}
    cache_metadata_path = RESULTS / "standardization_cache_provenance.json"
    if cache_path.is_file() and cache_metadata_path.is_file() and json.loads(cache_metadata_path.read_text()) == cache_key:
        cache = pd.read_csv(cache_path, keep_default_na=False)
        print(f"Loaded matching standardization cache for {len(cache):,} input SMILES", flush=True)
    else:
        prepared = []
        with ProcessPoolExecutor(max_workers=config.WORKER_THREADS, initializer=initialize_worker) as pool:
            for index, result in enumerate(pool.map(standardize_record, input_smiles, chunksize=128), start=1):
                prepared.append(result)
                if index % 10000 == 0 or index == len(input_smiles):
                    print(f"Standardized {index:,} / {len(input_smiles):,} input SMILES", flush=True)
        cache = pd.DataFrame(prepared, columns=["Input_Canonical_Isomeric_SMILES", "Standardized_SMILES", "Standardization_Status"])
        cache.insert(0, "SMILES", input_smiles)
        cache.to_csv(cache_path, index=False, compression={"method": "gzip", "mtime": 0})
        cache_metadata_path.write_text(json.dumps(cache_key, indent=2) + "\n")
    cache = cache.replace({"Input_Canonical_Isomeric_SMILES": {"": None}, "Standardized_SMILES": {"": None}})
    if len(cache) != len(input_smiles) or set(cache.SMILES) != set(input_smiles):
        raise ValueError("Standardization cache does not match all source SMILES")
    mapping = cache.set_index("SMILES")
    exclusions, changed = [], []
    for dataset, original in frames.items():
        frame = original.join(mapping, on="SMILES", validate="many_to_one")
        frame["Standardization_Status"] = frame.Standardization_Status.fillna("missing_smiles")
        frame["Structure_Changed_By_Standardization"] = frame.Input_Canonical_Isomeric_SMILES.ne(frame.Standardized_SMILES) & frame.Standardized_SMILES.notna()
        if dataset == "New102":
            for source, target in (("WT PGI", "WT_Hit_80"), ("∆6 PGI", "Delta6_Hit_80")):
                values = pd.to_numeric(frame[source], errors="raise")
                frame[target] = values.ge(config.GROWTH_INHIBITION_CUTOFF).astype(float).where(values.notna())
            if len(frame) != 102 or frame.Standardized_SMILES.isna().any():
                raise ValueError("All 102 new compounds must have valid standardized structures")
        frame.to_csv(RESULTS / f"{dataset}_standardized_records.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
        frames[dataset] = frame
        valid = frame.Standardized_SMILES.notna()
        summaries[dataset].update(valid_records=int(valid.sum()), excluded_records=int((~valid).sum()), standardized_unique_structures=frame.Standardized_SMILES.nunique(), records_changed_by_standardization=int(frame.Structure_Changed_By_Standardization.sum()))
        if dataset == "New102":
            summaries[dataset].update(wt_hit_records=int(frame.WT_Hit_80.eq(1).sum()), delta6_hit_records=int(frame.Delta6_Hit_80.eq(1).sum()), missing_wt_calls=int(frame.WT_Hit_80.isna().sum()))
        wanted = ["Source_Dataset", "Source_Row_Index", "SMILES", "Input_Canonical_Isomeric_SMILES", "Standardized_SMILES", "Standardization_Status"]
        exclusions.append(frame.loc[~valid, wanted])
        changed.append(frame.loc[frame.Structure_Changed_By_Standardization, wanted])
        print(f"{dataset}: {int(valid.sum()):,} valid records; {summaries[dataset]['standardized_unique_structures']:,} standardized structures; {int((~valid).sum())} excluded", flush=True)
    all_keys = sorted({value for frame in frames.values() for value in frame.Standardized_SMILES.dropna()})
    table = pd.DataFrame(index=pd.Index(all_keys, name="Standardized_SMILES"))
    for dataset, frame in frames.items():
        grouped = frame.loc[frame.Standardized_SMILES.notna()].groupby("Standardized_SMILES", sort=True)
        table[f"in_{dataset}"] = table.index.isin(grouped.groups)
        table[f"Records_{dataset}"] = grouped.size().reindex(table.index, fill_value=0)
        if dataset == "New102":
            for strain in ("WT", "Delta6"):
                table[f"{strain}_New102"] = grouped[f"{strain}_Hit_80"].max().reindex(table.index)
            table["New102_Compound_IDs_JSON"] = grouped["Batch Id"].agg(lambda values: json.dumps(values.tolist())).reindex(table.index).fillna("[]")
            table["New102_Source_Rows_JSON"] = grouped.Source_Row_Index.agg(lambda values: json.dumps(list(map(int, values)))).reindex(table.index).fillna("[]")
    table = table.reset_index()
    table.to_csv(RESULTS / "unique_standardized_molecules.csv", index=False)
    pd.concat(exclusions, ignore_index=True).to_csv(RESULTS / "excluded_source_records.csv", index=False)
    pd.concat(changed, ignore_index=True).to_csv(RESULTS / "standardization_changes.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    membership = table[[f"in_{dataset}" for dataset in SOURCES]].astype(int)
    (membership.T @ membership).to_csv(RESULTS / "dataset_structure_overlap.csv")
    background = table[[f"in_{dataset}" for dataset in SOURCES if dataset != "New102"]].any(axis=1)
    summary = {"datasets": summaries, "standardization": STANDARDIZATION, "joint_unique_structures": len(table), "gray_background_unique_structures": int(background.sum()), "new_structures_shared_with_background": int((table.in_New102 & background).sum()), "new_compounds_all_present": True, "background_activity_filter_applied": False}
    (RESULTS / "standardization_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return table, frames, summary


def fit_coordinates(table, frames, summary):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=config.FINGERPRINT_RADIUS, fpSize=config.FINGERPRINT_BITS, includeChirality=False)
    fingerprints = np.empty((len(table), config.FINGERPRINT_BITS), dtype=np.float32)
    for index, smiles in enumerate(table.Standardized_SMILES):
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(Chem.MolFromSmiles(smiles)), fingerprints[index])
        if (index + 1) % 50000 == 0:
            print(f"Fingerprints: {index+1:,} / {len(table):,}", flush=True)
    np.savez_compressed(RESULTS / "morgan_fingerprints.npz", fingerprints=fingerprints.astype(np.uint8))
    print(f"Fitting PCA for {len(table):,} structures", flush=True)
    with threadpool_limits(limits=config.WORKER_THREADS):
        pca = PCA(n_components=config.PCA_COMPONENTS, svd_solver="randomized", random_state=config.RANDOM_SEED)
        reduced = pca.fit_transform(fingerprints).astype(np.float32)
        if not np.isfinite(reduced).all() or not np.isfinite(pca.explained_variance_ratio_).all():
            raise ValueError("Nonfinite PCA projection")
        del fingerprints
        np.savez_compressed(RESULTS / "pca_projection.npz", reduced_coordinates=reduced, components=pca.components_, mean=pca.mean_, explained_variance_ratio=pca.explained_variance_ratio_)
        embedding = TSNE(n_components=2, perplexity=config.TSNE_PERPLEXITY, initialization="pca", metric="euclidean", neighbors="annoy", negative_gradient_method="fft", random_state=config.RANDOM_SEED, n_jobs=config.WORKER_THREADS, early_exaggeration_iter=config.EARLY_EXAGGERATION_ITERATIONS, n_iter=config.TSNE_ITERATIONS, verbose=True).fit(reduced)
    xy = np.asarray(embedding)
    if not np.isfinite(xy).all():
        raise ValueError("Nonfinite t-SNE coordinates")
    table[["tSNE_1", "tSNE_2"]] = xy
    table.to_csv(RESULTS / "joint_tsne_coordinates.csv", index=False)
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
    print(f"Saved {len(table):,} joint coordinates and all 102 new compound records", flush=True)


def main():
    table, frames, summary = prepare_records()
    fit_coordinates(table, frames, summary)


if __name__ == "__main__":
    main()
