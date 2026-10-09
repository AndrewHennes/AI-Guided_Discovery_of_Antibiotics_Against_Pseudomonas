"""Binary hit summaries, molecular representations, and chemical-space coordinates."""

from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import config

HIT_COLUMNS = {
    "8K": {
        "WT_8K": "WT_Binarized_at_80_Percent_Growth_Inhibition",
        "D6_8K": "Delta6_Binarized_at_80_Percent_Growth_Inhibition",
    },
    "39K": {"WT_39K": "PA_Hit_Binarized_at_80_Percent_Growth_Inhibition"},
    "ECBD": {"WT_ECBD": "Pseudomonas_aeruginosa_Binarized_Hit"},
    "COADD": {
        "WT_COADD": "PA_Binarized_at_80_Percent_Growth_Inhibition",
        "D5_COADD": "PA_MexB_Binarized_at_80_Percent_Growth_Inhibition",
    },
}
COADD_CONTINUOUS = {
    "WT_COADD": "PA_Growth_Inhibition",
    "D5_COADD": "PA_MexB_Growth_Inhibition",
}
PANEL_SOURCES = {"f": ("8K", "39K", "ECBD"), "g": ("8K", "COADD")}


def dataset_path(source):
    value = config.DATASET_PATHS[source]
    if not value:
        raise ValueError(f"Set the absolute {source} dataset path in config.py.")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise ValueError(f"The {source} dataset path must be absolute.")
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def binary_values(series, name):
    values = pd.to_numeric(series, errors="raise")
    if not values.dropna().isin([0, 1]).all():
        raise ValueError(f"{name} must contain 0, 1, or missing values.")
    return values.astype(float)


@lru_cache(maxsize=1)
def strain_hits():
    if config.DATA_MODE == "saved":
        frame = pd.read_csv(config.SAVED_DATA_DIRECTORY / "8K_binary_hits.csv.gz")
    elif config.DATA_MODE == "datasets":
        columns = [f"{strain}_Binarized_at_80_Percent_Growth_Inhibition" for strain in config.STRAINS]
        frame = pd.read_csv(dataset_path("8K"), usecols=columns)[columns]
        frame.columns = config.STRAINS
    else:
        raise ValueError('DATA_MODE must be "saved" or "datasets".')
    for strain in config.STRAINS:
        frame[strain] = binary_values(frame[strain], strain)
    return frame


def hit_summary():
    frame = strain_hits()
    counts = frame.sum().astype(int).to_dict()
    complete = frame.dropna().astype(int)
    intersections = {
        tuple(int(value) for value in pattern): int(count)
        for pattern, count in complete.value_counts(sort=False).items()
    }
    return counts, intersections


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@lru_cache(maxsize=1)
def prepare_molecules():
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    config.CACHE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    hashes = {source: file_hash(dataset_path(source)) for source in HIT_COLUMNS}
    records = []
    for source, mapping in HIT_COLUMNS.items():
        allowed = {"SMILES", "Standardized_SMILES", *mapping.values()}
        if source == "COADD":
            allowed.update(COADD_CONTINUOUS.values())
        frame = pd.read_csv(dataset_path(source), usecols=lambda name: name in allowed)
        smiles_column = "Standardized_SMILES" if "Standardized_SMILES" in frame else "SMILES"
        table = pd.DataFrame(index=frame.index)
        measured = np.zeros(len(frame), dtype=bool)
        for target, column in mapping.items():
            if column in frame:
                values = binary_values(frame[column], column)
            elif source == "COADD":
                values = pd.to_numeric(frame[COADD_CONTINUOUS[target]], errors="raise")
                values = values.ge(80).astype(float).where(values.notna())
            else:
                raise KeyError(f"Required column: {column}")
            table[target] = values
            measured |= values.notna().to_numpy()
        canonical = []
        for smiles in frame[smiles_column]:
            molecule = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
            canonical.append(Chem.MolToSmiles(molecule, isomericSmiles=True) if molecule else None)
        table["canonical_smiles"] = canonical
        table["source"] = source
        records.append(table.loc[table.canonical_smiles.notna() & measured])
        print(f"Prepared {source}", flush=True)
    merged = pd.concat(records, ignore_index=True)
    hit_columns = [name for mapping in HIT_COLUMNS.values() for name in mapping]
    molecules = merged.groupby("canonical_smiles", sort=True)[hit_columns].max()
    membership = pd.crosstab(merged.canonical_smiles, merged.source).gt(0)
    molecules = molecules.join(membership.add_prefix("in_"))
    (config.CACHE_DIRECTORY / "input_sha256.json").write_text(json.dumps(hashes, indent=2))
    return molecules


@lru_cache(maxsize=1)
def recompute_coordinates():
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from sklearn.decomposition import PCA
    from openTSNE import TSNE
    from threadpoolctl import threadpool_limits

    molecules = prepare_molecules()
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=config.FINGERPRINT_RADIUS, fpSize=config.FINGERPRINT_BITS
    )
    fingerprints = np.empty((len(molecules), config.FINGERPRINT_BITS), dtype=np.float32)
    for index, smiles in enumerate(molecules.index):
        molecule = Chem.MolFromSmiles(smiles)
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(molecule), fingerprints[index])
    with threadpool_limits(limits=config.WORKER_THREADS):
        pca = PCA(n_components=config.PCA_COMPONENTS, svd_solver="randomized", random_state=config.RANDOM_SEED)
        reduced = pca.fit_transform(fingerprints).astype(np.float32)
    del fingerprints
    result = {}
    for panel, sources in PANEL_SOURCES.items():
        selected = molecules[[f"in_{source}" for source in sources]].fillna(False).any(axis=1)
        with threadpool_limits(limits=config.WORKER_THREADS):
            embedding = TSNE(
                n_components=2,
                perplexity=config.TSNE_PERPLEXITY,
                initialization="pca",
                metric="euclidean",
                neighbors="annoy",
                negative_gradient_method="fft",
                random_state=config.RANDOM_SEED,
                n_jobs=config.WORKER_THREADS,
                early_exaggeration_iter=250,
                n_iter=500,
                verbose=True,
            ).fit(reduced[selected.to_numpy()])
        table = molecules.loc[selected].copy()
        table[["tSNE_1", "tSNE_2"]] = np.asarray(embedding)
        table.to_csv(config.CACHE_DIRECTORY / f"panel_{panel}_coordinates.csv")
        result[panel] = table.reset_index()
    return result


@lru_cache(maxsize=2)
def panel_coordinates(panel):
    saved = pd.read_csv(config.SAVED_DATA_DIRECTORY / f"panel_{panel}_coordinates.csv.gz")
    if config.DATA_MODE == "saved":
        return saved
    if config.DATA_MODE != "datasets":
        raise ValueError('DATA_MODE must be "saved" or "datasets".')
    if config.RECOMPUTE_EMBEDDINGS:
        return recompute_coordinates()[panel]
    molecules = prepare_molecules()
    selected = molecules[[f"in_{source}" for source in PANEL_SOURCES[panel]]].fillna(False).any(axis=1)
    table = molecules.loc[selected].copy()
    if set(table.index) != set(saved.canonical_smiles):
        raise ValueError("Set RECOMPUTE_EMBEDDINGS = True to calculate coordinates for these molecular structures.")
    coordinates = saved.set_index("canonical_smiles")[["tSNE_1", "tSNE_2"]]
    return table.join(coordinates).reset_index()


def write_summaries():
    config.OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    counts, intersections = hit_summary()
    summary = {
        "strain_counts": counts,
        "intersection_counts": [
            dict(zip(config.STRAINS, pattern), count=count)
            for pattern, count in sorted(intersections.items())
        ],
    }
    (config.OUTPUT_DIRECTORY / "hit_counts.json").write_text(json.dumps(summary, indent=2))
