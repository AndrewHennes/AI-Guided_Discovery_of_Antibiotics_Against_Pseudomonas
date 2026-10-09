"""Count and plot unique WT and mutant Pseudomonas hits by dataset."""

from functools import lru_cache
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

DIRECTORY = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision/dataset_hit_counts')
DATASET_PATHS = {'8K': '/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision/methods/saved_inputs/8K.csv.gz', '39K': '/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision/methods/saved_inputs/39K.csv.gz', 'ECBD': '/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision/methods/saved_inputs/ECBD.csv.gz', 'COADD': '/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision/methods/saved_inputs/COADD.csv.gz', 'PurchasedBroad': ''}
DATASET_LABELS = {"8K": "Broad 8K", "39K": "Broad 39K", "ECBD": "ECBD", "COADD": "CO-ADD", "PurchasedBroad": "Broad purchases"}
STRAIN_COLUMNS = {
    "8K": {
        "WT": "WT_Binarized_at_80_Percent_Growth_Inhibition",
        "Delta40": "Delta40_Binarized_at_80_Percent_Growth_Inhibition",
        "Delta3": "Delta3_Binarized_at_80_Percent_Growth_Inhibition",
        "Delta6": "Delta6_Binarized_at_80_Percent_Growth_Inhibition",
    },
    "39K": {"WT": "PA_Hit_Binarized_at_80_Percent_Growth_Inhibition"},
    "ECBD": {"WT": "Pseudomonas_aeruginosa_Binarized_Hit"},
    "COADD": {"WT": "PA_Binarized_at_80_Percent_Growth_Inhibition", "Delta5": "PA_MexB_Binarized_at_80_Percent_Growth_Inhibition"},
    "PurchasedBroad": {"WT": "WT_Hit", "Delta6": "Delta6_Hit"},
}
ID_COLUMNS = {"8K": "Catalog ID", "39K": None, "ECBD": "Inch_Key", "COADD": "COADD_Compounds_ID", "PurchasedBroad": "Name"}
STRAIN_ORDER = ["WT", "Delta40", "Delta3", "Delta6", "Delta5"]
STRAIN_LABELS = {"WT": "WT", "Delta40": "Δ40", "Delta3": "Δ3", "Delta6": "Δ6", "Delta5": "Δ5"}
STRAIN_COLORS = {"WT": "#00a650", "Delta40": "#faa21a", "Delta3": "#c03a2b", "Delta6": "#a860a5", "Delta5": "#5cc1df"}
FIGURE_STEM = "unique_wt_and_mutant_pseudomonas_hit_counts_by_dataset"


@lru_cache(maxsize=None)
def canonical_structure(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    molecule = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True) if molecule is not None else None


def load_records(dataset):
    path_text = DATASET_PATHS[dataset]
    if path_text:
        path = Path(path_text).expanduser()
        if not path.is_absolute():
            raise ValueError("Dataset paths must be absolute.")
        if dataset == "PurchasedBroad":
            frame = pd.read_excel(path, sheet_name="Data")
        else:
            wanted = {"SMILES", "Standardized_SMILES", ID_COLUMNS[dataset], *STRAIN_COLUMNS[dataset].values()}
            frame = pd.read_csv(path, usecols=lambda column: column in wanted, float_precision="round_trip")
    else:
        frame = pd.read_csv(DIRECTORY / "saved_inputs" / f"{dataset}_records.csv.gz", float_precision="round_trip")
    if "Source_Row_Index" not in frame:
        frame["Source_Row_Index"] = np.arange(len(frame))
    if dataset == "PurchasedBroad":
        for source, target in (("Percent_Growth_Inhibition_OD600_WT", "WT_Hit"), ("Percent_Growth_Inhibition_OD600_∆6_Strain", "Delta6_Hit")):
            values = pd.to_numeric(frame[source], errors="raise")
            frame[target] = values.gt(80).astype(float).where(values.notna())
    for column in STRAIN_COLUMNS[dataset].values():
        frame[column] = pd.to_numeric(frame[column], errors="raise")
        if not frame[column].dropna().isin([0, 1]).all():
            raise ValueError(f"{dataset}: {column} must contain binary calls or missing values.")
    smiles_column = "Standardized_SMILES" if "Standardized_SMILES" in frame else "SMILES"
    frame["Canonical_SMILES"] = frame[smiles_column].map(canonical_structure)
    return frame


def count_hits():
    RDLogger.DisableLog("rdApp.*")
    counts, compound_rows, exclusions = [], [], []
    wt_sets = {}
    for dataset in DATASET_LABELS:
        frame = load_records(dataset)
        valid = frame.loc[frame.Canonical_SMILES.notna()]
        groups = valid.groupby("Canonical_SMILES", sort=True)
        id_column = ID_COLUMNS[dataset]
        exclusions.append({"Dataset": dataset, "Source_Records": len(frame), "Invalid_or_Missing_Structure_Records": int(frame.Canonical_SMILES.isna().sum())})
        for strain, column in STRAIN_COLUMNS[dataset].items():
            hit_keys = groups[column].max().eq(1)
            keys = hit_keys.index[hit_keys]
            hit_records = frame.loc[frame[column].eq(1)]
            hit_ids = hit_records[id_column].dropna().astype(str).str.strip().nunique() if id_column else None
            counts.append({"Dataset": dataset, "Dataset_Label": DATASET_LABELS[dataset], "Strain": strain,
                           "Unique_Hit_Structures": len(keys), "Hit_Records": len(hit_records),
                           "Unique_Hit_Compound_IDs": hit_ids, "Compound_ID_Column": id_column,
                           "Measured_Unique_Structures": int(groups[column].count().gt(0).sum()),
                           "Structures_With_Mixed_Binary_Calls": int(groups[column].nunique().gt(1).sum()),
                           "Source_Hit_Column": column})
            for key in keys:
                rows = groups.get_group(key)
                active = rows.loc[rows[column].eq(1)]
                ids = sorted(active[id_column].dropna().astype(str).str.strip().unique()) if id_column else []
                compound_rows.append({"Dataset": dataset, "Strain": strain, "Canonical_SMILES": key,
                                      "Hit_Record_Count": len(active), "Source_Compound_IDs": " | ".join(ids),
                                      "Source_Hit_Record_Indices": " | ".join(map(str, active.Source_Row_Index.tolist()))})
            if strain == "WT":
                wt_sets[dataset] = set(keys)
        print(f"Counted {DATASET_LABELS[dataset]}", flush=True)
    table = pd.DataFrame(counts)
    table["Unique_Hit_Compound_IDs"] = table.Unique_Hit_Compound_IDs.astype("Int64")
    table.to_csv(DIRECTORY / "unique_hit_counts_by_dataset_and_strain.csv", index=False)
    pd.DataFrame(compound_rows).to_csv(DIRECTORY / "unique_hit_compounds.csv", index=False)
    pd.DataFrame(exclusions).to_csv(DIRECTORY / "input_record_summary.csv", index=False)
    wide = table.pivot(index="Dataset_Label", columns="Strain", values="Unique_Hit_Structures")
    wide = wide.reindex(index=list(DATASET_LABELS.values()), columns=STRAIN_ORDER).astype("Int64")
    wide.rename(columns=STRAIN_LABELS).to_csv(DIRECTORY / "WT_and_mutant_hit_counts.csv")
    overlap = pd.DataFrame({first: {second: len(wt_sets[first] & wt_sets[second]) for second in wt_sets} for first in wt_sets})
    overlap.to_csv(DIRECTORY / "WT_hit_overlap_between_datasets.csv")
    union = {"Counting_Unit": "canonical isomeric structure within each dataset", "Unique_WT_Hits_Across_All_Five_Datasets": len(set.union(*wt_sets.values())),
             "Unique_WT_Hits_Across_Four_Filtered_Datasets": len(set.union(*(value for key, value in wt_sets.items() if key != "PurchasedBroad")))}
    (DIRECTORY / "counting_definition.json").write_text(json.dumps(union, indent=2) + "\n")
    return table


def plot_counts(table):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()).resolve() / "pseudomonas_hit_counts_matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.linewidth": 0.8,
                         "pdf.fonttype": 42, "svg.fonttype": "path", "svg.hashsalt": "unique_pseudomonas_hit_counts"})
    figure, axes = plt.subplots(figsize=(7.3, 3.6))
    figure.subplots_adjust(left=0.11, right=0.985, bottom=0.23, top=0.82)
    width, step = 0.16, 0.19
    marker_positions = []
    for center, dataset in enumerate(DATASET_LABELS):
        selected = table.loc[table.Dataset.eq(dataset)].set_index("Strain")
        strains = [strain for strain in STRAIN_ORDER if strain in selected.index]
        positions = center + (np.arange(len(strains)) - (len(strains) - 1) / 2) * step
        for x, strain in zip(positions, strains):
            count = int(selected.loc[strain, "Unique_Hit_Structures"])
            axes.bar(x, count, width=width, color=STRAIN_COLORS[strain], edgecolor="black", linewidth=0.75, zorder=3)
            axes.text(x, count + 2.2, str(count), ha="center", va="bottom", fontsize=10)
            marker_positions.append({"Dataset": dataset, "Strain": strain, "X": float(x), "Count": count})
    axes.set_xlim(-0.58, len(DATASET_LABELS) - 0.45)
    axes.set_ylim(0, 140)
    axes.set_yticks(np.arange(0, 141, 20))
    axes.set_ylabel("Number of unique hits", fontsize=11)
    axes.set_xticks(np.arange(len(DATASET_LABELS)), ["Broad 8K", "Broad 39K", "ECBD", "CO-ADD", "Broad\npurchases"])
    axes.spines[["top", "right"]].set_visible(False)
    axes.tick_params(axis="x", length=0, pad=8)
    axes.tick_params(axis="y", direction="out", length=3.5, width=0.8)
    axes.legend(handles=[Patch(facecolor=STRAIN_COLORS[strain], edgecolor="black", linewidth=0.75, label=STRAIN_LABELS[strain]) for strain in STRAIN_ORDER],
                loc="lower center", bbox_to_anchor=(0.5, 1.08), ncol=5, frameon=False, handlelength=1.3, handletextpad=0.55, columnspacing=1.6)
    figure.text(0.11, 0.055, "Distinct structures per dataset; 80% growth inhibition hit calls", fontsize=8, color="#444444")
    for extension in ("svg", "pdf", "png"):
        options = {"dpi": 300} if extension == "png" else {}
        if extension == "svg":
            options["metadata"] = {"Date": None}
        figure.savefig(DIRECTORY / f"{FIGURE_STEM}.{extension}", facecolor="white", **options)
    pd.DataFrame(marker_positions).to_csv(DIRECTORY / "plotted_bar_coordinates.csv", index=False)
    plt.close(figure)


if __name__ == "__main__":
    plot_counts(count_hits())
