"""Build a self-contained molecular explorer of Pseudomonas hits."""

import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

from plot_unique_pseudomonas_hit_counts import (
    DATASET_LABELS, ID_COLUMNS, STRAIN_COLUMNS, STRAIN_COLORS,
    STRAIN_LABELS, STRAIN_ORDER, canonical_structure,
)

DIRECTORY = Path(__file__).resolve().parent
DATASET_PATHS = {"8K": "", "39K": "", "ECBD": "", "COADD": "", "PurchasedBroad": ""}
OUTPUT_PATH = DIRECTORY / "Pseudomonas_hit_explorer.html"
INHIBITION_COLUMNS = {
    "8K": {strain: f"Average_{strain}_Percent_Growth_Inhibition" for strain in ("WT", "Delta40", "Delta3", "Delta6")},
    "39K": {"WT": "PA_Percent_Growth_Inhibition_Mean"},
    "ECBD": {"WT": "Pseudomonas_aeruginosa_Percent_Growth_Inhibition"},
    "COADD": {"WT": "PA_Growth_Inhibition", "Delta5": "PA_MexB_Growth_Inhibition"},
    "PurchasedBroad": {"WT": "Percent_Growth_Inhibition_OD600_WT", "Delta6": "Percent_Growth_Inhibition_OD600_∆6_Strain"},
}


def strings(series):
    return sorted({str(value).strip() for value in series.dropna() if str(value).strip()})


def load_records(dataset):
    external_path = DATASET_PATHS[dataset]
    if external_path:
        path = Path(external_path).expanduser()
        if not path.is_absolute():
            raise ValueError("Dataset paths must be absolute.")
        frame = pd.read_excel(path, sheet_name="Data") if dataset == "PurchasedBroad" else pd.read_csv(path, low_memory=False, float_precision="round_trip")
        frame["Source_Row_Index"] = np.arange(len(frame))
    else:
        frame = pd.read_csv(DIRECTORY / "saved_inputs" / f"{dataset}_explorer_records.csv.gz", low_memory=False, float_precision="round_trip")
    if dataset == "PurchasedBroad":
        for strain, column in STRAIN_COLUMNS[dataset].items():
            values = pd.to_numeric(frame[INHIBITION_COLUMNS[dataset][strain]], errors="raise")
            frame[column] = values.gt(80).astype(float).where(values.notna())
    for column in STRAIN_COLUMNS[dataset].values():
        frame[column] = pd.to_numeric(frame[column], errors="raise")
        if not frame[column].dropna().isin([0, 1]).all():
            raise ValueError(f"{dataset}: nonbinary value in {column}")
    smiles_column = "Standardized_SMILES" if "Standardized_SMILES" in frame else "SMILES"
    frame["Canonical_SMILES"] = frame[smiles_column].map(canonical_structure)
    return frame


def molecule_data(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    drawer = rdMolDraw2D.MolDraw2DSVG(480, 300)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.padding = 0.09
    options.bondLineWidth = 1.8
    options.minFontSize = 12
    options.maxFontSize = 22
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, molecule)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = svg[svg.index("<svg"):]
    svg = re.sub(r"width='480px' height='300px'", 'width="100%" height="100%"', svg)
    return {
        "smiles": smiles, "svg": svg,
        "mw": round(Descriptors.MolWt(molecule), 2),
        "logp": round(Descriptors.MolLogP(molecule), 2),
        "tpsa": round(rdMolDescriptors.CalcTPSA(molecule), 2),
        "hbd": Lipinski.NumHDonors(molecule), "hba": Lipinski.NumHAcceptors(molecule),
    }


def assay_data(rows, dataset, strain):
    if strain not in STRAIN_COLUMNS[dataset]:
        return {"status": "not_assayed", "hitRecords": 0, "nonhitRecords": 0, "missingRecords": 0, "inhibitionValues": [], "sourceStrainNames": []}
    calls = rows[STRAIN_COLUMNS[dataset][strain]]
    hit_count, nonhit_count = int(calls.eq(1).sum()), int(calls.eq(0).sum())
    values = pd.to_numeric(rows[INHIBITION_COLUMNS[dataset][strain]], errors="raise")
    strain_column = {"WT": "PA_Strain", "Delta5": "PA_MexB_Strain"}.get(strain) if dataset == "COADD" else None
    return {
        "status": "hit" if hit_count else "nonhit" if nonhit_count else "unmeasured",
        "hitRecords": hit_count, "nonhitRecords": nonhit_count, "missingRecords": int(calls.isna().sum()),
        "inhibitionValues": sorted({float(value) for value in values if pd.notna(value) and np.isfinite(value)}),
        "sourceStrainNames": strings(rows[strain_column]) if strain_column and strain_column in rows else [],
    }


def build_data():
    RDLogger.DisableLog("rdApp.*")
    manifest = json.loads((DIRECTORY / "saved_inputs" / "source_manifest.json").read_text())
    descriptions = {dataset: "Filtered dataset used in Figure 1" for dataset in DATASET_LABELS}
    descriptions["PurchasedBroad"] = "183 purchased Broad compounds with measured WT and Δ6 activity"
    data = {
        "datasets": [],
        "strains": [{"id": strain, "label": STRAIN_LABELS[strain], "color": STRAIN_COLORS[strain]} for strain in STRAIN_ORDER],
        "molecules": {}, "entries": [],
        "scope": "Pseudomonas aeruginosa hits in the four filtered screening datasets and the purchased Broad screening set. CO-ADD WT includes the source PAO1 and ATCC 27853 assays.",
        "countingRule": "One entry per canonical isomeric structure per dataset. Any existing binary call of 1 designates a hit at the 80% cutoff; purchased Broad calls use inhibition >80%. Mutant totals include WT overlaps. Source row indices are zero-based. Molecular properties are calculated with RDKit.",
    }
    for dataset, label in DATASET_LABELS.items():
        frame = load_records(dataset)
        source = dict(manifest[dataset])
        if DATASET_PATHS[dataset]:
            source_path = Path(DATASET_PATHS[dataset]).expanduser().resolve()
            source.update(filename=source_path.name, sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(), source_rows=len(frame))
        valid = frame.loc[frame.Canonical_SMILES.notna()]
        counts = {strain: 0 for strain in STRAIN_COLUMNS[dataset]}
        dataset_entries = []
        for smiles, rows in valid.groupby("Canonical_SMILES", sort=True):
            assays = {strain: assay_data(rows, dataset, strain) for strain in STRAIN_ORDER}
            hits = [strain for strain in STRAIN_ORDER if assays[strain]["status"] == "hit"]
            if not hits:
                continue
            mol_key = "mol_" + hashlib.sha256(smiles.encode()).hexdigest()[:20]
            if mol_key not in data["molecules"]:
                data["molecules"][mol_key] = molecule_data(smiles)
            elif data["molecules"][mol_key]["smiles"] != smiles:
                raise ValueError("Molecule key collision")
            id_column = ID_COLUMNS[dataset]
            identifiers = strings(rows[id_column]) if id_column else []
            chemical_names = strings(rows["Common_Compound_Name"]) if "Common_Compound_Name" in rows else []
            name = chemical_names[0] if chemical_names else identifiers[0] if identifiers else f"{label} structure {len(dataset_entries) + 1:02d}"
            entry = {
                "id": f"{dataset}_{mol_key}", "molKey": mol_key, "datasetId": dataset,
                "name": name, "chemicalName": " | ".join(chemical_names), "identifiers": identifiers,
                "hitStrains": hits, "assays": assays,
                "sourceRows": sorted(map(int, rows.Source_Row_Index)),
                "sourceFilename": source["filename"],
            }
            dataset_entries.append(entry)
            for strain in hits:
                counts[strain] += 1
        data["entries"].extend(dataset_entries)
        data["datasets"].append({
            "id": dataset, "label": label, "shortLabel": label,
            "description": descriptions[dataset], "filename": source["filename"],
            "sourceRows": source["source_rows"], "sourceSha256": source["sha256"],
            "sourceWorksheet": source.get("worksheet", ""), "compoundIdColumn": ID_COLUMNS[dataset],
            "structureSourceColumn": "Standardized_SMILES" if "Standardized_SMILES" in frame else "SMILES",
            "totalHits": len(dataset_entries), "counts": counts,
        })
    return data


def main():
    data = build_data()
    serialized = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    (DIRECTORY / "hit_explorer_data.json").write_text(serialized + "\n")
    template = (DIRECTORY / "hit_explorer_template.html").read_text()
    if template.count("__HIT_EXPLORER_DATA__") != 1:
        raise ValueError("Expected exactly one data placeholder in HTML template")
    OUTPUT_PATH.write_text(template.replace("__HIT_EXPLORER_DATA__", serialized.replace("<", "\\u003c")))
    expected = pd.read_csv(DIRECTORY / "unique_hit_counts_by_dataset_and_strain.csv")
    validation = {"datasets": {}, "dataset_molecule_entries": len(data["entries"]), "globally_distinct_structures": len(data["molecules"])}
    for dataset in data["datasets"]:
        validation["datasets"][dataset["id"]] = {"any_strain": dataset["totalHits"], **dataset["counts"]}
        for strain, count in dataset["counts"].items():
            reference = expected.loc[expected.Dataset.eq(dataset["id"]) & expected.Strain.eq(strain), "Unique_Hit_Structures"]
            if not any(DATASET_PATHS.values()) and (len(reference) != 1 or int(reference.iloc[0]) != count):
                raise ValueError(f"Count differs from audited table: {dataset['id']} {strain}")
    validation["wt_dataset_entries"] = sum("WT" in entry["hitStrains"] for entry in data["entries"])
    validation["globally_distinct_wt_structures"] = len({entry["molKey"] for entry in data["entries"] if "WT" in entry["hitStrains"]})
    validation["missing_structure_drawings"] = sum(not molecule["svg"].startswith("<svg") for molecule in data["molecules"].values())
    validation["html_bytes"] = OUTPUT_PATH.stat().st_size
    (DIRECTORY / "hit_explorer_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps(validation, indent=2))
    print(OUTPUT_PATH)


if __name__ == "__main__":
    main()
