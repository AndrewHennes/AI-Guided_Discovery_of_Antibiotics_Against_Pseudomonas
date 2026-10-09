"""Canonical parent-structure preparation for chemical-space analysis."""

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

STANDARDIZATION = {
    "input_column": "SMILES in all four source datasets",
    "steps": ["RDKit sanitized SMILES parsing", "rdMolStandardize.Cleanup", "rdMolStandardize.FragmentParent with skipStandardize=True", "rdMolStandardize.Uncharger", "RemoveHs", "canonical isomeric SMILES serialization"],
    "tautomer_canonicalization": False,
    "stereochemistry": "Retained in standardized isomeric SMILES",
    "convergence": "Repeat cleanup and parent preparation to a fixed point; choose the lexicographically first canonical SMILES if charge redistribution cycles",
    "deduplication": "One joint embedding coordinate per standardized structure",
}


def initialize_worker():
    RDLogger.DisableLog("rdApp.*")


def standardize(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        return None, None, "missing_smiles"
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None, None, "invalid_smiles"
    original = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    current, states = original, []
    for _ in range(12):
        if current in states:
            return original, min(states[states.index(current):]), "ok"
        states.append(current)
        molecule = rdMolStandardize.Cleanup(Chem.MolFromSmiles(current))
        molecule = rdMolStandardize.FragmentParent(molecule, skipStandardize=True)
        molecule = rdMolStandardize.Uncharger().uncharge(molecule)
        molecule = Chem.RemoveHs(molecule)
        Chem.SanitizeMol(molecule)
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
        current = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    raise ValueError("Molecular standardization did not converge")


def standardize_record(smiles):
    try:
        return standardize(smiles)
    except Exception as error:
        return None, None, f"standardization_error: {type(error).__name__}: {error}"
