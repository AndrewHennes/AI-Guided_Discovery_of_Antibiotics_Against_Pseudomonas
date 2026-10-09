"""Count neutral alcohol OH groups bonded to saturated carbon atoms.

Phenols, enols, carboxylic acids, water, and alkoxide ions are not counted.
Each distinct hydroxyl oxygen is counted once, including multiple OH groups
on a single carbon atom.
"""
from rdkit import Chem

ALCOHOL_SMARTS = '[OX2H1;+0]-[CX4]'
ALCOHOL_QUERY = Chem.MolFromSmarts(ALCOHOL_SMARTS)
assert ALCOHOL_QUERY is not None


def alcohol_count(mol):
    return len({match[0] for match in mol.GetSubstructMatches(ALCOHOL_QUERY, uniquify=True)})


def verify_alcohol_examples():
    examples = {
        'CO': 1, 'CCO': 1, 'CC(O)C': 1, 'CC(C)(C)O': 1,
        'OCCO': 2, 'OCC(O)CO': 3, 'C(O)O': 2,
        'OCc1ccccc1': 1, 'C=CCO': 1, 'CC(O)C(=O)O': 1,
        'COCO': 1, 'O[CH2][CH2]O': 2,
        'Oc1ccccc1': 0, 'Oc1ccccc1CO': 1, 'CC(=O)O': 0,
        'O=C(O)C(=O)O': 0, 'C=CO': 0, 'COC': 0,
        'COC(=O)C': 0, 'COO': 0, 'NCO': 1, 'NO': 0,
        'C[O-]': 0, 'O': 0, '[OH-]': 0, 'C': 0,
    }
    for smiles, expected in examples.items():
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, smiles
        actual = alcohol_count(mol)
        assert actual == expected, (smiles, expected, actual)
    return len(examples)


if __name__ == '__main__':
    print(f'Passed {verify_alcohol_examples()} alcohol functional-group examples.')
