"""Count primary, secondary and tertiary amine nitrogens, including protonated forms.

Uses the exclusions in RDKit's functional-group Amine definition, with explicit
valence/charge restrictions to exclude quaternary ammonium and imine nitrogens.
"""
from rdkit import Chem

AMINE_SMARTS = ('[$([N;X3;+0]),$([N;X4;+1;!H0]);$(N-[#6]);'
                '!$(N-[!#6;!#1]);!$(N-C=[O,N,S]);!$(N-C#N);!$(N=*);!$(N#*)]')
AMINE_QUERY = Chem.MolFromSmarts(AMINE_SMARTS)
assert AMINE_QUERY is not None


def amine_counts(mol):
    counts = {'Primary_amines': 0, 'Secondary_amines': 0, 'Tertiary_amines': 0}
    labels = {1: 'Primary_amines', 2: 'Secondary_amines', 3: 'Tertiary_amines'}
    for (index,) in mol.GetSubstructMatches(AMINE_QUERY, uniquify=True):
        nitrogen = mol.GetAtomWithIdx(index)
        carbons = sum(a.GetAtomicNum() == 6 for a in nitrogen.GetNeighbors())
        assert carbons in labels, (index, carbons)
        counts[labels[carbons]] += 1
    counts['Amine_count'] = sum(counts.values())
    return counts


def verify_examples():
    examples = {
        'CN': 1, 'CNC': 1, 'CN(C)C': 1, 'C[NH3+]': 1,
        'C[NH2+]C': 1, 'C[NH+](C)C': 1, 'CNCCN': 2,
        'c1ccccc1N': 1, 'c1ccccc1[NH3+]': 1, 'N1CCCCC1': 1,
        'N1CCNCC1': 2, 'NCc1ccncc1': 1, 'CN(C)c1ccccc1': 1,
        'CC(=O)N': 0, 'CC(=O)NC': 0, 'COC(=O)N': 0, 'NC(=O)N': 0,
        'CC(=S)N': 0, 'CS(=O)(=O)N': 0, 'NC(=N)N': 0,
        'CC(=N)N': 0, 'N=CNC': 0, 'N#CN': 0, 'CC#N': 0,
        'c1ccncc1': 0, 'c1cc[nH]c1': 0, 'C[N+](C)(C)C': 0,
        'C[N+](C)(C)[O-]': 0, 'C[N+](=O)[O-]': 0, 'CNO': 0,
        'CNN': 0, 'N': 0, '[NH4+]': 0,
    }
    for smiles, expected in examples.items():
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, smiles
        actual = amine_counts(mol)['Amine_count']
        assert actual == expected, (smiles, expected, actual)
    return len(examples)


if __name__ == '__main__':
    print(f'Passed {verify_examples()} amine functional-group examples.')
