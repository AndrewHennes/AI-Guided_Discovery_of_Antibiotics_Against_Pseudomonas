"""Independently verify structure counts and the source activity-label cohorts."""
from pathlib import Path
import json
import pandas as pd
from rdkit import Chem
ROOT=Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
INPUT=Path(__file__).resolve().parent/'saved_inputs'
SOURCE_8K=Path('/Users/asselism/Desktop/Collins_Lab/Repositories/AI-Guided_Antipseudomonal_Discovery/AI-Guided_Antipseudomonal_Discovery/datasets/8K_Dataset/Broad_8K_Antibiotic_Activity_Dataset_Filtered_.csv')
counts=pd.read_csv(ROOT/'tables/updated_hit_counts.csv');checks=[]
for name,rows in counts.groupby('Dataset'):
    f=pd.read_csv(SOURCE_8K,usecols=lambda c:c!='Minimol_Representation',low_memory=False) if name=='8K' else pd.read_csv(INPUT/f'{name}.csv.gz',low_memory=False)
    column='Standardized_SMILES' if 'Standardized_SMILES' in f else 'SMILES'
    for row in rows.itertuples():
        keys=set()
        for s in f.loc[f[row.Source_column].eq(1),column]:
            m=Chem.MolFromSmiles(s);assert m is not None
            keys.add(Chem.MolToSmiles(m,canonical=True,isomericSmiles=True))
        assert len(keys)==row.Unique_hit_structures,(name,row.Strain)
        checks.append({'Dataset':name,'Strain':row.Strain,'Expected_from_positive_records':len(keys),'Exported':row.Unique_hit_structures,'Passed':True})
        if name=='8K' and row.Strain=='WT':wt_keys=keys
        if name=='8K' and row.Strain=='Delta6':d6_keys=keys
assert len(wt_keys)==8 and len(d6_keys)==124 and wt_keys.issubset(d6_keys)
assert len(d6_keys-wt_keys)==116
p=pd.read_csv(ROOT/'tables/all_8k_molecular_properties.csv')
assert not p.Canonical_SMILES.duplicated().any()
assert p.Category.value_counts().to_dict()=={'Inactive_in_WT_and_Delta6':7861,'Delta6_only_hit':116,'WT_hit':8}
assert set(p.loc[p.Category.eq('WT_hit'),'Canonical_SMILES'])==wt_keys
assert set(p.loc[p.Category.eq('Delta6_only_hit'),'Canonical_SMILES'])==d6_keys-wt_keys
assert p.loc[p.Category.eq('Delta6_only_hit'),'WT_max_measured_inhibition'].isna().sum()==3
assert len(pd.read_csv(ROOT/'tables/8K_WT_and_Delta6_only_hit_properties.csv'))==124
v=pd.read_csv(ROOT/'tables/inhibition_values_by_structure.csv')
summary=pd.read_csv(ROOT/'tables/inhibition_summary.csv')
for row in summary.itertuples():
    a=v[f'Average_{row.Strain}_Percent_Growth_Inhibition'].dropna()
    assert len(a)==row.Measured_unique_structures
    assert a.ge(80).sum()==row.Hits_at_80_percent==row.Binary_hit_structures
pd.DataFrame(checks).to_csv(ROOT/'tables/independent_count_validation.csv',index=False)
(ROOT/'tables/analysis_checks.json').write_text(json.dumps({'counts_match_independent_positive_record_check':True,'WT_hits':8,'Delta6_hits':124,'overlap':8,'Delta6_only_hits':116,'inactive_in_both':7861,'three_missing_continuous_WT_compounds_included_as_Delta6_only':True,'property_cohorts_cover_all_7985_structures':True,'older_property_table_used':False,'inhibition_counts_match_binary_hit_counts':True},indent=2))
print('Passed: WT 8; Delta6 124; overlap 8; Delta6-only 116; inactive 7861. All three compounds retained.')
