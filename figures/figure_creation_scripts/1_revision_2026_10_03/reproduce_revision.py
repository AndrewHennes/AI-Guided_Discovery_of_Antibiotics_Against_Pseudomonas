"""Reproduce Figure 1 additions and updated counts without command-line options."""
from pathlib import Path
from functools import lru_cache
import os,json,hashlib,sys
ROOT = Path('/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision')
INPUT = Path(__file__).resolve().parent / 'saved_inputs'
SOURCE_8K = Path('/Users/asselism/Desktop/Collins_Lab/Repositories/AI-Guided_Antipseudomonal_Discovery/AI-Guided_Antipseudomonal_Discovery/datasets/8K_Dataset/Broad_8K_Antibiotic_Activity_Dataset_Filtered_.csv')
TABLES = ROOT / 'tables'
PANELS = ROOT / 'panels'
for folder in [TABLES,PANELS]: folder.mkdir(parents=True,exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR','/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/work/matplotlib_cache')
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import rdkit
from rdkit import Chem,RDLogger
from rdkit.Chem import Descriptors,Crippen,rdMolDescriptors,GraphDescriptors
from violin_panels import plot_inhibition_violins
from property_boxplots import plot_property_boxplots
from amine_descriptor import amine_counts,verify_examples
from alcohol_descriptor import alcohol_count,verify_alcohol_examples
RDLogger.DisableLog('rdApp.warning')
# RGB appearance sampled from the supplied Illustrator file's existing bars.
COLORS={'WT':(0.07547,0.65324,0.31118),'Delta40':(0.97905,0.63760,0.11945),'Delta3':(0.75422,0.23545,0.17137),'Delta6':(0.66175,0.37743,0.65107),'Delta5':(0.36083,0.75999,0.87846),'Inactive':'#B9B9BD'}
STRAINS=['WT','Delta40','Delta3','Delta6']
LABELS={'WT':'WT','Delta40':'ΔOMP40','Delta3':'ΔE-pump3','Delta6':'ΔE-pump6','Delta5':'ΔE-pump5'}
COLUMNS={'8K':{s:f'{s}_Binarized_at_80_Percent_Growth_Inhibition' for s in STRAINS},'39K':{'WT':'PA_Hit_Binarized_at_80_Percent_Growth_Inhibition'},'ECBD':{'WT':'Pseudomonas_aeruginosa_Binarized_Hit'},'COADD':{'WT':'PA_Binarized_at_80_Percent_Growth_Inhibition','Delta5':'PA_MexB_Binarized_at_80_Percent_Growth_Inhibition'}}
PROPS=['Heavy_atoms','logP','TPSA','Log10_IPC','Amine_count','Aromatic_rings']
GROUPS=['WT_hit','Delta6_only_hit','Inactive_in_WT_and_Delta6']
GROUP_COLORS=[COLORS['WT'],COLORS['Delta6'],COLORS['Inactive']]
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.linewidth':.6,'xtick.major.width':.6,'ytick.major.width':.6,'xtick.major.size':2,'ytick.major.size':2,'pdf.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'figure1_revision_20261003','savefig.facecolor':'white'})

@lru_cache(None)
def canonical(smiles):
    if not isinstance(smiles,str) or not smiles.strip(): return None
    mol=Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol,isomericSmiles=True) if mol is not None else None

def read(name):
    if name == '8K':
        frame=pd.read_csv(SOURCE_8K,usecols=lambda c:c!='Minimol_Representation',low_memory=False)
        frame.insert(0,'Input_Row_Index',np.arange(len(frame)))
    else:
        frame=pd.read_csv(INPUT/f'{name}.csv.gz',low_memory=False)
    col='Standardized_SMILES' if 'Standardized_SMILES' in frame else 'SMILES'
    frame['Canonical_SMILES']=frame[col].map(canonical)
    if frame.Canonical_SMILES.isna().any():
        frame.loc[frame.Canonical_SMILES.isna()].to_csv(TABLES/f'{name}_invalid_structures.csv',index=False)
    return frame

def save(fig,stem):
    for ext in ['pdf','svg','png']:
        fig.savefig(PANELS/f'{stem}.{ext}',dpi=240,metadata={'Creator':'Figure 1 revision'} if ext=='pdf' else None)
    plt.close(fig)

def style(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(labelsize=6,pad=1.5)

def count_all():
    counts=[];hit_rows=[];audit=[]
    for name in COLUMNS:
        f=read(name); valid=f.loc[f.Canonical_SMILES.notna()]
        grouped=valid.groupby('Canonical_SMILES',sort=True)
        audit.append({'Dataset':name,'Source_records':len(f),'Unique_structures':valid.Canonical_SMILES.nunique(),'Invalid_structures':int(f.Canonical_SMILES.isna().sum())})
        for strain,col in COLUMNS[name].items():
            assert f[col].dropna().isin([0,1]).all(),(name,col)
            called=grouped[col].max();keys=called.index[called.eq(1)]
            counts.append({'Dataset':name,'Strain':strain,'Unique_hit_structures':len(keys),'Hit_records':int(f[col].eq(1).sum()),'Unique_structures_with_binary_call':int(called.notna().sum()),'Source_column':col})
            for key in keys:
                rows=valid.loc[valid.Canonical_SMILES.eq(key)&valid[col].eq(1)]
                hit_rows.append({'Dataset':name,'Strain':strain,'Canonical_SMILES':key,'Input_row_indices':';'.join(map(str,rows.Input_Row_Index))})
        print('Counted',name,flush=True)
    c=pd.DataFrame(counts);c.to_csv(TABLES/'updated_hit_counts.csv',index=False)
    pd.DataFrame(hit_rows).to_csv(TABLES/'counted_hit_structures.csv',index=False)
    pd.DataFrame(audit).to_csv(TABLES/'dataset_audit.csv',index=False)
    return c

def plot_counts(c):
    fig,ax=plt.subplots(figsize=(225/72,123/72));fig.subplots_adjust(left=.17,right=.99,bottom=.30,top=.79)
    ds=['8K','39K','ECBD','COADD']
    for i,d in enumerate(ds):
        z=c.loc[c.Dataset.eq(d)];width=.145
        for j,r in enumerate(z.itertuples()):
            x=i+(j-(len(z)-1)/2)*.17
            ax.bar(x,r.Unique_hit_structures,width,color=COLORS[r.Strain],edgecolor='black',linewidth=.4)
            ax.text(x,r.Unique_hit_structures+2,str(r.Unique_hit_structures),ha='center',va='bottom',fontsize=7)
    ax.set_ylim(0,140);ax.set_yticks([0,40,80,120]);ax.set_ylabel('Number of unique hits',fontsize=7,labelpad=2)
    ax.set_xticks(range(4),['Internal 8K','Internal 39K','ECBD','CO-ADD'],rotation=28,ha='right',fontsize=7)
    ax.tick_params(axis='x',length=0);style(ax)
    fig.legend([Patch(facecolor=COLORS[s],edgecolor='black',linewidth=.4) for s in [*STRAINS,'Delta5']], [LABELS[s] for s in [*STRAINS,'Delta5']],loc='upper center',bbox_to_anchor=(.57,1.01),ncol=3,frameon=False,fontsize=6,handlelength=1,columnspacing=.75,handletextpad=.3,labelspacing=.3)
    save(fig,'d_updated_compound_counts')

def inhibition(f):
    cols=[f'Average_{s}_Percent_Growth_Inhibition' for s in STRAINS]
    valid=f.loc[f.Canonical_SMILES.notna()]
    # Each strain uses its available continuous measurements. Never infer a
    # continuous inhibition value from a binary activity label.
    values=valid.groupby('Canonical_SMILES',sort=True)[cols].mean()
    values.to_csv(TABLES/'inhibition_values_by_structure.csv')
    summary=[]
    for strain,c in zip(STRAINS,cols):
        v=values[c].dropna()
        calls=valid.groupby('Canonical_SMILES')[COLUMNS['8K'][strain]].max()
        summary.append({'Strain':strain,'Measured_unique_structures':len(v),'Hits_at_80_percent':int(v.ge(80).sum()),'Hit_percent':100*v.ge(80).mean(),'Median':v.median(),'Available_records':int(f[c].notna().sum()),'Missing_continuous_records':int(f[c].isna().sum()),'Binary_hit_structures':int(calls.eq(1).sum())})
        assert int(v.ge(80).sum())==int(calls.eq(1).sum())
    summary=pd.DataFrame(summary);summary.to_csv(TABLES/'inhibition_summary.csv',index=False)
    plot_inhibition_violins(values,summary,STRAINS,LABELS,COLORS,save,TABLES)
    f[['Input_Row_Index','Canonical_SMILES',*cols]].to_csv(TABLES/'all_available_inhibition_values.csv',index=False)
    return summary

def properties(f):
    verify_examples()
    verify_alcohol_examples()
    wt='Average_WT_Percent_Growth_Inhibition';d6='Average_Delta6_Percent_Growth_Inhibition'
    binary=[COLUMNS['8K']['WT'],COLUMNS['8K']['Delta6']]
    assert f[binary].notna().all().all()
    assert f[binary].isin([0,1]).all().all()
    grouped=f.loc[f.Canonical_SMILES.notna()].groupby('Canonical_SMILES',sort=True)
    cohort=[]
    for key,rows in grouped:
        w=rows[wt].dropna();d=rows[d6].dropna()
        wt_call=int(rows[binary[0]].max()); d6_call=int(rows[binary[1]].max())
        if wt_call==1: category=GROUPS[0]
        elif d6_call==1: category=GROUPS[1]
        else: category=GROUPS[2]
        m=Chem.MolFromSmiles(key)
        ipc=GraphDescriptors.Ipc(m,avg=False)
        assert np.isfinite(ipc) and ipc>0,(key,ipc)
        cohort.append({'Canonical_SMILES':key,'Category':category,'WT_activity_label':wt_call,'Delta6_activity_label':d6_call,'Input_row_indices':';'.join(map(str,rows.Input_Row_Index)),'Catalog_IDs':';'.join(rows['Catalog ID'].astype(str)),'WT_max_measured_inhibition':w.max() if len(w) else np.nan,'Delta6_max_measured_inhibition':d.max() if len(d) else np.nan,'MW':Descriptors.MolWt(m),'Heavy_atoms':m.GetNumHeavyAtoms(),'Alcohol_count':alcohol_count(m),'logP':Crippen.MolLogP(m),'TPSA':rdMolDescriptors.CalcTPSA(m,includeSandP=False),'Formal_charge':Chem.GetFormalCharge(m),'IPC':ipc,'Log10_IPC':np.log10(ipc),'Rotatable_bonds':rdMolDescriptors.CalcNumRotatableBonds(m,rdMolDescriptors.NumRotatableBondsOptions.Strict),'Aromatic_rings':rdMolDescriptors.CalcNumAromaticRings(m),**amine_counts(m)})
    p=pd.DataFrame(cohort);p.to_csv(TABLES/'all_8k_molecular_properties.csv',index=False)
    p.groupby(['Category','Amine_count']).size().rename('Compound_count').reset_index().to_csv(TABLES/'amine_count_distribution.csv',index=False)
    p.groupby(['Category','Alcohol_count']).size().rename('Compound_count').reset_index().to_csv(TABLES/'alcohol_count_distribution.csv',index=False)
    hits=p.loc[p.Delta6_activity_label.eq(1)].copy()
    assert len(hits)==124
    assert p.Category.value_counts().to_dict()=={GROUPS[2]:7861,GROUPS[1]:116,GROUPS[0]:8}
    assert p.loc[p.WT_activity_label.eq(1),'Delta6_activity_label'].eq(1).all()
    hits.to_csv(TABLES/'8K_WT_and_Delta6_only_hit_properties.csv',index=False)
    p.loc[p.Category.eq(GROUPS[1])].to_csv(TABLES/'8K_Delta6_only_116_compounds.csv',index=False)
    p.loc[p.WT_max_measured_inhibition.isna()].to_csv(TABLES/'structures_with_missing_continuous_WT.csv',index=False)
    summary=[];tests=[]
    for prop in PROPS:
        for g in GROUPS:
            a=p.loc[p.Category.eq(g),prop]
            summary.append({'Property':prop,'Group':g,'N':len(a),'Median':a.median(),'Q1':a.quantile(.25),'Q3':a.quantile(.75),'Min':a.min(),'Max':a.max()})
        for first,second in [(0,1),(0,2),(1,2)]:
            a=p.loc[p.Category.eq(GROUPS[first]),prop];b=p.loc[p.Category.eq(GROUPS[second]),prop]
            u,pval=stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic')
            tests.append({'Property':prop,'Group_A':GROUPS[first],'Group_B':GROUPS[second],'N_A':len(a),'N_B':len(b),'Mann_Whitney_U':u,'p_value':pval,'Rank_biserial_A_minus_B':2*u/(len(a)*len(b))-1})
    pd.DataFrame(summary).to_csv(TABLES/'property_summaries.csv',index=False)
    t=pd.DataFrame(tests);t['BH_q_value']=stats.false_discovery_control(t.p_value.to_numpy(),method='bh');t.to_csv(TABLES/'exploratory_property_statistics.csv',index=False)
    counts=plot_property_boxplots(p,PROPS,GROUPS,GROUP_COLORS,save,TABLES)
    return counts,p

def main():
    counts=count_all();plot_counts(counts)
    f=read('8K');summary=inhibition(f)
    group_counts,p=properties(f)
    report={'version':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'rdkit':rdkit.__version__},'source_records_8K':len(f),'property_group_counts':dict(zip(GROUPS,group_counts)),'missing_continuous_WT_records':int(f.Average_WT_Percent_Growth_Inhibition.isna().sum()),'excluded_property_structures':0,'source_8K_path':str(SOURCE_8K),'source_8K_sha256':hashlib.sha256(SOURCE_8K.read_bytes()).hexdigest(),'older_hit_property_table_used':False,'colors_RGB':COLORS,'original_panel_b_annotation_n':7994,'panel_f_screen_annotation_n':7994,'panel_f_density_scale':'Shared log10 density width above 1e-5 floor','panel_g_fourth_descriptor':'log10(GraphDescriptors.Ipc(mol,avg=False))','panel_g_plot_type':'Box-and-whisker plots with 1.5 IQR whiskers','panel_g_first_descriptor':'Heavy atom count using Mol.GetNumHeavyAtoms()','panel_g_fifth_descriptor':'Primary, secondary, and tertiary amine nitrogens including protonated forms','original_panel_b_preserved':True,'hit_count_rule':'Canonical isomeric structures, source binary hit calls combined by OR.','property_rule':'Use supplied WT and Delta6 binary activity labels. WT=1 is WT hit; WT=0 and Delta6=1 is Delta6-only; both 0 is inactive.','inhibition_rule':'All available continuous measurements per strain, averaged within canonical structure; missing continuous values are not imputed.'}
    (TABLES/'validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
