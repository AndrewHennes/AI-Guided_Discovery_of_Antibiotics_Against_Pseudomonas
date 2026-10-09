"""Compare numerical score pairs between the workbook and scored CO-ADD."""
from pathlib import Path
import json
import pandas as pd,numpy as np
from scipy.spatial import cKDTree
from sklearn.metrics import roc_auc_score
DIRECTORY=Path(__file__).resolve().parent
SCORED_DATASET_PATH=''
WORKBOOK_PATH=''


def main():
    source=Path(SCORED_DATASET_PATH).expanduser() if SCORED_DATASET_PATH else DIRECTORY.parent/'2a_unfiltered/COADD_unfiltered_DeepFlux_GI60_GI80_scored.csv'
    workbook=Path(WORKBOOK_PATH).expanduser() if WORKBOOK_PATH else DIRECTORY/'saved_inputs/workbook_scores_and_labels.csv.gz'
    if not source.is_absolute() or not workbook.is_absolute():raise ValueError('Use absolute input paths.')
    w=pd.read_excel(workbook) if workbook.suffix=='.xlsx' else pd.read_csv(workbook)
    cols=['DeepFlux_GI60_mean','DeepFlux_GI80_mean']
    f=pd.read_csv(source,usecols=['Standardized_SMILES','COADD_Compounds_ID',*cols])
    u=f.drop_duplicates('Standardized_SMILES').dropna(subset=cols)
    x=u[cols].to_numpy();y=w[['GI_60_model_mean','GI_80_model_mean']].to_numpy()
    distances,indices=cKDTree(np.log(x+0.001)).query(np.log(y+0.001),k=2,p=np.inf)
    index=indices[:,0]
    match=np.isclose(x[index],y,rtol=1e-5,atol=6e-9).all(axis=1)
    w.insert(0,'Workbook_Row_Number',np.arange(len(w))+2)
    w['Close_Score_Pair']=match
    w['Candidate_Standardized_SMILES']=np.where(match,u.Standardized_SMILES.to_numpy()[index],'')
    for j,n in enumerate([60,80]):w[f'Score_Difference_GI{n}']=x[index,j]-y[:,j]
    w.to_csv(DIRECTORY/'workbook_numerical_score_pair_comparison.csv',index=False)
    report={'workbook_records':len(w),'close_score_pairs':int(match.sum()),'relative_tolerance':1e-5,'absolute_tolerance':6e-9,'matching_basis':'numerical score-pair comparison','models':{}}
    for j,n in enumerate([60,80]):
        old=y[match,j];new=x[index[match],j];labels=w.loc[match,'Binarized_at_80_PGI']
        report['models'][f'GI{n}']={'maximum_absolute_score_difference':float(np.abs(new-old).max()),'old_scores_workbook_80_label_auc':float(roc_auc_score(labels,old)),'new_scores_workbook_80_label_auc':float(roc_auc_score(labels,new))}
    (DIRECTORY/'score_comparison_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
