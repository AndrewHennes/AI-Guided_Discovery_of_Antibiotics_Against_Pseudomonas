"""Score CO-ADD with the supplied GI60 and GI80 ensembles."""
from pathlib import Path
import sys,json,hashlib
import pandas as pd,numpy as np,torch
DIRECTORY=Path(__file__).resolve().parent
METHODS=DIRECTORY.parent
sys.path.insert(0,str(METHODS))
from shared.figure2_models import DeepFlux
from shared.figure2_analysis import add_roc_columns
DATASET_PATH=''


def main():
    source=Path(DATASET_PATH).expanduser() if DATASET_PATH else DIRECTORY/'saved_inputs/COADD_input.csv.gz'
    if not source.is_absolute():raise ValueError('Use an absolute dataset path.')
    frame=pd.read_csv(source,dtype=str,keep_default_na=False)
    smiles=sorted(set(frame.Standardized_SMILES)-{''})
    index={s:i for i,s in enumerate(smiles)}
    features=np.full((len(smiles),512),np.nan,dtype=np.float32)
    for path in [METHODS/'2a/saved_inputs/COADD_minimol_features.npz',DIRECTORY/'saved_inputs/additional_minimol_features.npz']:
        with np.load(path,allow_pickle=False) as cache:
            for s,vector in zip(cache['smiles'],cache['features']):
                if s in index:features[index[s]]=vector
    valid=np.isfinite(features).all(axis=1)
    positions=np.array([index.get(s,-1) for s in frame.Standardized_SMILES])
    rowvalid=np.array([i>=0 and valid[i] for i in positions])
    frame['Source_Row_Index']=np.arange(len(frame))
    frame['DeepFlux_Scoring_Status']=np.where(rowvalid,'scored','missing_representation')
    torch.set_num_threads(8)
    x=torch.from_numpy(features[valid])
    checkpoints=[]
    for threshold in ['GI60','GI80']:
        member_scores=[]
        for number in range(10):
            path=METHODS/f'2a/checkpoints/{threshold}/ensemble_member_{number}.ckpt'
            checkpoint=torch.load(path,map_location='cpu',weights_only=True)
            model=DeepFlux(checkpoint['hyper_parameters'])
            model.load_state_dict(checkpoint['state_dict'],strict=True)
            model.eval()
            scores=np.full(len(smiles),np.nan)
            with torch.inference_mode():scores[valid]=np.concatenate([torch.sigmoid(model(batch)).numpy().ravel() for batch in x.split(4096)])
            member_scores.append(scores)
            values=np.full(len(frame),np.nan)
            values[rowvalid]=scores[positions[rowvalid]]
            frame[f'DeepFlux_{threshold}_member_{number}']=values
            checkpoints.append({'model':threshold,'member':number,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
        for suffix,values in [('mean',np.stack(member_scores).mean(axis=0)),('std',np.stack(member_scores).std(axis=0,ddof=0))]:
            rows=np.full(len(frame),np.nan);rows[rowvalid]=values[positions[rowvalid]]
            frame[f'DeepFlux_{threshold}_{suffix}']=rows
    frame=add_roc_columns(frame)
    frame.to_csv(DIRECTORY/'COADD_unfiltered_DeepFlux_GI60_GI80_scored.csv',index=False)
    cohort=frame.loc[frame.ROC_Included]
    if not np.isfinite(cohort[['DeepFlux_GI60_mean','DeepFlux_GI80_mean']]).all().all():raise ValueError('Missing scores in the ROC cohort.')
    cohort.to_csv(DIRECTORY/'COADD_unfiltered_Delta5_hits_and_inactive_ROC_cohort.csv',index=False)
    cohort.to_csv(DIRECTORY/'saved_inputs/COADD_roc_input_scores.csv.gz',index=False)
    report={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'rows':len(frame),'scored_rows':int(rowvalid.sum()),'checkpoints':checkpoints}
    (DIRECTORY/'rescoring_report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
