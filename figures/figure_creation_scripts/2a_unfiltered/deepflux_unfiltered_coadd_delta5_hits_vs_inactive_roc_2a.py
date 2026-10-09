"""Compute ROC curves for Delta5 hits and inactive CO-ADD records."""
from pathlib import Path
import sys,json
import pandas as pd
import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score, roc_curve

DIRECTORY = Path(__file__).resolve().parent
METHODS = DIRECTORY.parent
sys.path.insert(0, str(METHODS))
from shared.figure2_analysis import add_roc_columns, ROC_DEFINITION as MAIN_ROC_DEFINITION
from shared.figure2_roc_panel import render_roc_panel

SCORED_DATASET_PATH = ''
ROC_DEFINITION = {
    **MAIN_ROC_DEFINITION,
    'analysis_unit': 'input_dataset_record',
    'dataset_variant': 'unfiltered',
    'deduplication_key': None,
    'representative_selection': 'each eligible input record',
}


def main():
    source = Path(SCORED_DATASET_PATH).expanduser() if SCORED_DATASET_PATH else DIRECTORY / 'saved_inputs/COADD_roc_input_scores.csv.gz'
    if not source.is_absolute():
        raise ValueError('Use an absolute dataset path.')
    frame = add_roc_columns(pd.read_csv(source, low_memory=False))
    frame = frame.loc[frame.ROC_Included].copy()
    y = frame.ROC_Delta5_Hit_Positive_Label.to_numpy(dtype=int)
    curves, coordinates, metrics = {}, [], []
    for model in ['GI60', 'GI80']:
        scores = frame[f'DeepFlux_{model}_mean'].to_numpy(dtype=float)
        if not np.isfinite(scores).all():
            raise ValueError('Each included record must have a finite score.')
        fpr, tpr, thresholds = roc_curve(y, scores, drop_intermediate=False)
        pfpr, ptpr, pthresholds = roc_curve(y, scores, drop_intermediate=True)
        area = float(roc_auc_score(y, scores))
        n1, n0 = int(y.sum()), int((1-y).sum())
        rank_area = float((rankdata(scores)[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
        assert np.isclose(area, rank_area, atol=1e-12, rtol=0)
        assert np.isclose(area, np.trapezoid(ptpr, pfpr), atol=1e-12, rtol=0)
        curves[model] = {'fpr': pfpr, 'tpr': ptpr, 'threshold': pthresholds, 'auc': area}
        coordinates.extend({'Model': model, 'False_Positive_Rate': f, 'True_Positive_Rate': t, 'Score_Threshold': s} for f,t,s in zip(fpr,tpr,thresholds))
        metrics.append({'Model': model, 'AUC': area, 'Positive_Records': n1, 'Negative_Records': n0, 'Independent_Rank_AUC': rank_area, 'Score_Transformation': 'original ensemble_mean'})
    pd.DataFrame(metrics).to_csv(DIRECTORY / 'ROC_metrics.csv', index=False)
    pd.DataFrame(coordinates).to_csv(DIRECTORY / 'ROC_curve_coordinates.csv', index=False)
    (DIRECTORY / 'ROC_definition.json').write_text(json.dumps(ROC_DEFINITION, indent=2) + '\n')
    render_roc_panel(curves, DIRECTORY / 'deepflux_unfiltered_coadd_delta5_hits_vs_inactive_roc_2a')
    print(pd.DataFrame(metrics).to_string(index=False))


if __name__ == '__main__':
    main()
