"""Compare CO-ADD ROC curves at 60% and 80% growth inhibition."""
from pathlib import Path
import sys,json
import pandas as pd,numpy as np
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score,roc_curve
import pymupdf

DIRECTORY=Path(__file__).resolve().parent
METHODS=DIRECTORY.parent
sys.path.insert(0,str(METHODS))
from shared.figure2_roc_panel import render_roc_panel
from shared.rendering import compact_svg
WORKBOOK_PATH=''
FILTERED_SCORED_DATASET_PATH=''
UNFILTERED_SCORED_DATASET_PATH=''
REPEATED_STRUCTURE='c1ccc(C2=Nc3nonc3N=C(c3ccccc3)C2)cc1'


def source_path(value,default):
    path=Path(value).expanduser() if value else default
    if not path.is_absolute():
        raise ValueError('Use an absolute input path.')
    return path


def select_cohort(frame,cutoff):
    f=frame.copy()
    d5=pd.to_numeric(f.PA_MexB_Growth_Inhibition,errors='coerce')
    wt=pd.to_numeric(f.PA_Growth_Inhibition,errors='coerce')
    positive=d5.ge(cutoff)
    negative=d5.lt(cutoff)&wt.lt(cutoff)
    f['Label']=positive.astype(int)
    return f.loc[positive|negative].copy()


def curve(y,scores):
    y=np.asarray(y,dtype=int);scores=np.asarray(scores,dtype=float)
    if not np.isfinite(scores).all():raise ValueError('Missing ROC scores.')
    fpr,tpr,thresholds=roc_curve(y,scores,drop_intermediate=True)
    area=float(roc_auc_score(y,scores))
    n1=int(y.sum());n0=len(y)-n1
    rank_area=(rankdata(scores)[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
    assert abs(area-rank_area)<1e-12
    assert abs(area-np.trapezoid(tpr,fpr))<1e-12
    return {'fpr':fpr,'tpr':tpr,'threshold':thresholds,'auc':area}


def strip(panels,output_stem):
    document=pymupdf.open();page=document.new_page(width=176*len(panels),height=200)
    for i,(path,title,subtitle) in enumerate(panels):
        source=pymupdf.open(path)
        page.show_pdf_page(pymupdf.Rect(176*i+9,25,176*i+167,198),source,0)
        page.insert_text((176*i+(176-pymupdf.get_text_length(title,fontsize=8.5))/2,10),title,fontname='helv',fontsize=8.5)
        page.insert_text((176*i+(176-pymupdf.get_text_length(subtitle,fontsize=7))/2,21),subtitle,fontname='helv',fontsize=7)
    document.save(output_stem.with_suffix('.pdf'),garbage=3,deflate=True)
    output_stem.with_suffix('.svg').write_text(compact_svg(page.get_svg_image(text_as_path=True)),encoding='utf-8')
    page.get_pixmap(dpi=180,alpha=False).save(output_stem.with_suffix('.png'))
    document.close()


def main():
    wpath=source_path(WORKBOOK_PATH,DIRECTORY/'saved_inputs/workbook_scores_and_labels.csv.gz')
    workbook=pd.read_excel(wpath) if wpath.suffix=='.xlsx' else pd.read_csv(wpath)
    files={
        'filtered':source_path(FILTERED_SCORED_DATASET_PATH,METHODS/'2a/saved_inputs/COADD_roc_input_scores.csv.gz'),
        'unfiltered':source_path(UNFILTERED_SCORED_DATASET_PATH,METHODS/'2a_unfiltered/saved_inputs/COADD_roc_input_scores.csv.gz'),
    }
    frames={name:pd.read_csv(path,low_memory=False) for name,path in files.items()}
    metrics=[];coordinates=[];curves={};cohorts={}
    for threshold in [60,80]:
        for dataset,f in frames.items():
            selected=select_cohort(f,threshold)
            variants={'all_records':selected}
            for unit,key in [('one_per_COADD_ID','COADD_Compounds_ID'),('one_per_structure','Standardized_SMILES')]:
                if selected.groupby(key).Label.nunique().gt(1).any():
                    raise ValueError(f'Conflicting labels for {dataset}, {threshold}, {key}.')
                variants[unit]=selected.drop_duplicates(key)
            repeated=selected.Standardized_SMILES.eq(REPEATED_STRUCTURE)
            variants['repeated_hit_once_per_COADD_ID']=pd.concat([selected.loc[~repeated],selected.loc[repeated].drop_duplicates('COADD_Compounds_ID')])
            for unit,cohort in variants.items():
                key=(dataset,unit,threshold)
                cohorts[key]=cohort
                curves[key]={}
                for model in ['GI60','GI80']:
                    c=curve(cohort.Label,cohort[f'DeepFlux_{model}_mean'])
                    curves[key][model]=c
                    metrics.append({'Dataset':dataset,'Counting_Unit':unit,'Growth_Inhibition_Cutoff':threshold,'Model':model,'Positive_Records':int(cohort.Label.sum()),'Negative_Records':int(cohort.Label.eq(0).sum()),'AUC':c['auc']})
                    coordinates.extend({'Dataset':dataset,'Counting_Unit':unit,'Growth_Inhibition_Cutoff':threshold,'Model':model,'False_Positive_Rate':a,'True_Positive_Rate':b,'Score_Threshold':t} for a,b,t in zip(c['fpr'],c['tpr'],c['threshold']))
        key=('workbook','workbook_records',threshold)
        curves[key]={}
        labels=workbook[f'Binarized_at_{threshold}_PGI']
        for model,number in [('GI60',60),('GI80',80)]:
            c=curve(labels,workbook[f'GI_{number}_model_mean'])
            curves[key][model]=c
            metrics.append({'Dataset':'workbook','Counting_Unit':'workbook_records','Growth_Inhibition_Cutoff':threshold,'Model':model,'Positive_Records':int(labels.sum()),'Negative_Records':int(labels.eq(0).sum()),'AUC':c['auc']})
            coordinates.extend({'Dataset':'workbook','Counting_Unit':'workbook_records','Growth_Inhibition_Cutoff':threshold,'Model':model,'False_Positive_Rate':a,'True_Positive_Rate':b,'Score_Threshold':t} for a,b,t in zip(c['fpr'],c['tpr'],c['threshold']))
    metric_frame=pd.DataFrame(metrics)
    metric_frame.to_csv(DIRECTORY/'ROC_comparison_metrics.csv',index=False)
    metric_frame.loc[metric_frame.Model.eq('GI60')&metric_frame.Growth_Inhibition_Cutoff.eq(60)].to_csv(DIRECTORY/'GI60_at_60_percent_ROC_metrics.csv',index=False)
    pd.DataFrame(coordinates).to_csv(DIRECTORY/'ROC_comparison_plot_coordinates.csv',index=False)
    panel_cases=[
        (('workbook','workbook_records',60),'Earlier workbook','155 hits; 34,043 nonhits'),
        (('filtered','all_records',60),'Filtered CSV, all records','205 hits; 34,456 nonhits'),
        (('unfiltered','all_records',60),'Unfiltered CSV, all records','217 hits; 34,567 nonhits'),
        (('unfiltered','one_per_COADD_ID',60),'Unfiltered CSV, one per ID','148 hits; 32,275 nonhits'),
    ]
    panels=[]
    for key,title,subtitle in panel_cases:
        stem=DIRECTORY/('_'.join(map(str,key))+'_GI60_roc_2a')
        render_roc_panel(curves[key],stem,gi60_only=True)
        panels.append((stem.with_suffix('.pdf'),title,subtitle))
    strip(panels,DIRECTORY/'GI60_60_percent_cutoff_comparison_2a')
    panels=[]
    for key,title,subtitle in [
        (('workbook','workbook_records',80),'Earlier workbook, 80% cutoff','92 hits; 34,106 nonhits'),
        (('unfiltered','all_records',80),'Unfiltered CSV, all records','154 hits; 34,630 nonhits'),
        (('unfiltered','repeated_hit_once_per_COADD_ID',80),'Repeated hit counted once per ID','92 hits; 34,630 nonhits'),
    ]:
        stem=DIRECTORY/('_'.join(map(str,key))+'_roc_2a')
        render_roc_panel(curves[key],stem)
        panels.append((stem.with_suffix('.pdf'),title,subtitle))
    strip(panels,DIRECTORY/'workbook_and_unfiltered_80_percent_comparison_2a')
    for key in [('unfiltered','one_per_COADD_ID',60),('unfiltered','one_per_structure',60)]:
        cohorts[key].to_csv(DIRECTORY/('_'.join(map(str,key))+'_scored_ROC_cohort.csv'),index=False)
    definitions={'positive_class':'Delta5 growth inhibition >= selected cutoff','negative_class':'Delta5 and WT growth inhibition < selected cutoff','score':'original ensemble mean','cutoffs':[60,80],'row_units':['all_records','one_per_COADD_ID','one_per_structure','repeated_hit_once_per_COADD_ID'],'structure_key':'Standardized_SMILES','repeated_structure':REPEATED_STRUCTURE,'workbook_labels':['Binarized_at_60_PGI','Binarized_at_80_PGI']}
    (DIRECTORY/'comparison_definition.json').write_text(json.dumps(definitions,indent=2)+'\n')
    print(metric_frame.loc[metric_frame.Model.eq('GI60')&metric_frame.Growth_Inhibition_Cutoff.eq(60)].to_string(index=False))


if __name__=='__main__':
    main()
