from pathlib import Path
import pandas as pd
import numpy as np
import json, hashlib

P=Path(__file__).resolve().parent.parent
ROOT=P.parent.parent
T=P.parent/'table1'
DBS=['MIMIC','eICU','SICdb']
M2=['age','male','bmi','baseline_egfr_mdrd','hypertension','diabetes','ckd_or_renal_disease']
ADDED=['sepsis_sofa_score','spo2','rr','glucose','hemoglobin','wbc','platelets','magnesium']
# Plausibility rules fixed before model fitting, applied identically in all cohorts.
# Broad bounds do not winsorize unusual but clinically possible observations.
# eGFR is a derived covariate: >300 is treated as implausible for this adult
# regression analysis, NOT used to recompute frozen baseline Cr or MAKE30.
BOUNDS={'age':(18,120),'bmi':(5,150),'weight_kg':(5,500),
        'baseline_egfr_mdrd':(0,300),'sepsis_sofa_score':(0,24),
        'spo2':(0,100),'rr':(0,100),'glucose':(0,2500),
        'hemoglobin':(0,30),'wbc':(0,1000),'platelets':(0,5000),'magnesium':(0,30)}
EXCLUSIVE_ZERO={'baseline_egfr_mdrd','spo2','rr','glucose','hemoglobin','magnesium'}
old_paths={
 'MIMIC':('MIMIC_预后关联/V319_MIMIC_analysis_cohort_G3.csv','stay_id'),
 'eICU':('eICU_外部验证/V321_analysis_cohort_G3.csv','patientunitstayid'),
 'SICdb':('SICdb_外部验证/V322_analysis_cohort_G3.csv','case_id')}
flow=[];audit=[];removed=[];hashes=[];inventory=[]
for i,db in enumerate(DBS,1):
 source=T/'derived'/f'{db}_Table1_analysis.csv'
 rawpath=next((T/'data_and_sql').glob(f'{i:02d}_*.csv'))
 oldpath=ROOT/'STEP5_冻结分型与预后外部验证'/old_paths[db][0]
 for q in [source,rawpath,oldpath]:
  hashes.append({'dataset':db,'file':str(q),'sha256':hashlib.sha256(q.read_bytes()).hexdigest()})
 a=pd.read_csv(source).rename(columns=lambda c:c.removesuffix('_analysis'))
 assert a.cohort_unit_id.is_unique
 all_n=len(a);a=a[a.make30.notna()].copy()
 old=pd.read_csv(oldpath,low_memory=False)
 old=old[old.make30.notna()].copy()
 assert set(a.cohort_unit_id)==set(old[old_paths[db][1]])
 z=a[['cohort_unit_id','assigned_class','make30']].merge(old[[old_paths[db][1],'assigned_class','make30']],left_on='cohort_unit_id',right_on=old_paths[db][1],suffixes=('_new','_old'),validate='one_to_one')
 assert (z.assigned_class_new==z.assigned_class_old).all()
 assert (z.make30_new==z.make30_old).all()
 # Binary sex predictor: neither unknown nor other is silently coded female.
 sex=pd.read_csv(rawpath,usecols=['cohort_unit_id','sex']).set_index('cohort_unit_id').sex.astype(str).str.lower()
 a['male']=a.cohort_unit_id.map(sex.map({'m':1,'male':1,'f':0,'female':0}))
 for v,(lo,hi) in BOUNDS.items():
  s=pd.to_numeric(a[v],errors='coerce')
  bad=(s<lo)|(s>hi)|(~np.isfinite(s)&s.notna())
  if v in EXCLUSIVE_ZERO:bad |= s.eq(0)
  for idx in a.index[bad]:removed.append({'dataset':db,'cohort_unit_id':a.at[idx,'cohort_unit_id'],'variable':v,'original_value':s.at[idx],'reason':f'outside ({lo if v in EXCLUSIVE_ZERO else "["+str(lo)}, {hi}]'})
  a[v]=s.mask(bad)
  audit.append({'dataset':db,'variable':v,'n':len(a),'missing_before':int(s.isna().sum()),'values_set_missing':int(bad.sum()),'missing_after':int(a[v].isna().sum()),'missing_pct_after':100*a[v].isna().mean(),'lower':lo,'upper':hi,'zero_invalid':v in EXCLUSIVE_ZERO})
 for v in M2+ADDED+['weight_kg']:
  assert a[v].isna().mean()<.30,(db,v,a[v].isna().mean())
 for c in a.columns:
  if c not in ['cohort_unit_id','assigned_class']:
   inventory.append({'dataset':db,'variable':c,'missing_n':int(a[c].isna().sum()),'n':len(a),'missing_pct':100*a[c].isna().mean(),'used_in_model2':c in M2,'used_in_model3':c in M2+ADDED})
 out=a[['cohort_unit_id','assigned_class','make30']+M2+ADDED+['weight_kg']].sort_values('cohort_unit_id')
 out.to_csv(P/'data'/f'{db}_analysis_input.csv',index=False)
 flow.append({'dataset':db,'trajectory_n':all_n,'excluded_unknown_make30':all_n-len(a),'analysis_n':len(a),'events':int(a.make30.sum()),'model2_complete_case_n':int(a[M2].notna().all(axis=1).sum()),'model3_complete_case_n':int(a[M2+ADDED].notna().all(axis=1).sum()),'male_missing':int(a.male.isna().sum()),'frozen_ids_classes_outcomes_match':True})
pd.DataFrame(flow).to_csv(P/'cohort_flow.csv',index=False)
pd.DataFrame(audit).to_csv(P/'cleaning_and_missingness.csv',index=False)
pd.DataFrame(removed).to_csv(P/'diagnostics'/'values_set_missing.csv',index=False)
pd.DataFrame(inventory).to_csv(P/'variable_inventory.csv',index=False)
pd.DataFrame(hashes).to_csv(P/'input_sha256.csv',index=False)
(P/'model_specification.json').write_text(json.dumps({'model1':[],'model2':M2,'model3':M2+ADDED,'auxiliary_only':['weight_observed_centered','weight_missing'],'m':30,'iterations':30,'seed':20261018,'imputation':'within cohort MICE; PMM for continuous, logistic for binary sex; fixed class and observed outcome included as predictors','poisson':'log link, HC1 robust variance','pooling':'Rubin rules on log RR with Barnard-Rubin df','missing_rule':'strictly <30% in each cohort after cleaning','bounds':BOUNDS,'zero_invalid':sorted(EXCLUSIVE_ZERO)},ensure_ascii=False,indent=2),encoding='utf-8')
print(pd.DataFrame(flow).to_string(index=False))
print(pd.DataFrame(audit).query('values_set_missing>0').to_string(index=False))
