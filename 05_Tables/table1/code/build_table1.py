from pathlib import Path
import pandas as pd, numpy as np, json, hashlib
import warnings
warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
ROOT=Path(__file__).resolve().parent.parent; OUT=ROOT/'derived'
SRC=ROOT/'data_and_sql'; DBS=['MIMIC','eICU','SICdb']
data={}; dictionaries={}; provenance=[]; hashes=[]
for i,db in enumerate(DBS,1):
 p=next(SRC.glob(f'{i:02d}_*.csv')); raw=pd.read_csv(p,low_memory=False)
 assert len(raw)==[9804,5839,3831][i-1] and raw.cohort_unit_id.is_unique
 data[db]=raw.copy(); dictionaries[db]=pd.read_csv(next((SRC/'附加核查_可选').glob(f'{i+3:02d}_*')),keep_default_na=False)
 hashes.append({'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
def get(db,key,domain='lab',factor=1,specimen=None):
 d=dictionaries[db]; s=d[(d.domain==domain)&((d.source_key.astype(str)==str(key))|d.source_key.astype(str).str.startswith(str(key)+':'))]
 if specimen is not None:s=s[s.specimen==specimen]
 if len(s)!=1:raise ValueError((db,key,domain,len(s)))
 col=s.iloc[0].column_name
 provenance.append({'dataset':db,'source_key':key,'column':col,'factor':factor,'domain':domain,'unit':s.iloc[0].native_unit,'specimen':s.iloc[0].specimen})
 return pd.to_numeric(data[db][col],errors='coerce')*factor
base=json.loads((ROOT/'templates/variable_layout.json').read_text(encoding='utf-8'))
meta=[]
for r in base:
 if r['variable'].endswith('_t1'):continue
 if r['variable'] in ('diabetes_without_complication','diabetes_with_complication'):continue
 if r['variable']=='diabetes':r['label']='Diabetes'
 if r['variable']=='saps3':r['label']='SAPS 3 score (SICdb)'
 if r['variable']=='charlson_comorbidity_index':r['label']='Charlson Comorbidity Index'
 meta.append((r['section'],r['variable'],r['label'],r['section'].startswith(('Comorbidities','Treatment','Outcomes')) or r['variable']=='male'))
for db,a in data.items():
 a['male']=a.sex.astype(str).str.lower().map({'m':1,'male':1,'1':1,'f':0,'female':0,'0':0,'other':0})
 a['msi_window_mean']=a[[f'msi_m{i}' for i in range(1,13)]].mean(axis=1)
 if db=='SICdb':a['weight_kg']=get(db,'weight_on_admission','admission_anthropometry',.001)
 a.loc[a.weight_kg==0,'weight_kg']=np.nan
 bmi_existing=pd.to_numeric(a.bmi,errors='coerce')
 bmi_calculated=a.weight_kg/(pd.to_numeric(a.height_cm,errors='coerce').replace(0,np.nan)/100)**2
 a['bmi']=bmi_existing.fillna(bmi_calculated)
 for c in [m[1] for m in meta]:
  if c not in a:a[c]=np.nan
 for c in ['make30','make30_death','make30_new_krt','make30_persistent_kidney_dysfunction']:
  a.loc[a.make30_evaluable_window!=1,c]=np.nan
cci_path=SRC/'Charlson_补充材料/eICU_CCI_ICD分项及总分.csv'
cci=pd.read_csv(cci_path)
assert cci.patientunitstayid.is_unique and set(cci.patientunitstayid)==set(data['eICU'].cohort_unit_id)
data['eICU']['charlson_comorbidity_index']=data['eICU'].cohort_unit_id.map(cci.set_index('patientunitstayid').charlson_comorbidity_index)
assert data['eICU'].charlson_comorbidity_index.notna().sum()==5338
hashes.append({'file':str(cci_path),'sha256':hashlib.sha256(cci_path.read_bytes()).hexdigest()})
provenance.append({'dataset':'eICU','variable':'charlson_comorbidity_index','source_file':str(cci_path),'method':'Quan ICD-9/10 mapping, original Charlson weights plus age; no usable ICD codes = missing'})
labs=[
 ('albumin','Albumin, g/dL',50862,'albumin',287,1),
 ('sodium','Sodium, mmol/L',50983,'sodium',469,1),
 ('potassium','Potassium, mmol/L',50971,'potassium',463,1),
 ('chloride','Chloride, mmol/L',50902,'chloride',450,1),
 ('bicarbonate','Bicarbonate, mmol/L',50882,'bicarbonate',456,1),
 ('creatinine','Creatinine, mg/dL',50912,'creatinine',367,1),
 ('bun','Blood urea nitrogen, mg/dL',51006,'BUN',355,28.0134/60.0553),
 ('glucose','Glucose, mg/dL',50931,'glucose',348,1),
 ('hemoglobin','Hemoglobin, g/dL',51222,'Hgb',289,1),
 ('wbc','White blood cells, 10⁹/L',51301,'WBC x 1000',301,1),
 ('platelets','Platelets, 10⁹/L',51265,'platelets x 1000',314,1),
 ('bilirubin','Total bilirubin, mg/dL',50885,'total bilirubin',333,1),
 ('alt','ALT, U/L',50861,'ALT (SGPT)',617,1),
 ('ast','AST, U/L',50878,'AST (SGOT)',616,1),
 ('magnesium','Magnesium, mg/dL',50960,'magnesium',468,2.4305),
 ('phosphate','Phosphate, mg/dL',50970,'phosphate',471,3.0973762),
 ('inr','INR',51237,'PT - INR',3128,1),
]
for v,label,m,e,s,f in labs:
 for db,key in zip(DBS,[m,e,s]):data[db][v]=get(db,key,factor=f if db=='SICdb' else 1)
 meta.insert(next(i for i,m in enumerate(meta) if m[0]=='Comorbidities'),('Laboratory Values',v,label,False))
for v,k in [('sodium',686),('potassium',685),('chloride',683),('bicarbonate',666),('glucose',656)]:
 data['SICdb'][v]=data['SICdb'][v].fillna(get('SICdb',k,specimen='2296:Arterial'))
# Preserve one measurement channel per variable: no unweighted averaging of channel means.
vitals=[('heart_rate','Heart rate, beats/min',220045,'heartrate',707),('sbp','Systolic blood pressure, mmHg',220050,'systemicsystolic',701),('dbp','Diastolic blood pressure, mmHg',220051,'systemicdiastolic',702),('map','Mean arterial pressure, mmHg',220052,'systemicmean',703),('spo2','SpO₂, %',220277,'sao2',710),('rr','Respiratory rate, breaths/min',220210,'respiration',None),('temperature','Temperature, °C',223762,'temperature',709)]
extra=[]
for v,label,m,e,s in vitals:
 for db,key in zip(DBS,[m,e,s]):data[db][v]=get(db,key,'vital') if key else np.nan
 extra.append(('Vital Signs',v,label,False))
# Use a documented priority where an alternative measurement channel exists.
fallback={'sbp':(220179,'noninvasivesystolic',704),'dbp':(220180,'noninvasivediastolic',705),'map':(220181,'noninvasivemean',706)}
for v,ks in fallback.items():
 for db,k in zip(DBS,ks):data[db][v]=data[db][v].fillna(get(db,k,'vital'))
data['MIMIC']['temperature']=data['MIMIC'].temperature.fillna((get('MIMIC',223761,'vital')-32)*5/9)
supp=SRC/'SICdb_呼吸频率与体温补提取'
sp=supp/'01_SICdb_呼吸频率及体温_独立查询.csv'
sd=supp/'02_SICdb_候选字段字典_独立查询.csv'
new=pd.read_csv(sp);dictionary=pd.read_csv(sd)
assert len(new)==3831 and new.case_id.is_unique
assert set(new.case_id)==set(data['SICdb'].cohort_unit_id)
new=new.set_index('case_id').loc[data['SICdb'].cohort_unit_id].reset_index()
assert np.array_equal(new.sepsis_onset_offset_sec.to_numpy(),data['SICdb'].onset_offset.to_numpy())
assert np.array_equal(new.assigned_class.to_numpy(),data['SICdb'].assigned_class.to_numpy())
for sid,label in [(719,'RespRate'),(2274,'RespRateEKG'),(709,'Temperature')]:
 dr=dictionary.loc[dictionary.source_id==sid].iloc[0]
 assert dr.source_label==label and dr.source_group=='SignalFloat'
 assert dr.native_unit==('°C' if sid==709 else '/min')
 first='temperature' if sid==709 else 'respiratory'
 col=f'{first}_{sid}_mean';counts=new[f'{first}_{sid}_n'].fillna(0)
 assert ((counts>0)==new[col].notna()).all()
 provenance.append({'dataset':'SICdb','variable':'temperature' if sid==709 else 'rr','source_key':sid,'column':col,'domain':'vital','unit':dr.native_unit,'source_file':str(sp),'priority':1 if sid in (709,719) else 2})
rr=new.respiratory_719_mean.fillna(new.respiratory_2274_mean)
assert rr.notna().sum()==3796
assert np.allclose(data['SICdb'].temperature,new.temperature_709_mean,equal_nan=True)
data['SICdb']['rr']=rr.to_numpy()
data['SICdb']['temperature']=new.temperature_709_mean.to_numpy()
new['rr_selected_mean']=rr
new['rr_source_id']=np.where(new.respiratory_719_mean.notna(),719,np.where(new.respiratory_2274_mean.notna(),2274,np.nan))
new.to_csv(OUT/'SICdb_呼吸频率来源核查.csv',index=False,encoding='utf-8-sig')
for p in (sp,sd):hashes.append({'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
meta[4:4]=extra
for v,label,db,key,domain in [('apsiii','APS III (MIMIC)', 'MIMIC','apsiii','apsiii_admission'),('aps_eicu','Acute physiology score (eICU)','eICU','acutephysiologyscore','apache_result_admission'),('apache_eicu','APACHE score (eICU)','eICU','apachescore','apache_result_admission')]:
 for name in DBS:data[name][v]=np.nan
 data[db][v]=get(db,key,domain).replace(-1,np.nan)
 pos=next(i for i,m in enumerate(meta) if m[0]=='Laboratory Values')
 meta.insert(pos,('Severity Scores',v,label,False))
for db in DBS:
 data[db]['aps_score']=data[db]['apsiii'] if db=='MIMIC' else data[db]['aps_eicu'] if db=='eICU' else np.nan
meta=[(s,'aps_score','APS score',b) if v=='apsiii' else (s,v,l,b) for s,v,l,b in meta if v!='aps_eicu']
for db in DBS:
 a=data[db]
 if db!='SICdb':
  for var,prefix in [('hypertension','diagnosis_hypertension'),('ischemic_heart_disease','diagnosis_ischemicheart')]:
   cols=[c for c in a if c.startswith(prefix)];assert len(cols)==1
   a[var]=a[cols[0]];provenance.append({'dataset':db,'column':cols[0],'variable':var})
 else:
  a['hypertension']=a['history_preconditionarthyper_0no1yes_f285459dec'];a['ischemic_heart_disease']=np.nan
i=next(i for i,m in enumerate(meta) if m[0]=='Comorbidities')
meta[i:i]=[('Comorbidities','hypertension','Hypertension',True),('Comorbidities','ischemic_heart_disease','Ischemic heart disease',True)]
rows=[];audit=[]
for section,v,label,binary in meta:
 r=dict(section=section,variable=v,label=label); values={}; structural={}
 for db,a in data.items():
  x=pd.to_numeric(a[v],errors='coerce'); missing=int(x.isna().sum());structural[db]=missing==len(a)
  if section=='Comorbidities':x=x.fillna(0)
  values[db]=x.dropna();n=len(values[db]);q=values[db].quantile([.25,.5,.75]).tolist() if n else [np.nan]*3
  if not n:display='—'
  elif binary:display=f'{int((x==1).sum())} ({100*(x==1).sum()/n:.1f})'
  else:
   dec=2 if v in ['albumin','creatinine','baseline_creatinine','potassium','magnesium','phosphate','inr','msi_m1','msi_window_mean'] else 1
   display=f'{q[1]:.{dec}f} [{q[0]:.{dec}f}, {q[2]:.{dec}f}]'
  r[db]=display;r['missing_'+db]=f'{missing} ({missing/len(a)*100:.1f}%)'
  if section in ['Severity Scores','Comorbidities'] and structural[db]:r[db]='N.A.';r['missing_'+db]='N.A.'
  audit.append(dict(dataset=db,variable=v,n=len(a),missing_before_coding=missing,n_analyzed=n,events=int((x==1).sum()) if binary else None,q1=q[0],median=q[1],q3=q[2]))
  a[v+'_analysis']=x
 for db in DBS[1:]:
  x,y=values['MIMIC'],values[db];s=np.nan
  if len(x)>1 and len(y)>1 and not structural[db] and not structural['MIMIC']:
   den=np.sqrt((x.mean()*(1-x.mean())+y.mean()*(1-y.mean()))/2) if binary else np.sqrt((x.var()+y.var())/2)
   s=abs(x.mean()-y.mean())/den if den else (0 if x.mean()==y.mean() else np.nan)
  r['smd_'+db]=f'{s:.3f}' if np.isfinite(s) else '—'
 if label != 'Immunosuppression':rows.append(r)
OUT.mkdir(exist_ok=True)
(OUT/'Table1_rows.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
pd.DataFrame(rows).to_csv(OUT/'Table1_display.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(audit).to_csv(OUT/'Table1_numeric_audit.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(provenance).to_csv(OUT/'variable_sources.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(hashes).to_csv(OUT/'input_SHA256.csv',index=False)
for db,a in data.items():a[['cohort_unit_id','assigned_class']+[m[1]+'_analysis' for m in meta]].to_csv(OUT/(db+'_Table1_analysis.csv'),index=False)
assert data['SICdb'].albumin.notna().sum()==105
assert data['SICdb'].weight_kg.notna().sum()==3820
assert data['SICdb'].bmi.notna().sum()==3817
print('Table rows',len(rows));print(pd.DataFrame(audit).query("variable in ['albumin','weight_kg','bmi']").to_string(index=False))
