from pathlib import Path
import pandas as pd,numpy as np,re,sqlite3,json,hashlib
P=Path(__file__).parent
a=pd.read_csv(P/'01_eICU_Charlson诊断及病史_独立导出.csv')
base=pd.read_csv(P.parent/'02_eicu_独立查询.csv',low_memory=False)
assert len(a)==5839 and a.patientunitstayid.is_unique and set(a.patientunitstayid)==set(base.cohort_unit_id)
z=a.merge(base[['cohort_unit_id','patienthealthsystemstayid','age']],left_on='patientunitstayid',right_on='cohort_unit_id',validate='one_to_one',suffixes=('','_base'))
assert np.array_equal(z.age,z.age_base) and np.array_equal(z.patienthealthsystemstayid,z.patienthealthsystemstayid_base)
records=[];invalid=[]
for r in a.itertuples():
 for token in set(re.split(r'[,|;]',str(r.diagnosis_codes))):
  t=token.strip().upper()
  if t in ('','NAN'):continue
  if re.fullmatch(r'(?:[0-9]{3}|V[0-9]{2}|E[0-9]{3})(?:\.[0-9]{1,4})?',t):ver=9
  elif re.fullmatch(r'[A-Z][0-9]{2}(?:\.[A-Z0-9]{1,5})?',t):ver=10
  else:invalid.append({'patientunitstayid':r.patientunitstayid,'token':t});continue
  records.append((r.patientunitstayid,ver,t.replace('.','')))
db=sqlite3.connect(':memory:');db.create_function('GREATEST',-1,lambda *args:max([x for x in args if x is not None],default=None))
pd.DataFrame(records,columns=['hadm_id','icd_version','icd_code']).drop_duplicates().to_sql('diagnoses_icd',db,index=False)
pd.DataFrame({'subject_id':a.patientunitstayid,'hadm_id':a.patientunitstayid}).to_sql('admissions',db,index=False)
pd.DataFrame({'hadm_id':a.patientunitstayid,'age':a.age}).to_sql('age',db,index=False)
sql=(P/'reference_MIT_LCP_charlson.sql').read_text(encoding='utf-8')
for src,dst in [('physionet-data.mimiciv_hosp.diagnoses_icd','diagnoses_icd'),('physionet-data.mimiciv_hosp.admissions','admissions'),('physionet-data.mimiciv_derived.age','age')]:sql=sql.replace('`'+src+'`',dst)
# Known scoring cases verify age boundaries, duplicate-code handling and hierarchy.
test=sqlite3.connect(':memory:');test.create_function('GREATEST',-1,lambda *args:max([x for x in args if x is not None],default=None))
pd.DataFrame({'subject_id':[1,2,3,4,5],'hadm_id':[1,2,3,4,5]}).to_sql('admissions',test,index=False)
pd.DataFrame({'hadm_id':[1,2,3,4,5],'age':[49,50,65,80,60]}).to_sql('age',test,index=False)
pd.DataFrame([(1,10,'I219'),(1,9,'4109'),(2,10,'E119'),(2,10,'E112'),(3,10,'K703'),(3,10,'K721'),(4,10,'C349'),(4,10,'C780'),(5,10,'I10')],columns=['hadm_id','icd_version','icd_code']).to_sql('diagnoses_icd',test,index=False)
got=pd.read_sql_query(sql,test).set_index('hadm_id').charlson_comorbidity_index.to_dict()
assert got=={1:1,2:3,3:5,4:10,5:2},got
out=pd.read_sql_query(sql,db).rename(columns={'hadm_id':'patientunitstayid'}).drop(columns='subject_id')
valid=pd.DataFrame(records,columns=['patientunitstayid','version','code']).groupby('patientunitstayid').size()
out['n_parsed_code_tokens']=out.patientunitstayid.map(valid).fillna(0).astype(int)
out['cci_disease_only']=out.charlson_comorbidity_index-out.age_score
# No usable diagnostic codes is missing ascertainment, not an age-only CCI.
out.loc[out.n_parsed_code_tokens==0,['cci_disease_only','charlson_comorbidity_index']]=np.nan
assert out.patientunitstayid.is_unique and len(out)==5839
out.to_csv(P/'eICU_CCI_ICD分项及总分.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(invalid,columns=['patientunitstayid','token']).to_csv(P/'eICU_未识别诊断代码.csv',index=False,encoding='utf-8-sig')
summary={'n':len(out),'known':int(out.charlson_comorbidity_index.notna().sum()),'missing':int(out.charlson_comorbidity_index.isna().sum()),'q1_median_q3':out.charlson_comorbidity_index.quantile([.25,.5,.75]).tolist(),'unparsed_tokens':len(invalid),'unique_unparsed_tokens':sorted({r['token'] for r in invalid}),'definition':'MIMIC official Quan ICD-9/10 mappings and age-adjusted original Charlson weights; ICD records only; history text not used as automatic imputation','reference_sha256':hashlib.sha256((P/'reference_MIT_LCP_charlson.sql').read_bytes()).hexdigest()}
(P/'eICU_CCI核查.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False))
