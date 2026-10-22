from pathlib import Path
import pandas as pd, numpy as np, json, hashlib, shutil
from docx import Document
from docx.shared import Pt,Mm,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT,WD_ROW_HEIGHT_RULE
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
ROOT=Path('Z:/Research/MSI-1');R=ROOT/'results';P=R/'候选补充材料_20261018';T=P/'表';D=P/'数据'
DB=['MIMIC','eICU','SICdb'];names={'MIMIC':'MIMIC-IV','eICU':'eICU-CRD','SICdb':'SICdb'}
manifest=[]
def read(src,name):
 src=Path(src);shutil.copy2(src,D/name);manifest.append({'source':str(src),'snapshot':name,'sha256':hashlib.sha256(src.read_bytes()).hexdigest()});return pd.read_csv(src)
def new():
 d=Document();s=d.sections[0];s.page_width=Mm(297);s.page_height=Mm(210);s.left_margin=s.right_margin=Mm(13);s.top_margin=s.bottom_margin=Mm(12)
 for st in d.styles:
  for b in list(st.element.iter(qn('w:pBdr'))):b.getparent().remove(b)
 for name in ['Normal','Title']:
  st=d.styles[name];st.font.name='Calibri';st.font.color.rgb=RGBColor(0,0,0);st.font.size=Pt(9 if name=='Normal' else 12);st.paragraph_format.space_after=Pt(0);st.paragraph_format.line_spacing=1
 return d
def title(d,num,text,sub=None):
 p=d.add_paragraph(f'Supplementary Table S{num}  {text}','Title');p.runs[0].bold=True;p.paragraph_format.space_after=Pt(5)
 if sub:
  p=d.add_paragraph(sub);p.runs[0].bold=True;p.paragraph_format.space_after=Pt(7)
def note(d,s):
 p=d.add_paragraph(s);p.paragraph_format.space_before=Pt(7)
 for run in p.runs:run.font.size=Pt(9)
def tbl(d,headers,rows,weights,font=10,height=28,left=(),boldrows=(),merges=()):
 tb=d.add_table(rows=1+len(rows),cols=len(headers));tb.autofit=False;tb.alignment=WD_TABLE_ALIGNMENT.CENTER
 widths=[x/sum(weights)*768 for x in weights]
 for col,w in zip(tb.columns,widths):col.width=Pt(w)
 for i,vals in enumerate([headers]+rows):
  row=tb.rows[i];row.height=Pt(34 if i==0 else height);row.height_rule=WD_ROW_HEIGHT_RULE.AT_LEAST
  row._tr.get_or_add_trPr().append(OxmlElement('w:cantSplit'))
  if i==0:row._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))
  for j,(cell,w,v) in enumerate(zip(row.cells,widths,vals)):
   cell.width=Pt(w);cell.text=str(v);cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
   for para in cell.paragraphs:
    para.alignment=WD_ALIGN_PARAGRAPH.LEFT if j in left else WD_ALIGN_PARAGRAPH.CENTER
    para.paragraph_format.space_after=Pt(2);para.paragraph_format.space_before=Pt(2)
    for run in para.runs:run.font.size=Pt(font);run.bold=(i==0 or i-1 in boldrows)
 for col,start,end in merges:
  value=tb.cell(start+1,col).text;cell=tb.cell(start+1,col).merge(tb.cell(end+1,col));cell.text=value;cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
  for para in cell.paragraphs:
   para.alignment=WD_ALIGN_PARAGRAPH.CENTER
   for run in para.runs:run.font.size=Pt(font);run.bold=True
 borders=OxmlElement('w:tblBorders')
 for side in ['top','bottom','left','right','insideH','insideV']:
  e=OxmlElement('w:'+side);e.set(qn('w:val'),'single' if side in ['top','bottom'] else 'nil');e.set(qn('w:sz'),'6');e.set(qn('w:color'),'000000');borders.append(e)
 tb._tbl.tblPr.append(borders)
 for cell in tb.rows[0].cells:
  b=OxmlElement('w:tcBorders');e=OxmlElement('w:bottom');e.set(qn('w:val'),'single');e.set(qn('w:sz'),'4');e.set(qn('w:color'),'000000');b.append(e);cell._tc.get_or_add_tcPr().append(b)
 return tb
# S1 definitions are transcribed from actual extraction and analysis, not idealized harmonization.
d=new();title(d,1,'Key study definitions and database implementations','A  Cohort and trajectory measurements')
common=lambda x:[x,x,x]
rows=[['Database and role','MIMIC-IV v3.1; development','eICU-CRD v2.0; external validation','SICdb v1.0.8; external validation'],
['Eligible patient']+common('Age ≥18 years; sepsis onset during ICU stay; first qualifying stay/case per patient.'),
['Suspected infection','Derived suspicion-of-infection record based on antibiotic and culture pairing.','Recognized antibiotic administration/start as the infection proxy; no requirement for two administrations.','Recognized antibiotic administration as the infection proxy.'],
['SOFA criterion','Derived SOFA score ≥2, paired within −24 to +12 h of the infection anchor.','Locally derived hourly SOFA ≥2 within −24 to +12 h of the antibiotic anchor.','Recorded SOFA score (reference 3139) ≥2 within −24 to +12 h of the antibiotic anchor.'],
['Sepsis onset and landmark']+common('Onset is the earlier of the paired infection and qualifying SOFA times. Landmark is onset +24 h.'),
['MSI trajectory','Hourly median HR / selected hourly median MAP; maximum valid hourly MSI in each 2-h window.','Hourly HR and selected MAP from native vital records; maximum valid hourly MSI in each 2-h window.','Native hourly aggregate HR / selected MAP; maximum valid hourly MSI in each 2-h window.'],
['Electrolyte trajectories']+common('Sodium, potassium, chloride and bicarbonate: median within each of four 6-h windows after onset.'),
['Exposure eligibility']+common('Remain observed in ICU through onset +24 h; no KRT during the first 24 h; ≥9/12 MSI windows and ≥2/4 windows for each electrolyte.')]
tbl(d,['Definition','MIMIC-IV','eICU-CRD','SICdb'],rows,[100,222,222,222],font=9,height=35,left=(0,1,2,3))
note(d,'Baseline SOFA was assumed to be zero in the extraction rules. Infection proxies and native measurement aggregation differ by database. MSI, modified shock index (heart rate/mean arterial pressure); KRT, kidney replacement therapy.')
d.add_page_break();title(d,1,'Key study definitions and database implementations','B  MAKE30 definitions retained in the analysis')
rows=[['Risk window']+common('After sepsis onset +24 h, through the earlier of hospital discharge and 30 days after index ICU admission.'),
['Baseline creatinine','Median available creatinine in the 365 days before admission, or MDRD back-calculation assuming eGFR 75 if absent; take the lower of this value and first admission creatinine.','Lower of first admission creatinine and MDRD back-calculation assuming eGFR 75; a comparable pre-admission laboratory timeline was unavailable.','Median available creatinine in the preceding 365 days using linked patient/case offsets, or MDRD back-calculation assuming eGFR 75 if absent; take the lower of this value and first admission creatinine.'],
['Baseline eGFR','MDRD estimate from the constructed baseline creatinine.','MDRD estimate from the constructed baseline creatinine; this construction imposes an eGFR floor of 75.','MDRD estimate from the constructed baseline creatinine.'],
['MAKE30 composite']+common('Any of death, new KRT or persistent kidney dysfunction within the defined window. A positive component establishes an event even when another component is missing.'),
['Death component','Hospital death flag with death time in the window; discharge time substitutes when death time is absent.','Hospital discharge status Expired, with hospital discharge offset in the window; the frozen SQL assigns other statuses to zero.','Recorded death offset in the window, subject to the frozen discharge/administrative cutoff.'],
['New KRT']+common('First detected KRT after the landmark and by the cutoff, with no detected KRT at or before the landmark in the available linked records.'),
['Persistent kidney dysfunction']+common('Last available inpatient creatinine after the landmark and by the cutoff ≥2 times baseline creatinine. Unknown when the required creatinine is unavailable.'),
['Discharge cutoff','Recorded hospital discharge timestamp.','Recorded hospital discharge offset, aligned to the hospital timeline.','Frozen cutoff uses HospitalStayDays ×86400, otherwise HospitalDischargeDay ×86400, otherwise TimeOfStay; the administrative cutoff is ICUOffset +30 days.']]
tbl(d,['Definition','MIMIC-IV','eICU-CRD','SICdb'],rows,[100,222,222,222],font=8.5,height=32,left=(0,1,2,3))
note(d,'The 30-day administrative limit is indexed to ICU admission, not to the landmark. SICdb discharge fields retain day-level precision, origin and fallback limitations. These MAKE30 definitions were not revised using the later mortality adjudication for Figure 3. eGFR is expressed in mL/min/1.73 m².')
d.add_page_break();title(d,1,'Key study definitions and database implementations','C  Mortality analysis and covariate summaries')
rows=[['KM origin and censoring']+common('Time since sepsis onset +24 h. In-hospital death is the event; live hospital discharge is censored. Administrative limit: index ICU admission +30 days. No post-discharge mortality follow-up.'),
['Mortality ascertainment','Recorded hospital death time; discharge-time fallback retained for missing death times. Records with unresolved timing conflicts are excluded.','Hospital discharge status and supplementary discharge-location checks; hospital discharge offset serves as the death-time proxy.','Supplementary discharge-state/death-offset adjudication; deaths with confirmed in-window death timing are not truncated by the earlier day-level discharge proxy. Retained discharge-time limitations remain.'],
['KM sample size','9,706','5,804','3,811'],
['MAKE30 analysis sample size','9,758','5,709','3,635'],
['Trajectory / Table 1 sample size','9,804','5,839','3,831'],
['Dynamic baseline covariates']+common('Available measurements in the first 24 h after sepsis onset are summarized by their mean. Precomputed severity scores and baseline renal measures retain their source-specific definitions.'),
['BMI and stable characteristics']+common('Use available direct BMI; otherwise weight (kg)/height (m)². Demographics, height/weight and comorbidities use available admission/hospitalization information.'),
['Comorbidity coding']+common('Recorded diagnoses are mapped to study variables. Missing or unrecorded comorbidity indicators were coded as absent by the study rule; zero does not establish clinically verified absence.'),
['Covariate missingness handling']+common('Table 1 uses observed data. Adjusted regression uses separate cleaned analysis copies and multiple imputation; these do not overwrite Table 1 or frozen trajectory/outcome data.')]
tbl(d,['Definition','MIMIC-IV','eICU-CRD','SICdb'],rows,[100,222,222,222],font=8.8,height=29,left=(0,1,2,3))
note(d,'KM, Kaplan–Meier. The mortality analysis and MAKE30 use separately established evaluability rules and death ascertainment; their death counts and denominators are not interchangeable. Details of adjustment and imputation belong in the statistical methods.')
d.save(T/'TableS1_Definitions.docx')
# S2 initial screen and final classification
x=read(R/'补充表/V325_TableS5_G2G20_structure_G2G20.csv','S2_initial_screen.csv')
quality=read(R/'补充表/V325_TableS3_classification_quality_G3.csv','S2_final_classification.csv')
fit=read(ROOT/'STEP4_多起点精修与稳定性验证/多起点结果/V295_task_summary_G3_S01.csv','S2_final_fit.csv')
ms=read(ROOT/'STEP4_多起点精修与稳定性验证/多起点结果/V295_multistart_summary_G3.csv','S2_multistart.csv')
assert len(x)==19 and set(x.requested_g)==set(range(2,21))
d=new();title(d,2,'Trajectory model selection and classification quality','A  Initial MIMIC-IV candidate screen')
rows=[]
for a in x.itertuples():rows.append([int(a.requested_g),f'{a.final_logLik:,.1f}',f'{a.BIC:,.1f}',f'{a.min_class_pct:.2f}',f'{a.min_APPA:.3f}',int(a.n_class_lt_5pct),'Yes' if a.package_converged==1 else 'No'])
tbl(d,['Classes','Log likelihood','BIC','Smallest class (%)','Lowest class APPA','Classes <5%','Converged'],rows,[65,130,130,120,120,100,100],font=9,height=18,boldrows=(1,))
note(d,'BIC, Bayesian information criterion; APPA, average posterior probability of assigned class membership. The three-class row identifies the selected model family. BIC continued to decrease with additional classes; the three-class solution was not the BIC minimum. These initial-screen values must not be substituted for the final multistart-refined fit.')
d.add_page_break();title(d,2,'Trajectory model selection and classification quality','B  Final frozen three-class model and transported classifications')
rows=[]
for db in DB:
 z=quality[quality.dataset==db].sort_values('class');assert len(z)==3
 for a in z.itertuples():rows.append([names[db],a.phenotype,f'{a.n:,} ({a.pct:.1f})',f'{a.APPA:.3f}',f'{a.pct_max_posterior_ge_070:.1f}',f'{a.pct_max_posterior_ge_090:.1f}'])
tbl(d,['Cohort','Class','Patients, n (%)','APPA','Posterior ≥0.70 (%)','Posterior ≥0.90 (%)'],rows,[125,60,150,110,155,155],font=10,height=27,merges=[(0,0,2),(0,3,5),(0,6,8)])
f=fit.iloc[0]
note(d,f'The final MIMIC-IV G3 S01 fit had log likelihood {f.final_logLik:,.1f}, BIC {f.BIC:,.1f}, {int(f.npar)} parameters and {int(f.n_patients):,} patients. eICU and SICdb classifications use transported parameters without refitting. Posterior probabilities describe model assignment certainty, not externally verified classification accuracy.')
rows=[]
for a in ms.itertuples():rows.append([f'S{int(a.start_id):02d}'+(' (frozen)' if a.start_id==1 else ''),f'{a.BIC:,.1f}',f'{100*a.assignment_agreement:.2f}',f'{a.adjusted_rand_index:.3f}'])
d.add_paragraph('C  Five starts for the three-class model').runs[0].bold=True
tbl(d,['Start','BIC','Agreement with S01 (%)','Adjusted Rand index'],rows,[150,180,220,210],font=9,height=20)
note(d,'All five starts converged. Agreement statistics compare assignments after label alignment with the frozen S01 reference. Start S04 reached a different partition; convergence alone does not establish identical class membership across starts.')
d.save(T/'TableS2_Model_selection.docx')
# S3 frozen components on the exact outcome-analysis cohort; no outcome redefinition
paths={'MIMIC':ROOT/'STEP5_冻结分型与预后外部验证/MIMIC_预后关联/V319_MIMIC_analysis_cohort_G3.csv','eICU':ROOT/'STEP5_冻结分型与预后外部验证/eICU_外部验证/V321_analysis_cohort_G3.csv','SICdb':ROOT/'STEP5_冻结分型与预后外部验证/SICdb_外部验证/V322_analysis_cohort_G3.csv'}
endpoint={'make30':'MAKE30','make30_death':'Death','make30_new_krt':'New KRT','make30_persistent_kidney_dysfunction':'Persistent kidney dysfunction'}
rows=[];agg=[]
for db in DB:
 a=pd.read_csv(paths[db]);assert (len(a),int(a.make30.sum()))=={'MIMIC':(9758,2442),'eICU':(5709,2060),'SICdb':(3635,949)}[db]
 cols=list(endpoint)[1:];assert np.array_equal(a[cols].eq(1).any(axis=1).astype(int),a.make30.astype(int))
 manifest.append({'source':str(paths[db]),'snapshot':'aggregated only; no patient data copied','sha256':hashlib.sha256(paths[db].read_bytes()).hexdigest()})
 for v,label in endpoint.items():
  vals=[];missing=[]
  for cl in ['C1','C2','C3']:
   z=a.loc[a.phenotype==cl,v];n=int(z.notna().sum());events=int(z.sum());miss=int(z.isna().sum());vals.append(f'{events:,}/{n:,} ({events/n*100:.1f})');missing.append(str(miss));agg.append([db,cl,v,len(z),n,events,miss,events/n*100])
  rows.append([names[db],label]+vals+[' / '.join(missing)])
pd.DataFrame(agg,columns=['dataset','class','endpoint','cohort_n','known_n','events','missing_n','percent']).to_csv(D/'S3_components.csv',index=False)
d=new();title(d,3,'MAKE30 and component outcomes by trajectory class','Patients included in the MAKE30 association analyses')
tbl(d,['Cohort','Outcome','C1\nEvents/N (%)','C2\nEvents/N (%)','C3\nEvents/N (%)','Unknown, n\nC1 / C2 / C3'],rows,[95,170,135,135,135,110],font=9.5,height=30,left=(1,),merges=[(0,0,3),(0,4,7),(0,8,11)])
note(d,'N is the number with an evaluable value for the indicated component within the MAKE30 analysis cohort. Components may overlap and must not be summed. The frozen MAKE30 death component is shown; it is not the separately adjudicated mortality endpoint used for the KM curves in Figure 3. KRT, kidney replacement therapy. Persistent kidney dysfunction denotes last inpatient creatinine ≥2 times baseline within the study window.')
d.save(T/'TableS3_MAKE30_components.docx')
# S4/S5 approved external subgroup tables, number change only
for no,db in [(4,'eICU'),(5,'SICdb')]:
 src=R/f'补充表_亚组与交互作用/Supplementary_Table_{db}.docx';d=Document(src)
 for p in d.paragraphs:
  for run in p.runs:run.text=run.text.replace('Supplementary Table Sx',f'Supplementary Table S{no}')
 d.save(T/f'TableS{no}_Subgroups_{db}.docx')
 manifest.append({'source':str(src),'snapshot':f'TableS{no}; title number changed only','sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
# S6 all saved scenarios, classification and archived prognostic estimates distinguished from Model 3
z=read(R/'V331_MAR_MNAR敏感性_G3/V331_分类汇总_G3.csv','S6_classification.csv')
rr=read(R/'V331_MAR_MNAR敏感性_G3/V331_RR合并结果_G3.csv','S6_archived_RR.csv')
shifts=read(R/'V331_MAR_MNAR敏感性_G3/V331_情景标准差偏移_G3.csv','S6_delta_scenarios.csv')
scenarios=[('Original','Original classification'),('MAR','MAR'),('MNAR_down025','All −0.25 SD'),('MNAR_down050','All −0.50 SD'),('MNAR_up025','All +0.25 SD'),('MNAR_up050','All +0.50 SD'),('MNAR_pattern025','Pattern +0.25 SD'),('MNAR_pattern050','Pattern +0.50 SD'),('MNAR_reverse025','Pattern −0.25 SD'),('MNAR_reverse050','Pattern −0.50 SD')]
assert len(z)==30 and len(rr)==60
outrows=[];d=new()
for page,db in enumerate(DB):
 if page:d.add_page_break()
 title(d,6,'Sensitivity to missing trajectory measurements',names[db]+' | Fixed classifier; archived eight-covariate association model')
 rows=[]
 for sn,label in scenarios:
  a=z[(z.dataset==db)&(z.scenario==sn)].iloc[0]
  vals=[label,f'{100*a.agreement_mean:.1f} ({100*a.agreement_min:.1f}–{100*a.agreement_max:.1f})',f'{a.ARI_mean:.3f}',f'{a.n_C1_mean:.0f} / {a.n_C2_mean:.0f} / {a.n_C3_mean:.0f}']
  for co in ['C2 vs C1','C3 vs C1']:
   b=rr[(rr.dataset==db)&(rr.scenario==sn)&(rr.comparison==co)].iloc[0];vals.append(f'{b.RR:.2f} ({b.lower:.2f}–{b.upper:.2f})')
  rows.append(vals);outrows.append([db,sn]+vals)
 tbl(d,['Scenario','Same-class agreement\nMean (range), %','Mean ARI','Mean class counts\nC1 / C2 / C3','C2 vs C1\nRR (95% CI)','C3 vs C1\nRR (95% CI)'],rows,[145,160,70,165,125,125],font=9,height=25,boldrows=(0,1))
 note(d,'Classification summaries are means across 20 completed datasets; class counts are rounded means and ranges are across imputations, not confidence intervals. Only originally missing trajectory measurements were imputed or shifted; the classifier and observed values were fixed. MAR uses PMM imputation without a shift. MNAR shifts use frozen MIMIC-IV standard deviations; the pattern is MSI+, sodium unchanged, potassium+, chloride+, bicarbonate−, with the reverse pattern for negative shifts. Values beyond observed cohort ranges were clipped.')
 note(d,'RRs use the archived adjustment set: age, sex, sepsis-onset SOFA, baseline eGFR, chronic kidney/renal disease, diabetes, mechanical ventilation and vasopressor use. This is not the current Model 3 and must not be presented as its sensitivity analysis. The Original row uses this same archived specification. ARI, adjusted Rand index; PMM, predictive mean matching; MAR/MNAR, missing at random/not at random. These are assumed scenarios, not established missingness mechanisms.')
d.save(T/'TableS6_Missing_measurement_sensitivity.docx')
pd.DataFrame(outrows,columns=['dataset','scenario','label','agreement','ARI','class_counts','RR_C2','RR_C3']).to_csv(D/'S6_display_values.csv',index=False)
(P/'核查/source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('Created S1–S6 Word tables. S3 OR/denominator checks passed; S6 uses archived specification explicitly.')
