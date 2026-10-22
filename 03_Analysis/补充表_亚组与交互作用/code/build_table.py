from pathlib import Path
import pandas as pd,numpy as np,html,json,hashlib
from reportlab.platypus import SimpleDocTemplate,Table,TableStyle,Paragraph,Spacer,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4,landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from docx import Document
from docx.shared import Pt,Mm
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
import pypdfium2 as pdfium
P=Path(__file__).resolve().parent.parent
DB=['MIMIC','eICU','SICdb'];SG=['Age','Sex','Diabetes']
for db in DB:
 assert f'COMPLETE {db}' in (P/'diagnostics'/f'{db}_safe_run.log').read_text(encoding='utf-8',errors='replace')
def combine(s):
 return pd.concat([pd.read_csv(P/f'{db}_{s}.csv').assign(dataset=db) for db in DB],ignore_index=True)
r=combine('RR');it=combine('interaction');counts=combine('observed_counts');ic=combine('imputed_counts');audit=combine('model_audit');draws=combine('draws');conv=combine('convergence')
assert len(r)==48 and len(it)==12 and len(audit)==330
assert audit.converged.all() and len(it[it.p_interaction.isna()])==1
# Publication selection revised by user: remove eGFR stratification; retain original fits for audit.
r=r[r.subgroup.isin(SG)].copy();it=it[it.subgroup.isin(SG)].reset_index(drop=True);counts=counts[counts.subgroup.isin(SG)].copy();ic=ic[ic.subgroup.isin(SG)].copy()
p=it.p_interaction.to_numpy();ix=np.where(np.isfinite(p))[0];order=ix[np.argsort(p[ix])];q=np.minimum.accumulate((p[order]*9/np.arange(1,len(order)+1))[::-1])[::-1];it['q_BH']=np.nan;it.loc[order,'q_BH']=np.minimum(q,1)
# Independently reproduce every scalar pooled estimate and observed count total.
for row in r[r.RR.notna()].itertuples():
 z=draws[(draws.dataset==row.dataset)&(draws.subgroup==row.subgroup)&(draws.stratum==row.stratum)&(draws['class']==row.comparison[:2])]
 assert len(z)==30 and abs(np.exp(z.Q.mean())-row.RR)<1e-10
for db,n,event in [('MIMIC',9758,2442),('eICU',5709,2060),('SICdb',3635,949)]:
 for sg in SG:
  z=counts[(counts.dataset==db)&(counts.subgroup==sg)];assert z.n.sum()==n and z.events.sum()==event
r.to_csv(P/'pooled_subgroup_RR.csv',index=False);it.to_csv(P/'pooled_interaction_FDR.csv',index=False);counts.to_csv(P/'observed_events_and_N.csv',index=False)
ic.groupby(['dataset','subgroup','stratum','class']).agg(n_min=('n','min'),n_max=('n','max'),events_min=('events','min'),events_max=('events','max')).to_csv(P/'imputed_counts_ranges.csv')
audit.to_csv(P/'model_fit_audit.csv',index=False);conv.to_csv(P/'imputation_convergence.csv',index=False)
fmtp=lambda p: '<0.001' if p<.001 else f'{p:.3f}'
labels={'Age':['<65 years','≥65 years'],'Sex':['Female','Male'],'Diabetes':['No','Yes']}
headers=['Subgroup','Stratum','C1\nEvents/N','C2\nEvents/N','C3\nEvents/N','C2 vs C1\nRR (95% CI)','C3 vs C1\nRR (95% CI)','P for\ninteraction','FDR q']
note='Data are presented as events/patients and risk ratios (RRs) with 95% confidence intervals (CIs), with C1 as the reference within each stratum. RRs were estimated using modified Poisson regression adjusted for Model 3 covariates. Interaction P values test the interaction between trajectory class and subgroup. FDR q values were calculated using the Benjamini–Hochberg procedure across nine interaction tests.'
pdfmetrics.registerFont(TTFont('Calibri','C:/Windows/Fonts/calibri.ttf'))
pdfmetrics.registerFont(TTFont('CalibriBold','C:/Windows/Fonts/calibrib.ttf'))
sty=ParagraphStyle('cell',fontName='Calibri',fontSize=10,leading=13,alignment=1)
bold=ParagraphStyle('bold',parent=sty,fontName='CalibriBold')
ns=ParagraphStyle('note',parent=sty,fontSize=9,leading=12,alignment=0)
title=ParagraphStyle('title',parent=bold,fontSize=12,leading=15,alignment=0)
subtitle_style=ParagraphStyle('subtitle',parent=bold,alignment=0)
widths=[68,82,76,76,76,115,115,79,67]
story=[]
htmlparts=['<html><meta charset="utf-8"><style>body{font-family:Calibri,Arial;margin:30px;color:black}section{page-break-after:always}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:12px 5px;text-align:center;vertical-align:middle}thead{border-top:2px solid;border-bottom:1px solid}tbody{border-bottom:2px solid}p{font-size:12px;line-height:1.4}</style>']
for page,db in enumerate(DB):
 if page:story.append(PageBreak())
 titletext='Supplementary Table Sx  Subgroup associations of trajectory classes with MAKE30'
 subtitle=f'{dict(MIMIC="MIMIC-IV",eICU="eICU-CRD",SICdb="SICdb")[db]}  |  Model 3 adjusted analyses'
 story.extend([Paragraph(titletext,title),Spacer(1,4),Paragraph(subtitle,subtitle_style),Spacer(1,12)])
 rows=[headers]
 for sg in SG:
  tt=it[(it.dataset==db)&(it.subgroup==sg)].iloc[0]
  for g in [0,1]:
   vals=[sg if g==0 else '',labels[sg][g]]
   for cl in ['C1','C2','C3']:
    a=counts[(counts.dataset==db)&(counts.subgroup==sg)&(counts.stratum.astype(str)==str(g))&(counts['class']==cl)].iloc[0]
    vals.append(f'{int(a.events):,}/{int(a.n):,}')
   for cl in ['C2','C3']:
    a=r[(r.dataset==db)&(r.subgroup==sg)&(r.stratum==g)&(r.comparison==f'{cl} vs C1')].iloc[0]
    vals.append(f'{a.RR:.2f} ({a.low:.2f}–{a.high:.2f})')
   vals.extend([fmtp(tt.p_interaction),fmtp(tt.q_BH)] if g==0 else ['',''])
   rows.append(vals)
 cells=[[Paragraph(html.escape(str(v)).replace('\n','<br/>'),bold if i==0 or j==0 else sty) for j,v in enumerate(row)] for i,row in enumerate(rows)]
 cmds=[('VALIGN',(0,0),(-1,-1),'MIDDLE'),('LEFTPADDING',(0,0),(-1,-1),3),('RIGHTPADDING',(0,0),(-1,-1),3),('LINEABOVE',(0,0),(-1,0),.8,'black'),('LINEBELOW',(0,0),(-1,0),.5,'black'),('LINEBELOW',(0,-1),(-1,-1),.8,'black')]
 for i in [1,3,5]:
  for col in [0,7,8]:cmds.append(('SPAN',(col,i),(col,i+1)))
 t=Table(cells,colWidths=widths,rowHeights=[39]+[33]*6,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle(cmds))
 story.extend([t,Spacer(1,9),Paragraph(html.escape(note),ns)])
 htmlparts.append(f'<section><h3>{titletext}</h3><b>{subtitle}</b><table><thead><tr>'+''.join(f'<th>{html.escape(v).replace(chr(10),"<br>")}</th>' for v in headers)+'</tr></thead><tbody>')
 for i,row in enumerate(rows[1:],1):
  htmlparts.append('<tr>')
  for j,v in enumerate(row):
   if j in [0,7,8] and i%2==0:continue
   span=' rowspan="2"' if j in [0,7,8] else ''
   htmlparts.append(f'<td{span}>{html.escape(v)}</td>')
  htmlparts.append('</tr>')
 htmlparts.append('</tbody></table><p>'+html.escape(note)+'</p></section>')
htmlparts.append('</html>')
(P/'Supplementary_Table_subgroups.html').write_text(''.join(htmlparts),encoding='utf-8')
SimpleDocTemplate(str(P/'Supplementary_Table_subgroups.pdf'),pagesize=landscape(A4),leftMargin=33,rightMargin=33,topMargin=25,bottomMargin=22).build(story)
pdf=pdfium.PdfDocument(str(P/'Supplementary_Table_subgroups.pdf'));assert len(pdf)==3
for i,db in enumerate(DB):
 page=pdf[i];im=page.render(scale=2).to_pil();im.save(P/f'预览_{db}.png')
 im.crop((0,0,im.width,750)).save(P/f'预览_{db}_紧凑.png')
 assert 'Unknown subgroup' not in page.get_textpage().get_text_range()
(P/'table_legend.txt').write_text(note+'\n',encoding='utf-8')
print('Updated three cohort tables with merged subgroup and P/q cells; validated 36 RR rows and nine tests.')
