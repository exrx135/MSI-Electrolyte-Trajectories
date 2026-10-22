from pathlib import Path
import pandas as pd, numpy as np, html, hashlib, shutil
from PIL import Image, ImageDraw, ImageFont
import sys
ROOT=Path('Z:/Research/MSI-1')
OUT=ROOT/'results/候选补充材料_20261018'
P=OUT/'图'
SOURCE=ROOT/'results/figure5'
DB=sys.argv[1]
PREFIX={'eICU':'FigureS1_eICU','SICdb':'FigureS2_SICdb'}[DB]
files=['nodes.csv','flows.csv','counts.csv','thresholds.csv']
nd=pd.read_csv(SOURCE/'nodes.csv');ed=pd.read_csv(SOURCE/'flows.csv');cn=pd.read_csv(SOURCE/'counts.csv')
nd=nd[nd.dataset==DB];ed=ed[ed.dataset==DB];cn=cn[cn.dataset==DB].sort_values('phenotype')
assert len(cn)==3
W,H=2400,2300
im=Image.new('RGB',(W,H),'white');d=ImageDraw.Draw(im,'RGBA')
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><title>{DB} within-class state transitions</title><rect width="100%" height="100%" fill="white"/>']
ink='#243F53'
states=['Low','Reference','High','Missing'];colors=dict(zip(states,['#00468B','#277938','#DD1215','#D0D5DA']))
def txt(x,y,t,size=26,bold=False,anchor='middle',color=ink):
 f=ImageFont.truetype('C:/Windows/Fonts/arialbd.ttf' if bold else 'C:/Windows/Fonts/arial.ttf',size)
 d.text((x,y),t,font=f,fill=color,anchor={'middle':'ms','start':'ls','end':'rs'}[anchor])
 svg.append(f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{700 if bold else 400}" fill="{color}" text-anchor="{anchor}">{html.escape(t)}</text>')
def rect(x,y,w,h,fill):
 d.rectangle((x,y,x+w,y+h),fill=fill);svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}"/>')
def line(pts,color,width=1):
 d.line(pts,fill=color,width=width);svg.append('<polyline points="'+' '.join(f'{x:.3f},{y:.3f}' for x,y in pts)+f'" fill="none" stroke="{color}" stroke-width="{width}"/>')
def ribbon(pts,col):
 rgb=tuple(int(col[k:k+2],16) for k in (1,3,5))
 alpha=76 if col!='#D0D5DA' else 65
 d.polygon(pts,fill=rgb+(alpha,));svg.append('<polygon points="'+' '.join(f'{x:.3f},{y:.3f}' for x,y in pts)+f'" fill="{col}" fill-opacity="{alpha/255:.6f}"/>')
# Legend colors encode clinical states, not trajectory classes.
legend=[('Low','Low / lower MSI'),('Reference','Reference / middle MSI'),('High','High / upper MSI'),('Missing','Missing')]
font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',29)
widths=[46+font.getlength(t) for s,t in legend];x=(W-sum(widths)-3*60)/2
for (s,t),width in zip(legend,widths):
 rect(x,30,28,28,colors[s]);txt(x+43,56,t,29,anchor='start');x+=width+60
for j,r in enumerate(cn.itertuples()):
 txt(475+j*750,115,f'{r.phenotype}  |  n = {r.n:,} ({r.pct:.1f}%)',31,True)
marks=[('MSI','MSI (bpm/mmHg)'),('Na','Sodium (mmol/L)'),('K','Potassium (mmol/L)'),('Cl','Chloride (mmol/L)'),('HCO3','Bicarbonate (mmol/L)')]
checks=[]
for i,(m,label) in enumerate(marks):
 top=145+i*420;rect(30,top,2340,405,'#F6F9FB' if i%2 else '#FFFFFF')
 txt(68,top+37,label,29,True,anchor='start')
 for j,cl in enumerate(['C1','C2','C3']):
  x=150+j*750;y=top+62;pw=650;ph=280
  xx=lambda h:x+h/24*pw
  yy=lambda v:y+ph-v/100*ph
  z=nd[(nd.indicator==m)&(nd.phenotype==cl)];flows=ed[(ed.indicator==m)&(ed.phenotype==cl)]
  hours=sorted(z.time_hour.unique());hw=.60 if m=='MSI' else .9
  for tick in [0,50,100]:
   line([(x,yy(tick)),(x+pw,yy(tick))],'#DFE6EB')
   if j==0:txt(x-18,yy(tick)+8,f'{tick}%',23,anchor='end')
  for h,h2 in zip(hours[:-1],hours[1:]):
   a=z[z.time_hour==h].set_index('state').loc[states];b=z[z.time_hour==h2].set_index('state').loc[states]
   f=flows[flows.from_hour==h].set_index(['from_state','to_state'])
   assert f.n.sum()==int(a.class_n.iloc[0])
   src=a.ymin.to_numpy().copy();dst=b.ymin.to_numpy().copy()
   for si,s in enumerate(states):
    for ti,t in enumerate(states):
     row=f.loc[(s,t)];v=float(row.pct)
     assert abs(v-100*row.n/a.class_n.iloc[0])<1e-8
     if v:
      u=np.linspace(0,1,50);ease=u*u*(3-2*u);xs=xx(h+hw+(h2-h-2*hw)*u);lo=src[si]+(dst[ti]-src[si])*ease
      ribbon(list(zip(xs,yy(lo)))+list(zip(xs[::-1],yy((lo+v)[::-1]))),colors[s])
     src[si]+=v;dst[ti]+=v
   assert np.allclose(src,a.ymax,atol=1e-8) and np.allclose(dst,b.ymax,atol=1e-8)
  for h in hours:
   a=z[z.time_hour==h];assert a.n.sum()==a.class_n.iloc[0]
   for r in a.itertuples():
    left,right=xx(h-hw),xx(h+hw);yt,yb=yy(r.ymax),yy(r.ymin)
    rect(left,yt,right-left,yb-yt,colors[r.state]);line([(left,yt),(right,yt)],'white')
    if m!='MSI' and r.pct>=8:
     txt(xx(h),(yt+yb)/2+8,f'{r.pct:.0f}%',22,True,color=ink if r.state=='Missing' else 'white')
  centers=list(range(2,24,4)) if m=='MSI' else [3,9,15,21]
  half=2 if m=='MSI' else 3
  for h in centers:txt(xx(h),y+ph+35,f'{h-half}–{h+half}',24)
  checks.append(f'{m} {cl}: {len(hours)} windows; all nodes and adjacent flows reconcile to class n={int(z.class_n.iloc[0])}.')
txt(1200,2280,'Time window after sepsis onset (h)',29)
svg.append('</svg>');(P/(PREFIX+'.svg')).write_text('\n'.join(svg),encoding='utf-8');im.save(P/(PREFIX+'.png'),dpi=(300,300))
im.resize((1440,1380),Image.Resampling.LANCZOS).save(P/(PREFIX+'_预览.png'))
(OUT/'核查'/(PREFIX+'_绘图核查.txt')).write_text('\n'.join(checks)+'\nSaved V330 proportions and flows unchanged; no refitting, regrouping or imputation.\n'+ '\n'.join(f'{f}: {hashlib.sha256((SOURCE/f).read_bytes()).hexdigest()}' for f in files),encoding='utf-8')
print(PREFIX, 'all 15 panels passed flow conservation checks.')




