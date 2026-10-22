from pathlib import Path
import pandas as pd,numpy as np,html,hashlib
from PIL import Image,ImageDraw,ImageFont
P=Path(__file__).parent
cv=pd.read_csv(P/'curve_data.csv');cn=pd.read_csv(P/'class_counts.csv')
W,H=2400,2180
im=Image.new('RGB',(W,H),'white');d=ImageDraw.Draw(im,'RGBA')
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><title>Five-variable trajectories across three cohorts</title><rect width="100%" height="100%" fill="white"/>']
ink='#243F53';colors={'C1':'#0D4483','C2':'#DD1215','C3':'#39A08F'}
def txt(x,y,t,size=26,bold=False,anchor='middle',color=ink):
 f=ImageFont.truetype('C:/Windows/Fonts/arialbd.ttf' if bold else 'C:/Windows/Fonts/arial.ttf',size)
 d.text((x,y),t,font=f,fill=color,anchor={'middle':'ms','start':'ls','end':'rs'}[anchor])
 svg.append(f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{700 if bold else 400}" fill="{color}" text-anchor="{anchor}">{html.escape(t)}</text>')
def rect(x,y,w,h,fill):
 d.rectangle((x,y,x+w,y+h),fill=fill);svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}"/>')
def line(pts,color,width=2):
 d.line(pts,fill=color,width=width);svg.append('<polyline points="'+' '.join(f'{x:.2f},{y:.2f}' for x,y in pts)+f'" fill="none" stroke="{color}" stroke-width="{width}" stroke-linejoin="round"/>')
def band(pts,color):
 rgb=tuple(int(color[k:k+2],16) for k in (1,3,5));d.polygon(pts,fill=rgb+(35,));svg.append('<polygon points="'+' '.join(f'{x:.2f},{y:.2f}' for x,y in pts)+f'" fill="{color}" fill-opacity="0.14"/>')
DB=['MIMIC','eICU','SICdb'];titles=['MIMIC-IV','eICU-CRD','SICdb']
for j,db in enumerate(DB):
 x=150+750*j
 txt(x+325,65,titles[j],31,True)
 a=cn[cn.dataset==db]
 font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',29)
 entries=[(r.phenotype,f'{r.phenotype}: {r.n:,}') for r in a.itertuples()]
 widths=[42+font.getlength(t) for cl,t in entries];total=sum(widths)+48
 pos=x+325-total/2
 for (cl,t),width in zip(entries,widths):
  line([(pos,102),(pos+30,102)],colors[cl],5)
  txt(pos+42,112,t,29,anchor='start')
  pos+=width+24
marks=[('MSI','MSI (bpm/mmHg)',[1.15,1.25,1.35,1.45]),('Na','Sodium (mmol/L)',[137,138,139,140,141]),('K','Potassium (mmol/L)',[3.8,4.0,4.2,4.4]),('Cl','Chloride (mmol/L)',[100,102.5,105,107.5]),('HCO3','Bicarbonate (mmol/L)',[17.5,20,22.5,25,27.5])]
for i,(m,label,ticks) in enumerate(marks):
 top=150+i*400
 rect(35,top,2330,370,'#F6F9FB' if i%2==1 else '#FFFFFF')
 txt(75,top+37,label,29,True,anchor='start')
 z=cv[cv.indicator==m];low=min(z.lower.min(),min(ticks));high=max(z.upper.max(),max(ticks));pad=(high-low)*.09;low-=pad;high+=pad
 first,last=(1,23) if m=='MSI' else (3,21)
 centers=list(range(2,24,4)) if m=='MSI' else [3,9,15,21]
 labels=[f'{h-2}–{h+2}' for h in centers] if m=='MSI' else [f'{h-3}–{h+3}' for h in centers]
 for j,db in enumerate(DB):
  x=150+750*j;y=top+65;ph=240;pw=650
  xx=lambda h:x+15+(np.asarray(h)-first)/(last-first)*(pw-30)
  yy=lambda v:y+ph-(np.asarray(v)-low)/(high-low)*ph
  for t in ticks:
   yt=float(yy(t));line([(x,yt),(x+pw,yt)],'#DDE5EB',1)
   if j==0:txt(x-18,yt+8,f'{t:g}',23,anchor='end')
  line([(x,y+ph),(x+pw,y+ph)],'#9DAFBC',1)
  for h,l in zip(centers,labels):
   xt=float(xx(h));line([(xt,y+ph),(xt,y+ph+5)],'#9DAFBC',1);txt(xt,y+ph+31,l,24)
  for cl,col in colors.items():
   a=z[(z.dataset==db)&(z.phenotype==cl)].sort_values('hour')
   assert len(a)==201 and a.hour.min()==first and a.hour.max()==last
   xp=xx(a.hour.to_numpy());lo=yy(a.lower.to_numpy());hi=yy(a.upper.to_numpy());mu=yy(a['mean'].to_numpy())
   band(list(zip(xp,lo))+list(zip(xp[::-1],hi[::-1])),col)
   line(list(zip(xp,mu)),col,5)
txt(1200,2150,'Time window after sepsis onset (h)',29)
svg.append('</svg>');(P/'Figure2.svg').write_text('\n'.join(svg),encoding='utf-8');im.save(P/'Figure2.png')
(P/'绘图核查.txt').write_text('45 curves, 201 points each. All original mean/lower/upper values retained. Row-specific midpoint ranges mapped linearly to equal panel widths; no extrapolation or refitting.\nSource SHA256: '+hashlib.sha256((P/'curve_data.csv').read_bytes()).hexdigest(),encoding='utf-8')
print('Created 5 x 3 figure with unchanged curve estimates.')
