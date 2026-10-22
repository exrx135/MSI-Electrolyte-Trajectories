from pathlib import Path
import pandas as pd,numpy as np,html
from PIL import Image,ImageDraw,ImageFont
P=Path(__file__).parent
W,H=2700,1140
im=Image.new('RGB',(W,H),'white');d=ImageDraw.Draw(im)
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><title>Kaplan–Meier curves across three cohorts</title><rect width="100%" height="100%" fill="white"/>']
ink='#243F53';pal={'C1':'#0D4483','C2':'#DD1215','C3':'#39A08F'}
def txt(x,y,t,size=27,bold=False,anchor='middle',color=ink):
 font=ImageFont.truetype('C:/Windows/Fonts/arialbd.ttf' if bold else 'C:/Windows/Fonts/arial.ttf',size)
 d.text((x,y),t,font=font,fill=color,anchor={'middle':'ms','start':'ls','end':'rs'}[anchor])
 svg.append(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{700 if bold else 400}" fill="{color}">{html.escape(t)}</text>')
def line(pts,color=ink,width=2):
 d.line(pts,fill=color,width=width);svg.append('<polyline points="'+' '.join(f'{x:.3f},{y:.3f}' for x,y in pts)+f'" fill="none" stroke="{color}" stroke-width="{width}" stroke-linejoin="round"/>')
def rect(x,y,w,h):
 d.rectangle((x,y,x+w,y+h),fill='white',outline='#A6B2BB',width=1)
 svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="white" stroke="#A6B2BB"/>')
tests=pd.read_csv(P/'KM_总体logrank检验.csv')
for j,(db,title) in enumerate([('MIMIC','MIMIC-IV'),('eICU','eICU-CRD'),('SICdb','SICdb')]):
 x=125+890*j;y=120;pw=730;ph=600
 xx=lambda v:x+10+np.asarray(v)/30*(pw-20)
 yy=lambda v:y+ph-(np.asarray(v)-.4)/.6*(ph-12)
 txt(x+pw/2,65,title,35,True)
 line([(x,y),(x,y+ph),(x+pw,y+ph)],width=2)
 for t in np.arange(.4,1.01,.1):
  yt=float(yy(t));line([(x-7,yt),(x,yt)]);txt(x-17,yt+9,f'{t:.1f}',26,anchor='end')
 for t in range(0,31,5):
  xt=float(xx(t));line([(xt,y+ph),(xt,y+ph+8)]);txt(xt,y+ph+39,str(t),26)
 cv=pd.read_csv(P/f'{db}_KM曲线及95CI.csv');ri=pd.read_csv(P/f'{db}_KM风险人数.csv')
 for cl,col in pal.items():
  a=cv[cv.phenotype==cl].sort_values('time',kind='stable');assert a.iloc[0].survival==1
  pts=[];last=1
  for row in a.itertuples():
   xt=float(xx(row.time));pts.append((xt,float(yy(last))));pts.append((xt,float(yy(row.survival))));last=row.survival
  line(pts,col,4)
 pv=float(tests.loc[tests.dataset==db,'p'].iloc[0]);label='Log-rank P < 0.001' if pv<.001 else f'Log-rank P = {pv:.3f}'
 txt(x+pw-15,y+53,label,28,anchor='end')
 legend_font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',26)
 legend_labels=[f"{cl} (n = {int(ri[(ri.phenotype==cl)&(ri.time==0)].n_risk.iloc[0]):,})" for cl in pal]
 legend_width=93+max(legend_font.getlength(t) for t in legend_labels)+20
 rect(x+55,y+ph-185,legend_width,150)
 for k,(cl,col) in enumerate(pal.items()):
  ly=y+ph-149+45*k;n=int(ri[(ri.phenotype==cl)&(ri.time==0)].n_risk.iloc[0])
  line([(x+78,ly-8),(x+127,ly-8)],col,5);txt(x+148,ly,f'{cl} (n = {n:,})',26,anchor='start',color='#000000')
 # Equal-width table columns, independent from the continuous plot axis.
 txt(x,815,'Number at risk',28,anchor='start',color='#000000')
 tx=np.linspace(x+90,x+pw-35,5);times=[0,5,15,20,25]
 txt(x,865,'Day',25,anchor='start',color='#000000')
 for xp,t in zip(tx,times):txt(float(xp),865,str(t),25,True,color='#000000')
 line([(x,880),(x+pw,880)],'#DDE5EB',1)
 for k,(cl,col) in enumerate(pal.items()):
  ty=927+50*k;txt(x,ty,cl,26,anchor='start',color='#000000')
  for xp,t in zip(tx,times):
   val=int(ri[(ri.phenotype==cl)&(ri.time==t)].n_risk.iloc[0]);txt(float(xp),ty,str(val),25,color='#000000')
 txt(x+pw/2,1090,'Days after the 24-h landmark',27)
# Rotated y-axis title in native SVG, matched in preview.
svg.append(f'<text x="35" y="420" transform="rotate(-90 35 420)" text-anchor="middle" font-family="Arial" font-size="29" fill="{ink}">Survival probability</text>')
font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',29)
layer=Image.new('RGBA',(350,45),(255,255,255,0));ld=ImageDraw.Draw(layer);ld.text((175,22),'Survival probability',font=font,fill=ink,anchor='mm');layer=layer.rotate(90,expand=True);im.paste(layer,(13,245),layer)
svg.append('</svg>')
(P/'Figure3.svg').write_text('\n'.join(svg),encoding='utf-8');im.save(P/'Figure3.png')
print('Reused saved KM estimates, log-rank tests and risk counts. No survival analysis rerun.')
