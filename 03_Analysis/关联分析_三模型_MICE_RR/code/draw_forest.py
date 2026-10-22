from pathlib import Path
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
import math, html

P=Path(__file__).resolve().parent.parent
a=pd.read_csv(P/'pooled_RR_results.csv')
assert len(a)==18
W,H,S=3360,820,2
im=Image.new('RGB',(W*S,H*S),'white');d=ImageDraw.Draw(im)
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><title>MAKE30 risk ratios across three adjustment models</title><rect width="100%" height="100%" fill="white"/>']
ink='#111111';pal={'MIMIC':'#0D4483','eICU':'#DD1215','SICdb':'#39A08F'}
titles={'MIMIC':'MIMIC-IV','eICU':'eICU-CRD','SICdb':'SICdb'}
def text(x,y,t,size=29,bold=False,anchor='middle',color=ink):
 f=ImageFont.truetype('C:/Windows/Fonts/arialbd.ttf' if bold else 'C:/Windows/Fonts/arial.ttf',size*S)
 d.text((x*S,y*S),t,font=f,fill=color,anchor={'middle':'ms','start':'ls','end':'rs'}[anchor])
 svg.append(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{700 if bold else 400}" fill="{color}">{html.escape(t)}</text>')
def line(pts,color=ink,width=2):
 d.line([(x*S,y*S) for x,y in pts],fill=color,width=width*S)
 svg.append('<polyline points="'+' '.join(f'{x:.3f},{y:.3f}' for x,y in pts)+f'" fill="none" stroke="{color}" stroke-width="{width}"/>')
def circle(x,y,r,color):
 d.ellipse(((x-r)*S,(y-r)*S,(x+r)*S,(y+r)*S),fill=color)
 svg.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{color}"/>')
def marker(x,y,db,color):
 if db=='MIMIC':circle(x,y,7,color);return
 pts=([(x,y-9),(x+9,y),(x,y+9),(x-9,y)] if db=='eICU'
      else [(x,y-11),(x+10,y+5.5),(x-10,y+5.5)])
 d.polygon([(px*S,py*S) for px,py in pts],fill=color)
 svg.append('<polygon points="'+' '.join(f'{px},{py}' for px,py in pts)+f'" fill="{color}"/>')
# Alternating comparison bands; these are visual guides, not statistical groups.
d.rectangle((30*S,245*S,3330*S,425*S),fill='#FAFAFA')
svg.append('<rect x="30" y="245" width="3300" height="180" fill="#FAFAFA"/>')
for j,(db,col) in enumerate(pal.items()):
 x=1080+j*400;line([(x,60),(x+52,60)],col,4);marker(x+26,60,db,col);text(x+70,71,titles[db],32,anchor='start')
lower=min(.9,float(a.CI95_low.min())*.94)
upper=max(4.5,math.ceil(float(a.CI95_high.max())*1.02*2)/2)
ticks=[v for v in [.5,.75,1,1.5,2,3,4,4.5,5,6,8,10,15,20] if lower<=v<=upper]
group_ys=[[300,335,370],[535,570,605]]
for j,model in enumerate(['Model 1','Model 2','Model 3']):
 x=230+j*1030;pw=650
 xx=lambda v:x+math.log(v/lower)/math.log(upper/lower)*pw
 text(x+pw/2,155,model,39,True)
 text(x+pw+185,215,'RR (95% CI)',29,True)
 for y in range(245,660,18):line([(xx(1),y),(xx(1),min(y+9,660))],'#AAB5BC',2)
 for comparison,ys in zip(['C2 vs C1','C3 vs C1'],group_ys):
  if j==0:text(x-35,ys[1]+10,comparison,32,True,anchor='end')
  for (db,col),y in zip(pal.items(),ys):
   z=a[(a.dataset==db)&(a.model==model)&(a.comparison==comparison)]
   assert len(z)==1
   r=z.iloc[0];lo,hi,mid=xx(r.CI95_low),xx(r.CI95_high),xx(r.RR)
   line([(lo,y),(hi,y)],col,4)
   for xp in [lo,hi]:line([(xp,y-6),(xp,y+6)],col,2)
   marker(mid,y,db,col)
   text(x+pw+185,y+10,f'{r.RR:.2f} ({r.CI95_low:.2f}–{r.CI95_high:.2f})',29)
 # Joined left and bottom axes form a clear right angle in each panel.
 line([(x,245),(x,660),(x+pw,660)],'#555555',2)
 for ys in group_ys:line([(x-7,ys[1]),(x,ys[1])],'#555555',2)
 for t in ticks:
  xp=xx(t);line([(xp,660),(xp,669)],ink,2);text(xp,706,f'{t:g}',26)
 text(x+pw/2,765,'Risk ratio (95% CI)',30)
svg.append('</svg>')
(P/'Fig3_三模型Poisson_RR.svg').write_text('\n'.join(svg),encoding='utf-8')
im.save(P/'Fig3_三模型Poisson_RR.png',dpi=(300,300))
# Compact preview uses the same artwork; full-resolution PNG remains publication output.
im.resize((W,H),Image.Resampling.LANCZOS).save(P/'Fig3_三模型Poisson_RR_预览.png')
print('Created native editable SVG and 6720 x 1640 PNG; shared logarithmic RR axes, compact comparison groups, and database-specific markers.')
