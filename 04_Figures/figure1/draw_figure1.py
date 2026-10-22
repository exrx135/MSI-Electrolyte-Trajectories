from pathlib import Path
import html,json,csv
from PIL import Image,ImageDraw,ImageFont
P=Path(__file__).parent
W,H=2400,1370
im=Image.new('RGB',(W,H),'white');draw=ImageDraw.Draw(im)
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><title>Selection of the three study cohorts</title><rect width="100%" height="100%" fill="white"/>']
ink='#314958';blue='#214F6B';light='#E4EFF6';red='#F7E3E1';green='#DDEEE4';purple='#EBE6F2'
def line(points,color=ink,width=2,arrow=False):
 draw.line(points,fill=color,width=width)
 svg.append('<polyline points="'+' '.join(f'{x},{y}' for x,y in points)+f'" fill="none" stroke="{color}" stroke-width="{width}"/>')
 if arrow:
  x,y=points[-1];px,py=points[-2]
  pts=[(x,y),(x-6,y-10),(x+6,y-10)] if x==px else [(x,y),(x-10,y-6),(x-10,y+6)]
  draw.polygon(pts,fill=color);svg.append('<polygon points="'+' '.join(f'{a},{b}' for a,b in pts)+f'" fill="{color}"/>')
def text(x,y,lines,size=21,bold=False,color=ink,anchor='middle'):
 font=ImageFont.truetype('C:/Windows/Fonts/arialbd.ttf' if bold else 'C:/Windows/Fonts/arial.ttf',size)
 for i,t in enumerate(lines):
  yy=y+i*(size+8)
  svg.append(f'<text x="{x}" y="{yy}" text-anchor="{anchor}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{700 if bold else 400}" fill="{color}">{html.escape(t)}</text>')
  draw.text((x,yy),t,font=font,fill=color,anchor='ms' if anchor=='middle' else 'ls')
def box(x,y,w,h,lines,fill=light,size=21,bold=False,color=ink):
 draw.rounded_rectangle((x,y,x+w,y+h),radius=7,fill=fill,outline=ink,width=2)
 svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="7" fill="{fill}" stroke="{ink}" stroke-width="2"/>')
 block=len(lines)*(size+8)-8
 text(x+w/2,y+(h-block)/2+size-3,lines,size,bold,color)
D=[dict(db='MIMIC-IV v3.1',role='Trajectory derivation',start=65366,adult=16788,land=15354,pool=14529,sql=9804,final=9804,make=9758,km=9706,exclude=[48578,1434,825,4725]),dict(db='eICU-CRD v2.0',role='External application',start=139367,adult=24274,land=19759,pool=18824,sql=5867,final=5839,make=5709,km=5804,exclude=[115093,4515,935,12957]),dict(db='SICdb v1.0.8',role='External application',start=21566,adult=5344,land=4393,pool=4139,sql=3831,final=3831,make=3635,km=3811,exclude=[16222,951,254,308])]
def curve(start,control1,control2,end):
 color='#8561A2'
 svg.append(f'<path d="M {start[0]} {start[1]} C {control1[0]} {control1[1]}, {control2[0]} {control2[1]}, {end[0]} {end[1]}" fill="none" stroke="{color}" stroke-width="3"/>')
 pts=[]
 for j in range(101):
  t=j/100;u=1-t
  pts.append((u**3*start[0]+3*u*u*t*control1[0]+3*u*t*t*control2[0]+t**3*end[0],u**3*start[1]+3*u*u*t*control1[1]+3*u*t*t*control2[1]+t**3*end[1]))
 draw.line(pts,fill=color,width=3)
 x,y=end
 import math
 dx=x-control2[0];dy=y-control2[1];length=math.hypot(dx,dy);ux=dx/length;uy=dy/length
 tri=[(x,y),(x-13*ux-6*uy,y-13*uy+6*ux),(x-13*ux+6*uy,y-13*uy-6*ux)]
 draw.polygon(tri,fill=color);svg.append('<polygon points="'+' '.join(f'{a},{b}' for a,b in tri)+f'" fill="{color}"/>')
for i,d in enumerate(D):
 x=30+800*i;c=x+210
 role='Development cohort' if i==0 else 'External validation cohort'
 box(x,30,420,105,[d['db'],f'({role})'],blue,25,True,'white')
 line([(c,135),(c,200)],arrow=True)
 box(x,200,420,90,['Unique patients',f"n = {d['start']:,}"],light,25)
 line([(c,290),(c,560)],arrow=True)
 box(x+280,335,460,180,[f"Excluded (n = {d['exclude'][0]:,})",'Age <18 years','Study-defined sepsis criteria not met','Sepsis onset outside the ICU stay'],red,21)
 line([(c,420),(x+280,420)],arrow=True)
 box(x,560,420,100,['Adults with ICU-onset sepsis',f"n = {d['adult']:,}"],light,25)
 line([(c,660),(c,1010)],arrow=True)
 extra=28 if i==1 else 0;missing=d['final']-d['make']
 total=sum(d['exclude'][1:])+extra
 assert d['adult']-total==d['final']
 labels=[f'Excluded (n = {total:,})',
 f"ICU stay after sepsis <24 h (n = {d['exclude'][1]:,})",
 f"KRT within 24 h of sepsis (n = {d['exclude'][2]:,})",
 f"Insufficient trajectory data (n = {d['exclude'][3]:,})"]
 if extra:labels+=['Invalid time records (n = 28)']
 box(x+280,715,460,240,labels,red,21)
 line([(c,835),(x+280,835)],arrow=True)
 box(x,1010,420,100,['Final trajectory cohort',f"n = {d['final']:,}"],green,26,True)
 line([(c,1110),(c,1180)],arrow=True)
 labels=['Five-variable trajectory','model development'] if i==0 else ['Parameter transport','validation']
 box(x,1180,420,100,labels,purple,25,True)
curve((450,1230),(570,1205),(720,1205),(830,1230))
curve((450,1230),(820,1360),(1390,1360),(1630,1230))
text(640,1200,['Parameter transport'],21,color='#8561A2')
text(1080,1313,['Parameter transport'],21,color='#8561A2')
svg.append('</svg>')
(P/'Figure1.svg').write_text('\n'.join(svg),encoding='utf-8')
im.save(P/'Figure1.png')
print('Shen-style SVG and PNG generated.')
