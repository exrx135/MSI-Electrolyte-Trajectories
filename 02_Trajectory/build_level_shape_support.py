from pathlib import Path
import pandas as pd,numpy as np,itertools,hashlib,json
from PIL import Image,ImageDraw,ImageFont
import html
from docx import Document
from docx.shared import Mm,Pt,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
ROOT=Path('Z:/Research/MSI-1');P=ROOT/'analysis_updates/G3_level_shape_support_20261019';P.mkdir(parents=True,exist_ok=True)
src = Path.home() / 'Desktop' / 'Research' / 'MSI' / '提取轨迹文件' / 'GBMTM_stage3_output' / 'V297_joint_portrait_G3G6' / 'V297_plot_data_G3G6.csv'; x=pd.read_csv(src); rows=[]
for g in [3,4,5,6]:
 z=x[x.requested_g==g]
 for a,b in itertools.combinations(range(1,g+1),2):
  j=z[z['class']==a].merge(z[z['class']==b],on=['indicator','time_hour_center'],suffixes=('_a','_b'))
  assert len(j)==28
  vals=[np.sqrt(((j[c+'_a']-j[c+'_b'])**2).groupby(j.indicator).mean().mean()) for c in ['fitted_z','clinical_centered_z']]
  rows.append([g,a,b,*vals,vals[1]/vals[0]])
t=pd.DataFrame(rows,columns=['classes','class_a','class_b','absolute_RMS','centered_RMS','centered_to_absolute_ratio']);assert len(t)==34
t.to_csv(P/'All_pairwise_distances.csv',index=False)
W,H=2400,2100
im=Image.new('RGB',(W,H),'white');dr=ImageDraw.Draw(im)
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><rect width="100%" height="100%" fill="white"/>']
def txt(x,y,ss,size=24,bold=False,anchor='middle'):
 f=ImageFont.truetype('C:/Windows/Fonts/arialbd.ttf' if bold else 'C:/Windows/Fonts/arial.ttf',size)
 dr.text((x,y),str(ss),font=f,fill='#20262B',anchor={'middle':'ms','start':'ls','end':'rs'}[anchor]);svg.append(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="Arial" font-size="{size}" font-weight="{700 if bold else 400}">{html.escape(str(ss))}</text>')
def line(x1,y1,x2,y2,color='#CDD2D6',width=3):
 dr.line((x1,y1,x2,y2),fill=color,width=width);svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"/>')
def mark(x,y,color,diamond=False):
 r=8
 if diamond:
  pts=[(x,y-r),(x+r,y),(x,y+r),(x-r,y)];dr.polygon(pts,fill=color);svg.append('<polygon points="'+' '.join(f'{a},{b}' for a,b in pts)+f'" fill="{color}"/>')
 else:dr.ellipse((x-r,y-r,x+r,y+r),fill=color);svg.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{color}"/>')
mark(530,60,'#DD1215');txt(550,70,'Absolute trajectories',29,anchor='start');mark(1200,60,'#0D4483',True);txt(1220,70,'Baseline-centered trajectories',29,anchor='start')
for ix,g in enumerate([3,4,5,6]):
 x0=200+(ix%2)*1180;y0=180+(ix//2)*970;pw=900;ph=790
 q=t[t.classes==g].sort_values(['class_a','class_b']);xx=lambda v:x0+v/.4*pw
 txt(x0+pw/2,y0-35,f'{g}-class solution',32,True)
 for tick in [0,.1,.2,.3,.4]:
  x=xx(tick);line(x,y0,x,y0+ph,'#EDF0F2',2);txt(x,y0+ph+38,f'{tick:.1f}',25)
 line(x0,y0+ph,x0+pw,y0+ph,'#333333',2)
 for j,r in enumerate(q.itertuples()):
  y=y0+(j+.5)/len(q)*ph
  txt(x0-25,y+8,f'C{r.class_a} vs C{r.class_b}',25,anchor='end')
  line(xx(r.centered_RMS),y,xx(r.absolute_RMS),y,width=4);mark(xx(r.absolute_RMS),y,'#DD1215');mark(xx(r.centered_RMS),y,'#0D4483',True)
 txt(x0+pw/2,y0+ph+85,'Pairwise RMS difference (standardized units)',26)
im.save(P/'Level and shape differences across candidate models.png',dpi=(300,300));(P/'Level and shape differences across candidate models.svg').write_text('\n'.join(svg)+ '</svg>',encoding='utf-8')
summary=[]
for g in [3,4,5,6]:
 q=t[t.classes==g];c=q.loc[q.centered_RMS.idxmin()]
 summary.append([g,len(q),q.absolute_RMS.median(),q.centered_RMS.median(),q.centered_to_absolute_ratio.median(),f'C{int(c.class_a)} vs C{int(c.class_b)}',c.centered_RMS])
s=pd.DataFrame(summary,columns=['classes','pairs','median_absolute_RMS','median_centered_RMS','median_ratio','closest_shape_pair','minimum_centered_RMS']);s.to_csv(P/'Model_level_summary.csv',index=False)
# Reuse the already approved three-line table layout.
helper=(ROOT/'历史归档/results_整理前_20261019/results/候选补充材料_20261018/代码/build_candidate_tables.py').read_text(encoding='utf-8').split('# S1 definitions')[0];ns={};exec(helper,ns)
d=ns['new']();p=d.add_paragraph('Level and shape separation across candidate trajectory models','Title');p.runs[0].bold=True
p=d.add_paragraph('Supporting material for model selection');p.paragraph_format.space_after=Pt(10)
values=[[str(g),str(n),f'{a:.4f}',f'{b:.4f}',f'{r:.3f}',pair,f'{mi:.4f}'] for g,n,a,b,r,pair,mi in summary]
ns['tbl'](d,['Classes','Class pairs','Median absolute\nRMS difference','Median centered\nRMS difference','Median centered /\nabsolute ratio','Closest pair\nafter centering','Minimum centered\nRMS difference'],values,[65,75,135,135,135,135,135],font=10,height=38)
ns['note'](d,'Distances were calculated from saved candidate-model fitted values at the original observation-window anchors. Each variable was standardized using its saved scale; its first fitted value was subtracted for baseline centering. Squared differences were averaged over time within each variable, averaged equally across five variables, and square-rooted. All within-model class pairs are included. Ratios are descriptive distance ratios, not explained-variance fractions. Class labels are specific to each candidate model and must not be matched by number across models. No hypothesis tests or new trajectory fits were performed.')
d.save(P/'Level and shape separation across candidate trajectory models.docx')
legend='Each panel includes all within-model class pairs from the saved three- to six-class candidate solutions. Red circles show RMS differences between absolute standardized trajectories; blue diamonds show differences after subtracting the first fitted value separately for each variable and class. Gray lines connect the two distances for the same pair. Distances give equal weight to each of the five variables and use only the original window anchors (12 for MSI and four for each electrolyte). Smaller centered distances indicate more similar changes from the first window despite differences in absolute levels. Class numbers are model-specific. These are descriptive comparisons of the initial candidate fits, not tests of an optimal class number; no trajectory model was refitted.'
(P/'Figure legend.txt').write_text(legend,encoding='utf-8')
(P/'Interpretation.md').write_text('这是G3选择的支持材料，不自动加入当前正式补充编号。全部34个类别对均展示，避免仅挑选支持性例子。G3至G6的中心化/绝对距离比中位数分别为0.525、0.344、0.282、0.276，提示复杂模型中更多类间分离来自起始水平。G3内部也存在水平差异，故中心化后缩小不是G3独有优越性的检验。总体支持“新增细分主要体现水平差异”，不能单独证明G3为唯一正确类别数，不能把不同G的同编号类别当作同一类。最终选择结合V298当时的最小充分原则、临床解释及分类质量。\n',encoding='utf-8')
(P/'source_manifest.json').write_text(json.dumps({'source':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'pairs':34,'no_refitting':True},indent=2),encoding='utf-8')
print(s.to_string(index=False))
