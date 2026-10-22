from pathlib import Path
from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
P=Path('Z:/Research/MSI-1/results/候选补充材料_20261018')
for name, indexes in [('TableS3_MAKE30_components',range(5,9)),('TableS4_Subgroups_eICU',range(3,5)),('TableS5_Subgroups_SICdb',range(3,5))]:
 p=P/'表'/(name+'.docx');d=Document(p);t=d.tables[0]
 print(name,[[c.text for c in r.cells[:2]] for r in t.rows])
 for i in indexes:
  for tc in t.rows[i]._tr.findall(qn('w:tc')):
   pr=tc.get_or_add_tcPr()
   for old in pr.findall(qn('w:shd')):pr.remove(old)
   sh=OxmlElement('w:shd');sh.set(qn('w:fill'),'F5F5F5');sh.set(qn('w:val'),'clear');pr.append(sh)
 d.save(p)
