"""Copy code into a NEW private execution workspace; never copy any patient data.

Usage: python tools/prepare_private_workspace.py /absolute/new/private/workspace
This only prepares source files. It never runs SQL, fits a model or opens a database.
"""
from pathlib import Path
import argparse
import hashlib
import json

REPO = Path(__file__).resolve().parents[1]

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('workspace', type=Path)
    args = ap.parse_args()
    target = args.workspace.resolve()
    if target == REPO or target.is_relative_to(REPO):
        ap.error('The private workspace must be outside this public-code repository.')
    if target.exists():
        ap.error('Choose a new directory. Existing directories are never overwritten.')
    entries = json.loads((REPO/'manifest.json').read_text(encoding='utf-8'))
    for e in entries:
        src = REPO / e['file']
        if hashlib.sha256(src.read_bytes()).hexdigest() != e['sha256']:
            raise ValueError('Source file has changed since packaging: '+e['file'])
        if not (target/e['execution_relative_path']).resolve().is_relative_to(target):
            raise ValueError('Invalid execution path in manifest.')
    target.mkdir(parents=True)
    project = target.as_posix()
    source_msi = Path.home() / 'Desktop' / 'Research' / 'MSI'
    replacements = [
        ('Z:/Research/MSI-1/历史归档/results_整理前_20261019/results',project+'/results'),
        ('Z:/Research/MSI-1',project),
        ('Z:\\Research\\MSI-1',str(target)),
        (source_msi.as_posix(),project+'/MSI'),
        (str(source_msi),str(target/'MSI')),
        ('历史归档/results_整理前_20261019/results','results'),
        (r'\u5386\u53f2\u5f52\u6863/results_\u6574\u7406\u524d_20261019/results','results'),
    ]
    copied=[]
    for e in entries:
        s=(REPO/e['file']).read_text(encoding='utf-8-sig')
        for old,new in replacements:s=s.replace(old,new)
        dst=target/e['execution_relative_path'];dst.parent.mkdir(parents=True,exist_ok=True)
        if dst.exists():raise FileExistsError(dst)
        dst.write_text(s,encoding='utf-8')
        copied.append(e['execution_relative_path'])
    # Data-free configuration and final subgroup note required by the renderers.
    for src,rel in [
        ('03_Analysis/model_specification.json','results/关联分析_三模型_MICE_RR/model_specification.json'),
        ('05_Tables/table_legend.txt','results/补充表_亚组与交互作用/table_legend.txt'),
        ('01_SQL/00_create_private_input_tables.sql','SQL/00_create_private_input_tables.sql'),
    ]:
        dst=target/rel;dst.parent.mkdir(parents=True,exist_ok=True)
        if dst.exists():raise FileExistsError(dst)
        dst.write_bytes((REPO/src).read_bytes())
    (target/'private_inputs').mkdir()
    (target/'private_workspace_manifest.json').write_text(json.dumps(copied,ensure_ascii=False,indent=2),encoding='utf-8')
    print(f'Prepared {len(entries)} source files in {target}. No analysis was run.')
    print('Supply authorised private inputs, check module-specific paths and read docs/WORKFLOW.md before execution.')

if __name__=='__main__':main()
