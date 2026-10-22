from pathlib import Path
import pandas as pd
import numpy as np
import json

P=Path(__file__).resolve().parent.parent
dbs=['MIMIC','eICU','SICdb']
for db in dbs:
 log=(P/'diagnostics'/f'{db}_final_run_log.txt').read_text(encoding='utf-8',errors='replace')
 assert 'ALL ANALYSES COMPLETE' in log,(db,'final run not complete')
 assert 'Execution halted' not in log,(db,'final R error')
targets=['pooled_RR_results.csv','diagnostics/model_fit_audit.csv',
         'diagnostics/per_imputation_coefficients.csv','diagnostics/imputed_distributions.csv',
         'diagnostics/crude_direct_ratio_check.csv']
for target in targets:
 q=P/target
 parts=[q.parent/(db+'_'+q.name) for db in dbs]
 if all(f.exists() for f in parts):
  pd.concat([pd.read_csv(f) for f in parts],ignore_index=True).to_csv(q,index=False)
results=pd.read_csv(P/'pooled_RR_results.csv')
draws=pd.read_csv(P/'diagnostics/per_imputation_coefficients.csv')
fits=pd.read_csv(P/'diagnostics/model_fit_audit.csv')
assert len(results)==18 and len(fits)==3*(1+30*2)
assert not results.duplicated(['dataset','model','comparison']).any()
assert fits.converged.all() and fits.warning.fillna('').eq('').all()
expected={'MIMIC':(9758,2442),'eICU':(5709,2060),'SICdb':(3635,949)}
for row in results.itertuples():
 assert (row.n,row.events)==expected[row.dataset]
 z=draws[(draws.dataset==row.dataset)&(draws.model==row.model)&(draws.comparison==row.comparison)]
 assert len(z)==(1 if row.model=='Model 1' else 30)
 q=z.Q.mean();u=z.U.mean();b=0 if len(z)==1 else z.Q.var(ddof=1)
 t=u+(1+1/len(z))*b
 assert abs(np.log(row.RR)-q)<1e-10
 assert abs(row.SE_log_RR**2-t)<1e-10
 assert 0<row.CI95_low<row.RR<row.CI95_high and 0<=row.p<=1
direct=pd.read_csv(P/'diagnostics/crude_direct_ratio_check.csv')
assert direct.difference.abs().max()<1e-7
summary=fits.groupby(['dataset','model']).agg(n=('n','first'),max_fitted_risk=('max_fitted_risk','max'),min_risks_gt1=('fitted_risk_gt1','min'),max_risks_gt1=('fitted_risk_gt1','max'),max_vif=('max_vif','max'),max_hat=('max_hat','max')).reset_index()
summary['max_pct_risks_gt1']=100*summary.max_risks_gt1/summary.n
summary.to_csv(P/'diagnostics/model_fit_summary.csv',index=False)
results['RR_95CI']=results.apply(lambda r:f'{r.RR:.2f} ({r.CI95_low:.2f}–{r.CI95_high:.2f})',axis=1)
display=results.pivot(index=['comparison','dataset'],columns='model',values='RR_95CI').reindex(pd.MultiIndex.from_product([['C2 vs C1','C3 vs C1'],dbs],names=['comparison','dataset']))
display.to_csv(P/'RR_三模型对照表.csv',encoding='utf-8-sig')
convpath=P/'diagnostics/imputation_convergence_summary.csv'
conv=pd.read_csv(convpath) if convpath.exists() else pd.DataFrame()
if not conv.empty:
 assert set(conv.dataset)==set(dbs)
 assert conv.iterations.min()>=30
mc=(results.MCSE_log_RR/results.SE_log_RR).max()
text=['# 本次运行核查','', '各库纳入人数与冻结MAKE30分析完全一致；未改变分型或结局。',
      'Model 1直接发生率比校验通过；18个结果的Rubin合并均值、总方差经Python独立核对通过。',
      f'共{len(fits)}个回归拟合全部收敛，无秩缺失；30份插补之间的类别系数Monte Carlo标准误/总标准误最大为{mc:.4f}。','',
      '## 拟合检查（不隐去Poisson模型的限制）',summary.to_string(index=False),'',
      '大于1的拟合值是modified Poisson的已知局限，提示相应个体的概率均值模型不能按真实预测概率解释；HC1不消除此问题。这些结果是关联RR，不是风险预测器。未因该诊断擅自换模型、删人或根据显著性筛选变量。',
      '', '## 插补链诊断',conv.to_string(index=False) if not conv.empty else '尚未生成',
      '', '经典split-Rhat用于链均值/标准差的描述性检查；检查通过不证明MAR成立。',
      '', '## 三模型结果',display.to_string()]
(P/'运行核查与限制.md').write_text('\n'.join(text),encoding='utf-8')
print(display.to_string());print(summary.to_string(index=False));print('Maximum coefficient MCSE / SE:',mc)
