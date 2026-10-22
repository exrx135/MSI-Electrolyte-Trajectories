from pathlib import Path
import pandas as pd,re
root=Path(__file__).resolve().parent.parent;src=root/'data_and_sql'
out=src/'SICdb_呼吸频率与体温补提取';out.mkdir(exist_ok=True)
old=(src/'03_sicdb_独立查询.sql').read_text(encoding='utf-8')
helper=old[old.index('CREATE OR REPLACE FUNCTION pg_temp.t1_sic_native'):old.index('DO $load$ BEGIN PERFORM pg_temp.t1_sic_native')]
helper=helper.replace('case_col text; time_col text;', 'case_col text; time_col text; data_col text;')
helper=helper.replace("IF f->>'name'='offset_sec' THEN time_col:=col; END IF;", "IF f->>'name'='offset_sec' THEN time_col:=col; END IF;\n IF f->>'name'='data_id' THEN data_col:=col; END IF;")
helper=helper.replace('IF use_window THEN',"IF native_table='data_float_h' THEN sql_text:=sql_text||format(' JOIN sic_signal_candidates c ON d.%I=c.ref_id',data_col); END IF;\n IF use_window THEN")
loads={}
for name in ['sic_src_refs','sic_src_data_float_h']:
 loads[name]=next(line for line in old.splitlines() if line.startswith('DO $load$') and "('"+name+"'" in line)
setup='''-- Standalone PostgreSQL query. Native SICdb tables only; no earlier study tables required.
-- Only frozen 3831 cases; [sepsis onset, onset + 86400 seconds).
-- Read-only native sources; all generated objects are temporary.
CREATE TEMP TABLE IF NOT EXISTS sic_vital_session(dummy integer);
'''+helper+loads['sic_src_refs']+'''
DROP TABLE IF EXISTS pg_temp.sic_signal_candidates;
CREATE TEMP TABLE sic_signal_candidates AS
SELECT *, CASE WHEN concat_ws(' ',ref_value,ref_name) ~* '(temp|temperatur)' THEN 'temperature' ELSE 'respiratory' END AS candidate_family
FROM sic_src_refs
WHERE concat_ws(' ',ref_value,ref_name) ~* '(resp|breath|atem|frequen|temp|ventilat)' OR ref_id=709;
-- Broad matching intentionally retains measured, spontaneous, total and SET channels separately.
'''
a=pd.read_csv(src/'03_sicdb_独立查询.csv',usecols=['cohort_unit_id','onset_offset','assigned_class'])
assert len(a)==3831 and a.cohort_unit_id.is_unique and a.onset_offset.notna().all()
seed='DROP TABLE IF EXISTS pg_temp.t1_seed;\nCREATE TEMP TABLE t1_seed(unit_id bigint PRIMARY KEY,onset_offset numeric,assigned_class integer);\nINSERT INTO t1_seed VALUES\n'
seed+=',\n'.join(f'({int(r.cohort_unit_id)},{r.onset_offset},{int(r.assigned_class)})' for r in a.itertuples())+';\n'
tail='''
DROP TABLE IF EXISTS pg_temp.sic_vital_agg;
CREATE TEMP TABLE sic_vital_agg AS
SELECT case_id,data_id,avg(val) FILTER (WHERE val::text NOT IN ('NaN','Infinity','-Infinity')) AS mean_value,
count(val) FILTER (WHERE val::text NOT IN ('NaN','Infinity','-Infinity')) AS n_values
FROM sic_src_data_float_h GROUP BY case_id,data_id;
CREATE INDEX ON sic_vital_agg(case_id);
DO $wide$
DECLARE r record; cols text:='';
BEGIN
 FOR r IN SELECT c.* FROM sic_signal_candidates c
 WHERE EXISTS(SELECT 1 FROM sic_vital_agg a WHERE a.data_id=c.ref_id)
 ORDER BY c.candidate_family,c.ref_id LOOP
 cols:=cols||format(',max(a.mean_value) FILTER (WHERE a.data_id=%s) AS %I,max(a.n_values) FILTER (WHERE a.data_id=%s) AS %I',
 r.ref_id,r.candidate_family||'_'||r.ref_id||'_mean',r.ref_id,r.candidate_family||'_'||r.ref_id||'_n');
 END LOOP;
 EXECUTE 'DROP TABLE IF EXISTS pg_temp.sic_vital_export';
 EXECUTE 'CREATE TEMP TABLE sic_vital_export AS SELECT s.unit_id AS case_id,s.onset_offset AS sepsis_onset_offset_sec,s.assigned_class'||cols||
 ' FROM t1_seed s LEFT JOIN sic_vital_agg a ON a.case_id=s.unit_id GROUP BY s.unit_id,s.onset_offset,s.assigned_class';
END $wide$;
-- EXPORT THIS FINAL RESULT: one patient per row; each candidate has mean and count columns.
SELECT * FROM sic_vital_export ORDER BY case_id;
'''
(out/'01_SICdb_呼吸频率及体温_独立查询.sql').write_text(setup+seed+loads['sic_src_data_float_h']+tail,encoding='utf-8')
(out/'02_SICdb_候选字段字典_独立查询.sql').write_text(setup+'''
SELECT ref_id AS source_id,candidate_family,ref_value AS source_label,ref_name AS source_group,ref_unit AS native_unit,
candidate_family||'_'||ref_id||'_mean' AS mean_column,
candidate_family||'_'||ref_id||'_n' AS count_column
FROM sic_signal_candidates ORDER BY candidate_family,ref_id;
''',encoding='utf-8')
(out/'运行说明.txt').write_text('''请在SICdb数据库分别完整运行01和02，导出各自最后的查询结果为同名CSV，放回本文件夹。
两个SQL都独立运行，不依赖之前查询生成的表或同一个会话。
01内置3831名最终纳入患者的ID和已冻结sepsis onset，不重选队列。结果每人一行，每个候选指标分别导出均值和有效记录数，没有JSON字段。
02只查询字典，很快；其候选项可能比01实际出现的列多，这是正常情况。
呼吸相关检索覆盖resp、breath、atem、frequen、ventilat；体温覆盖temp、temperatur及已有709项目。
为避免误将呼吸机设置频率、分钟通气量等当成患者呼吸频率，候选项分列导出，拿到结果后结合字典选择。体温不同部位也保留分列。
本次使用原生data_float_h小时数据，均值是24小时窗口内存储值的算术均值；不解压分钟数据、不自行替换为呼吸机设置值、不做额外极端值过滤。此结果不代表已经穷尽所有原始数据表。
仅创建临时表，不修改原始数据库。当前环境无数据库连接，SQL尚未在你的数据库执行验证。
''',encoding='utf-8')
print(out)
