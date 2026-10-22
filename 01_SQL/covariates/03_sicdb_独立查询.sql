-- STANDALONE: run only this file in the corresponding native database.
-- Requires privately populated local_inputs.sicdb_baseline_seed; no patient records are distributed.
-- Intermediate objects are created within this run. Only the FINAL SELECT returns rows.
-- If the connection is still in a failed transaction, execute ROLLBACK first.
BEGIN;

DROP TABLE IF EXISTS pg_temp.t1_seed;
DROP TABLE IF EXISTS pg_temp.t1_payload;
DROP TABLE IF EXISTS pg_temp.t1_audit;
CREATE TEMP TABLE t1_seed (
 unit_id bigint PRIMARY KEY, hospital_id bigint, patient_id text,
 onset_timestamp timestamp, onset_offset numeric,
 assigned_class integer, frozen_baseline jsonb NOT NULL
);
CREATE TEMP TABLE t1_payload(unit_id bigint, source text, records jsonb,
 PRIMARY KEY(unit_id,source));
CREATE TEMP TABLE t1_audit(source text PRIMARY KEY, status text,
 patients_with_records bigint, record_count bigint, detail text);
CREATE OR REPLACE FUNCTION pg_temp.t1_collect(
 source_name text, dependencies text[], query_text text, required boolean DEFAULT false
) RETURNS void LANGUAGE plpgsql AS $fn$
DECLARE dep text; n bigint; missing text := '';
BEGIN
 FOREACH dep IN ARRAY dependencies LOOP
  IF to_regclass(dep) IS NULL THEN missing := missing || dep || '; '; END IF;
 END LOOP;
 IF missing <> '' THEN
  IF required THEN RAISE EXCEPTION 'Required source missing: %',missing; END IF;
  INSERT INTO t1_audit VALUES(source_name,'SOURCE_MISSING',0,0,missing);
  RAISE WARNING 'Optional source % missing: %. This is missing data, NOT a negative finding.', source_name,missing;
  RETURN;
 END IF;
 EXECUTE format('INSERT INTO t1_payload SELECT unit_id, %L, jsonb_agg(payload) FROM (%s) q GROUP BY unit_id',source_name,query_text);
 GET DIAGNOSTICS n = ROW_COUNT;
 IF source_name='identity' AND n<>(SELECT count(*) FROM t1_seed) THEN
  RAISE EXCEPTION 'Frozen cohort IDs do not all match database identity: % matched',n;
 END IF;
 INSERT INTO t1_audit
 SELECT source_name, CASE WHEN n=0 THEN 'AVAILABLE_BUT_NO_MATCHES' ELSE 'AVAILABLE' END,
 n,COALESCE(sum(jsonb_array_length(records)),0),array_to_string(dependencies,', ')
 FROM t1_payload WHERE source=source_name;
 RAISE NOTICE 'Collected %: % patients',source_name,n;
END $fn$;

-- Public-code edition: patient rows are loaded from a private local table.
INSERT INTO t1_seed SELECT * FROM local_inputs.sicdb_baseline_seed;
DO $seed_check$ BEGIN IF NOT EXISTS (SELECT 1 FROM t1_seed) THEN
RAISE EXCEPTION 'Populate local_inputs.sicdb_baseline_seed with authorised local data first';
END IF; END $seed_check$;
ANALYZE t1_seed;

-- Build current-session normalized sources DIRECTLY from native SICdb tables.
-- This does not read any study_make30_final or previous Table 1 output table.
CREATE OR REPLACE FUNCTION pg_temp.t1_sic_native(out_table text, native_table text,
 fields jsonb, required boolean DEFAULT true, use_window boolean DEFAULT false)
RETURNS void LANGUAGE plpgsql AS $native$
DECLARE ns text; tbl text; candidates integer; f jsonb; col text; expr text := '';
 case_col text; time_col text; sql_text text;
BEGIN
 SELECT count(DISTINCT n.nspname),min(n.nspname) INTO candidates,ns
 FROM pg_namespace n WHERE n.nspname NOT LIKE 'pg_%' AND n.nspname <> 'information_schema'
 AND EXISTS(SELECT 1 FROM pg_class c WHERE c.relnamespace=n.oid AND lower(c.relname)='cases' AND c.relkind IN ('r','p','v','m','f'))
 AND EXISTS(SELECT 1 FROM pg_class c WHERE c.relnamespace=n.oid AND lower(c.relname)='d_references' AND c.relkind IN ('r','p','v','m','f'));
 IF candidates<>1 THEN RAISE EXCEPTION 'Expected one native SICdb schema containing cases+d_references; found %',candidates; END IF;
 SELECT c.relname INTO tbl FROM pg_class c JOIN pg_namespace n ON n.oid=c.relnamespace
 WHERE n.nspname=ns AND lower(c.relname)=native_table AND c.relkind IN ('r','p','v','m','f');
 IF tbl IS NULL THEN
  IF required THEN RAISE EXCEPTION 'Native SICdb source %.% missing',ns,native_table; END IF;
  RETURN;
 END IF;
 FOR f IN SELECT value FROM jsonb_array_elements(fields) LOOP
  SELECT a.attname INTO col FROM pg_attribute a
  WHERE a.attrelid=format('%I.%I',ns,tbl)::regclass AND a.attnum>0 AND NOT a.attisdropped
  AND regexp_replace(lower(a.attname),'[^a-z0-9]','','g') IN (SELECT jsonb_array_elements_text(f->'aliases'))
  ORDER BY (SELECT min(z.ordinality) FROM jsonb_array_elements_text(f->'aliases') WITH ORDINALITY z(value,ordinality)
    WHERE z.value=regexp_replace(lower(a.attname),'[^a-z0-9]','','g')) LIMIT 1;
  IF col IS NULL AND COALESCE((f->>'required')::boolean,false) THEN
   RAISE EXCEPTION 'Native %.% missing required field %',ns,tbl,f->>'name'; END IF;
  IF expr<>'' THEN expr:=expr||','; END IF;
  expr:=expr||CASE WHEN col IS NULL THEN format('NULL::%s AS %I',f->>'type',f->>'name')
    ELSE format('d.%I::%s AS %I',col,f->>'type',f->>'name') END;
  IF f->>'name'='case_id' THEN case_col:=col; END IF;
  IF f->>'name'='offset_sec' THEN time_col:=col; END IF;
 END LOOP;
 IF native_table='cases' THEN expr:=expr||',to_jsonb(d) AS raw_record'; END IF;
 EXECUTE format('DROP TABLE IF EXISTS pg_temp.%I',out_table);
 sql_text:=format('CREATE TEMP TABLE %I AS SELECT %s FROM %I.%I d',out_table,expr,ns,tbl);
 IF case_col IS NOT NULL THEN sql_text:=sql_text||format(' JOIN t1_seed s ON d.%I::bigint=s.unit_id',case_col); END IF;
 IF use_window THEN
  IF case_col IS NULL OR time_col IS NULL THEN RAISE EXCEPTION 'Time window columns not resolved'; END IF;
  sql_text:=sql_text||format(' WHERE d.%I::numeric>=s.onset_offset AND d.%I::numeric<s.onset_offset+86400',time_col,time_col);
 END IF;
 EXECUTE sql_text;
 IF case_col IS NOT NULL THEN EXECUTE format('CREATE INDEX ON %I(case_id)',out_table); END IF;
 RAISE NOTICE 'Native source: %.% -> current session %',ns,tbl,out_table;
END $native$;

DO $load$ BEGIN PERFORM pg_temp.t1_sic_native('sic_src_cases','cases','[{"name": "case_id", "aliases": ["caseid", "case", "admissionid", "stayid"], "type": "bigint", "required": true}, {"name": "patient_id", "aliases": ["patientid"], "type": "bigint", "required": true}, {"name": "time_of_stay_sec", "aliases": ["timeofstay"], "type": "bigint", "required": false}, {"name": "icu_offset_sec", "aliases": ["icuoffset"], "type": "bigint", "required": false}, {"name": "age_years", "aliases": ["ageonadmission", "age"], "type": "numeric", "required": false}, {"name": "offset_after_first_admission_sec", "aliases": ["offsetafterfirstadmission"], "type": "bigint", "required": false}, {"name": "offset_of_death_sec", "aliases": ["offsetofdeath"], "type": "bigint", "required": false}, {"name": "hospital_discharge_day", "aliases": ["hospitaldischargeday"], "type": "integer", "required": false}, {"name": "hospital_stay_days", "aliases": ["hospitalstaydays"], "type": "integer", "required": false}, {"name": "weight_kg", "aliases": ["weightonadmission"], "type": "numeric", "required": false}, {"name": "height_cm", "aliases": ["heightonadmission"], "type": "numeric", "required": false}, {"name": "saps3", "aliases": ["saps3"], "type": "numeric", "required": false}, {"name": "sex_ref_id", "aliases": ["sex"], "type": "integer", "required": false}, {"name": "hospital_discharge_type_ref", "aliases": ["hospitaldischargetype"], "type": "integer", "required": false}, {"name": "icd10_main", "aliases": ["icd10main"], "type": "text", "required": false}, {"name": "icd10_main_text", "aliases": ["icd10maintext"], "type": "text", "required": false}]'::jsonb,true,false); END $load$;

DO $load$ BEGIN PERFORM pg_temp.t1_sic_native('sic_src_refs','d_references','[{"name": "ref_id", "aliases": ["referenceglobalid"], "type": "integer", "required": true}, {"name": "ref_value", "aliases": ["referencevalue"], "type": "text", "required": true}, {"name": "ref_name", "aliases": ["referencename"], "type": "text", "required": false}, {"name": "ref_unit", "aliases": ["referenceunit"], "type": "text", "required": false}]'::jsonb,true,false); END $load$;

DO $load$ BEGIN PERFORM pg_temp.t1_sic_native('sic_src_laboratory','laboratory','[{"name": "case_id", "aliases": ["caseid", "case", "admissionid", "stayid"], "type": "bigint", "required": true}, {"name": "patient_id", "aliases": ["patientid"], "type": "bigint", "required": false}, {"name": "lab_id", "aliases": ["laboratoryid", "drugid", "labid", "referenceglobalid", "referenceid", "itemid"], "type": "integer", "required": true}, {"name": "offset_sec", "aliases": ["offset", "offsetvalue", "offsetsec", "timeoffset", "chartoffset", "labresultoffset"], "type": "bigint", "required": true}, {"name": "lab_value", "aliases": ["laboratoryvalue", "labvalue", "labresult", "value", "val"], "type": "double precision", "required": true}, {"name": "lab_type", "aliases": ["laboratorytype", "labtype"], "type": "integer", "required": false}]'::jsonb,true,true); END $load$;

DO $load$ BEGIN PERFORM pg_temp.t1_sic_native('sic_src_data_float_h','data_float_h','[{"name": "case_id", "aliases": ["caseid", "case", "admissionid", "stayid"], "type": "bigint", "required": true}, {"name": "data_id", "aliases": ["dataid", "fieldid", "referenceglobalid", "referenceid", "variableid", "itemid"], "type": "integer", "required": true}, {"name": "offset_sec", "aliases": ["offset", "offsetvalue", "offsetsec", "timeoffset", "chartoffset", "labresultoffset"], "type": "bigint", "required": true}, {"name": "val", "aliases": ["val", "value", "valuefloat", "floatvalue", "mean"], "type": "double precision", "required": true}]'::jsonb,true,true); END $load$;

DO $load$ BEGIN PERFORM pg_temp.t1_sic_native('sic_src_data_ref','data_ref','[{"name": "case_id", "aliases": ["caseid", "case", "admissionid", "stayid"], "type": "bigint", "required": true}, {"name": "field_id", "aliases": ["fieldid"], "type": "integer", "required": true}, {"name": "ref_id", "aliases": ["refid"], "type": "integer", "required": true}]'::jsonb,true,false); END $load$;
DROP TABLE IF EXISTS pg_temp.sic_case_extra;
CREATE TEMP TABLE sic_case_extra AS SELECT case_id,sex_ref_id,weight_kg,height_cm,saps3,hospital_discharge_type_ref,icd10_main,icd10_main_text FROM sic_src_cases;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('identity', ARRAY['pg_temp.sic_src_cases'], $query$
SELECT s.unit_id,to_jsonb(d) payload FROM t1_seed s JOIN pg_temp.sic_src_cases d ON d.case_id=s.unit_id AND d.patient_id::text=s.patient_id
$query$, true); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('case_extra', ARRAY['pg_temp.sic_case_extra'], $query$
SELECT s.unit_id,to_jsonb(d) payload FROM t1_seed s JOIN pg_temp.sic_case_extra d ON d.case_id=s.unit_id
$query$, false); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('labs_24h', ARRAY['pg_temp.sic_src_laboratory','pg_temp.sic_src_refs'], $query$
SELECT s.unit_id,to_jsonb(d)||jsonb_build_object('dictionary',to_jsonb(r),'specimen_dictionary',to_jsonb(t)) payload
FROM t1_seed s JOIN pg_temp.sic_src_laboratory d ON d.case_id=s.unit_id
LEFT JOIN pg_temp.sic_src_refs r ON r.ref_id=d.lab_id
LEFT JOIN pg_temp.sic_src_refs t ON t.ref_id=d.lab_type
WHERE d.offset_sec>=s.onset_offset AND d.offset_sec<s.onset_offset+86400
$query$, true); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('hourly_signals_24h', ARRAY['pg_temp.sic_src_data_float_h','pg_temp.sic_src_refs'], $query$
SELECT s.unit_id,to_jsonb(d)||jsonb_build_object('dictionary',to_jsonb(r)) payload
FROM t1_seed s JOIN pg_temp.sic_src_data_float_h d ON d.case_id=s.unit_id LEFT JOIN pg_temp.sic_src_refs r ON r.ref_id=d.data_id
WHERE d.offset_sec>=s.onset_offset AND d.offset_sec<s.onset_offset+86400
$query$, true); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('structured_history', ARRAY['pg_temp.sic_src_data_ref','pg_temp.sic_src_refs'], $query$
SELECT s.unit_id,to_jsonb(d)||jsonb_build_object('field_dictionary',to_jsonb(f),'value_dictionary',to_jsonb(r)) payload
FROM t1_seed s JOIN pg_temp.sic_src_data_ref d ON d.case_id=s.unit_id
LEFT JOIN pg_temp.sic_src_refs f ON f.ref_id=d.field_id LEFT JOIN pg_temp.sic_src_refs r ON r.ref_id=d.ref_id
$query$, true); END $collect$;
DO $collect$ BEGIN PERFORM pg_temp.t1_collect('raw_cases:native',ARRAY['pg_temp.sic_src_cases'],
        'SELECT s.case_id AS unit_id,s.raw_record AS payload FROM pg_temp.sic_src_cases s',true); END $collect$;
DROP TABLE IF EXISTS pg_temp.t1_current_data;
CREATE TEMP TABLE t1_current_data AS SELECT s.unit_id,s.hospital_id,s.patient_id,
s.onset_timestamp AS sepsis_onset_timestamp,s.onset_offset AS sepsis_onset_offset,
'seconds'::text AS onset_offset_unit,
s.assigned_class,s.frozen_baseline,COALESCE(p.data,'{}'::jsonb) AS collected_data
FROM t1_seed s LEFT JOIN (SELECT unit_id,jsonb_object_agg(source,records) data FROM t1_payload GROUP BY unit_id) p USING(unit_id);
CREATE UNIQUE INDEX ON t1_current_data(unit_id);
CREATE SCHEMA IF NOT EXISTS table1_direct_20261016;
DROP TABLE IF EXISTS table1_direct_20261016.sicdb_source_audit;
CREATE TABLE table1_direct_20261016.sicdb_source_audit AS SELECT * FROM t1_audit;

CREATE TEMP TABLE IF NOT EXISTS t1_wide_session_init(dummy integer);
CREATE OR REPLACE FUNCTION pg_temp.t1_num(v text) RETURNS numeric
LANGUAGE sql IMMUTABLE AS $$
 SELECT CASE WHEN btrim(v) ~ '^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$'
 THEN btrim(v)::numeric ELSE NULL END
$$;
CREATE OR REPLACE FUNCTION pg_temp.t1_norm(v text) RETURNS text
LANGUAGE sql IMMUTABLE AS $$ SELECT regexp_replace(lower(COALESCE(v,'')),'[^a-z0-9]','','g') $$;
CREATE OR REPLACE FUNCTION pg_temp.t1_lab(label text, code text, db text) RETURNS text
LANGUAGE plpgsql IMMUTABLE AS $$
DECLARE n text := pg_temp.t1_norm(label);
BEGIN
 IF db='mimic' THEN
  RETURN CASE code
   WHEN '51301' THEN 'wbc' WHEN '51222' THEN 'hemoglobin' WHEN '51265' THEN 'platelets'
   WHEN '50912' THEN 'creatinine' WHEN '51006' THEN 'bun' WHEN '50861' THEN 'alt'
   WHEN '50878' THEN 'ast' WHEN '50862' THEN 'albumin' WHEN '50885' THEN 'bilirubin'
   WHEN '50813' THEN 'lactate' WHEN '50931' THEN 'glucose' WHEN '50983' THEN 'sodium'
   WHEN '50971' THEN 'potassium' WHEN '50902' THEN 'chloride' WHEN '50882' THEN 'bicarbonate'
   WHEN '50868' THEN 'aniongap' WHEN '51237' THEN 'inr' WHEN '51274' THEN 'pt'
   WHEN '51275' THEN 'ptt' WHEN '50820' THEN 'ph' WHEN '50818' THEN 'pco2'
   WHEN '50821' THEN 'po2' WHEN '50816' THEN 'fio2' WHEN '50893' THEN 'calcium'
   WHEN '50960' THEN 'magnesium' WHEN '50970' THEN 'phosphate' END;
 END IF;
 IF db='sicdb' THEN
  IF code='367' THEN RETURN 'creatinine'; END IF;
  IF code IN ('469','686') THEN RETURN 'sodium'; END IF;
  IF code IN ('463','685') THEN RETURN 'potassium'; END IF;
  IF code IN ('450','683') THEN RETURN 'chloride'; END IF;
  IF code IN ('456','666') THEN RETURN 'bicarbonate'; END IF;
 END IF;
 RETURN CASE
  WHEN n IN ('wbc','wbcx1000','whitebloodcells','whitebloodcellcount','leukocytes','leukozyten') THEN 'wbc'
  WHEN n IN ('hgb','hb','hemoglobin','haemoglobin','hmoglobin') THEN 'hemoglobin'
  WHEN n IN ('platelets','plateletcount','plt','thrombocytes','thrombozyten') THEN 'platelets'
  WHEN n IN ('creatinine','creatinin','kreatinin','crea') THEN 'creatinine'
  WHEN n IN ('bun','ureanitrogen','bloodureanitrogen') THEN 'bun'
  WHEN n IN ('urea','harnstoff') THEN 'urea'
  WHEN n IN ('alt','sgpt','altsgpt','alanineaminotransferase','gpt') THEN 'alt'
  WHEN n IN ('ast','sgot','astsgot','aspartateaminotransferase','got') THEN 'ast'
  WHEN n IN ('albumin','albuminlevel') THEN 'albumin'
  WHEN n IN ('totalbilirubin','bilirubintotal','bilirubin','gesamtbilirubin') THEN 'bilirubin'
  WHEN n IN ('lactate','lacticacid','laktat') THEN 'lactate'
  WHEN n IN ('glucose','glukose','bloodglucose','bedsideglucose') THEN 'glucose'
  WHEN n IN ('sodium','na','natrium') THEN 'sodium'
  WHEN n IN ('potassium','k','kalium') THEN 'potassium'
  WHEN n IN ('chloride','cl','chlorid') THEN 'chloride'
  WHEN n IN ('bicarbonate','hco3','hco3actual','bicarbonateactual','totalco2','co2','hco3act') THEN 'bicarbonate'
  WHEN n IN ('aniongap','anionenlcke') THEN 'aniongap'
  WHEN n IN ('inr','ptinr') THEN 'inr'
  WHEN n IN ('pt','prothrombintime') THEN 'pt'
  WHEN n IN ('ptt','aptt','partialthromboplastintime') THEN 'ptt'
  WHEN n IN ('ph','pharterial') THEN 'ph'
  WHEN n IN ('paco2','pco2','paco2arterial') THEN 'pco2'
  WHEN n IN ('pao2','po2','pao2arterial') THEN 'po2'
  WHEN n IN ('fio2','fractionofinspiredoxygen') THEN 'fio2'
  WHEN n IN ('calcium','totalcalcium','calciumtotal') THEN 'calcium'
  WHEN n IN ('magnesium','mg') THEN 'magnesium'
  WHEN n IN ('phosphate','phosphorus','phosphat') THEN 'phosphate' END;
END $$;
DROP TABLE IF EXISTS pg_temp.t1_features;
CREATE TEMP TABLE t1_features(unit_id bigint, domain text, variable text,
 source_key text, source_label text, native_unit text, specimen text,
 value numeric, weight numeric, method text);
DO $$ BEGIN
 IF to_regclass('pg_temp.t1_current_data') IS NULL THEN RAISE EXCEPTION 'Collected table pg_temp.t1_current_data missing. Stop and report; do not rerun old extraction.'; END IF;
 IF (SELECT count(*) FROM pg_temp.t1_current_data) <> 3831 THEN RAISE EXCEPTION 'Unexpected collected cohort count'; END IF;
END $$;
INSERT INTO t1_features
SELECT d.unit_id,'lab',COALESCE(pg_temp.t1_lab(r#>>'{dictionary,ref_value}',r->>'lab_id','sicdb'),'native_'||COALESCE(r#>>'{dictionary,ref_value}',r->>'lab_id','unknown')),r->>'lab_id',r#>>'{dictionary,ref_value}',r#>>'{dictionary,ref_unit}',concat_ws(':',r->>'lab_type',r#>>'{specimen_dictionary,ref_value}'),
 pg_temp.t1_num(r->>'lab_value'),1,'arithmetic_mean_of_stored_records_0_24h'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'labs_24h','[]'::jsonb)) r;
INSERT INTO t1_features
SELECT d.unit_id,'vital',CASE r->>'data_id'
 WHEN '707' THEN 'heart_rate_ecg' WHEN '708' THEN 'heart_rate_spo2' WHEN '724' THEN 'heart_rate_abp'
 WHEN '703' THEN 'map_arterial' WHEN '706' THEN 'map_noninvasive' WHEN '3139' THEN 'sofa'
 ELSE r#>>'{dictionary,ref_value}' END,
 r->>'data_id',r#>>'{dictionary,ref_value}',r#>>'{dictionary,ref_unit}',NULL,
 pg_temp.t1_num(r->>'val'),1,'arithmetic_mean_of_hourly_source_values_0_24h'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'hourly_signals_24h','[]'::jsonb)) r
WHERE r->>'data_id' IN ('707','708','724','703','706','3139')
 OR COALESCE(r#>>'{dictionary,ref_value}','') ~* '(heart.?rate|blood.?pressure|temperature|respirat.*rate|spo2|sao2|glasgow|(^|[^a-z])gcs|sofa|fio2)';
INSERT INTO t1_features
SELECT d.unit_id,'case_extra',e.key,e.key,e.key,'score_or_source_unit',NULL,pg_temp.t1_num(e.value),1,'unique_value_else_null'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'case_extra','[]'::jsonb)) r
CROSS JOIN LATERAL jsonb_each_text(r) e WHERE e.key IN ('weight_kg','height_cm','saps3');
INSERT INTO t1_features
SELECT d.unit_id,'history',COALESCE(r#>>'{field_dictionary,ref_value}',r->>'field_id'),
 r->>'field_id',r#>>'{field_dictionary,ref_value}','0=no;1=yes',NULL,
 CASE lower(btrim(COALESCE(r#>>'{value_dictionary,ref_value}','')))
 WHEN 'yes' THEN 1 WHEN 'ja' THEN 1 WHEN 'true' THEN 1 WHEN '1' THEN 1 WHEN 'present' THEN 1 WHEN 'positive' THEN 1
 WHEN 'no' THEN 0 WHEN 'nein' THEN 0 WHEN 'false' THEN 0 WHEN '0' THEN 0 WHEN 'absent' THEN 0 WHEN 'negative' THEN 0 END,
 1,'explicit_structured_yes_no_conflict_to_null'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'structured_history','[]'::jsonb)) r
WHERE concat_ws(' ',r#>>'{field_dictionary,ref_value}',r#>>'{field_dictionary,ref_name}') ~* '(precondition|comorbid|history)';
INSERT INTO t1_features
SELECT d.unit_id,'admission_anthropometry',pg_temp.t1_norm(e.key),e.key,e.key,
 CASE WHEN pg_temp.t1_norm(e.key)='weightonadmission' THEN 'kg' ELSE 'cm' END,NULL,
 pg_temp.t1_num(e.value),1,'unique_value_else_null'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_each(d.collected_data) a
CROSS JOIN LATERAL jsonb_array_elements(a.value) r CROSS JOIN LATERAL jsonb_each_text(r) e
WHERE a.key LIKE 'raw_cases:%' AND pg_temp.t1_norm(e.key) IN ('weightonadmission','heightonadmission');
DROP TABLE IF EXISTS table1_direct_20261016.sicdb_history_flat_v2;
CREATE TABLE table1_direct_20261016.sicdb_history_flat_v2 AS
SELECT DISTINCT d.unit_id,r->>'field_id' field_id,r#>>'{field_dictionary,ref_value}' field_label,
 r#>>'{field_dictionary,ref_name}' field_group,r->>'ref_id' value_id,
 r#>>'{value_dictionary,ref_value}' value_label
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'structured_history','[]'::jsonb)) r;

DROP TABLE IF EXISTS pg_temp.t1_agg;
CREATE TEMP TABLE t1_agg AS
SELECT unit_id,domain,variable,source_key,source_label,native_unit,specimen,method,
 CASE WHEN method IN ('unique_value_else_null','explicit_structured_yes_no_conflict_to_null')
 THEN CASE WHEN count(DISTINCT value)=1 THEN min(value) END
 ELSE sum(value*weight) FILTER(WHERE value IS NOT NULL AND weight>0) /
 NULLIF(sum(weight) FILTER(WHERE value IS NOT NULL AND weight>0),0) END AS feature_value,
 COALESCE(sum(weight) FILTER(WHERE value IS NOT NULL AND weight>0),0) AS n_measurements,
 count(DISTINCT value)>1 AND method IN ('unique_value_else_null','explicit_structured_yes_no_conflict_to_null') AS conflicting_values,
 left(domain,18)||'_'||left(COALESCE(NULLIF(pg_temp.t1_norm(variable),''),'unknown'),20)||'_'||
 left(COALESCE(NULLIF(pg_temp.t1_norm(native_unit),''),'unitna'),8)||'_'||
 substr(md5(jsonb_build_array(domain,source_key,source_label,native_unit,specimen,method)::text),1,10) AS column_name
FROM t1_features GROUP BY unit_id,domain,variable,source_key,source_label,native_unit,specimen,method;
DROP TABLE IF EXISTS table1_direct_20261016.sicdb_field_map_v2;
CREATE TABLE table1_direct_20261016.sicdb_field_map_v2 AS
SELECT column_name,domain,variable,source_key,source_label,native_unit,specimen,method,
 variable IS NOT NULL AS included_in_wide,
 count(*) AS patients_with_source_records,count(feature_value) AS patients_with_numeric_value,
 sum(n_measurements) AS numeric_measurements,count(*) FILTER(WHERE conflicting_values) AS patients_with_conflicting_values
FROM t1_agg GROUP BY column_name,domain,variable,source_key,source_label,native_unit,specimen,method;
DO $$ BEGIN
 IF EXISTS(SELECT column_name FROM table1_direct_20261016.sicdb_field_map_v2 GROUP BY column_name HAVING count(*)>1)
 THEN RAISE EXCEPTION 'Feature column naming collision; report before exporting'; END IF;
END $$;
DO $pivot$
DECLARE feature_select text; k integer;
BEGIN
 SELECT count(*) INTO k FROM table1_direct_20261016.sicdb_field_map_v2 WHERE included_in_wide;
 IF k+73>1500 THEN RAISE EXCEPTION 'Too many columns (%). Stop and report; do not silently truncate.',k; END IF;
 SELECT string_agg(format('max(a.feature_value) FILTER (WHERE a.column_name=%L) AS %I',column_name,column_name),',' ORDER BY domain,variable,column_name)
 INTO feature_select FROM table1_direct_20261016.sicdb_field_map_v2 WHERE included_in_wide;
 DROP TABLE IF EXISTS table1_direct_20261016.sicdb_wide_v2;
 IF feature_select IS NULL THEN RAISE EXCEPTION 'No mapped features; report before export'; END IF;
 EXECUTE 'CREATE TEMP TABLE t1_pivot_values AS SELECT unit_id,'||feature_select||' FROM t1_agg a WHERE variable IS NOT NULL GROUP BY unit_id';
 EXECUTE $base$CREATE TABLE table1_direct_20261016.sicdb_wide_v2 AS
 SELECT d.unit_id AS cohort_unit_id,
d.assigned_class,
d.sepsis_onset_timestamp AS onset_timestamp_iso,
d.sepsis_onset_offset AS onset_offset,
d.onset_offset_unit,
pg_temp.t1_num(d.frozen_baseline->>'patient_id') AS "patient_id",
pg_temp.t1_num(d.frozen_baseline->>'case_id') AS "case_id",
pg_temp.t1_num(d.frozen_baseline->>'sepsis_offset_sec') AS "sepsis_offset_sec",
pg_temp.t1_num(d.frozen_baseline->>'landmark_offset_sec') AS "landmark_offset_sec",
pg_temp.t1_num(d.frozen_baseline->>'msi_m1') AS "msi_m1",
pg_temp.t1_num(d.frozen_baseline->>'msi_m2') AS "msi_m2",
pg_temp.t1_num(d.frozen_baseline->>'msi_m3') AS "msi_m3",
pg_temp.t1_num(d.frozen_baseline->>'msi_m4') AS "msi_m4",
pg_temp.t1_num(d.frozen_baseline->>'msi_m5') AS "msi_m5",
pg_temp.t1_num(d.frozen_baseline->>'msi_m6') AS "msi_m6",
pg_temp.t1_num(d.frozen_baseline->>'msi_m7') AS "msi_m7",
pg_temp.t1_num(d.frozen_baseline->>'msi_m8') AS "msi_m8",
pg_temp.t1_num(d.frozen_baseline->>'msi_m9') AS "msi_m9",
pg_temp.t1_num(d.frozen_baseline->>'msi_m10') AS "msi_m10",
pg_temp.t1_num(d.frozen_baseline->>'msi_m11') AS "msi_m11",
pg_temp.t1_num(d.frozen_baseline->>'msi_m12') AS "msi_m12",
pg_temp.t1_num(d.frozen_baseline->>'sodium_t1') AS "sodium_t1",
pg_temp.t1_num(d.frozen_baseline->>'sodium_t2') AS "sodium_t2",
pg_temp.t1_num(d.frozen_baseline->>'sodium_t3') AS "sodium_t3",
pg_temp.t1_num(d.frozen_baseline->>'sodium_t4') AS "sodium_t4",
pg_temp.t1_num(d.frozen_baseline->>'potassium_t1') AS "potassium_t1",
pg_temp.t1_num(d.frozen_baseline->>'potassium_t2') AS "potassium_t2",
pg_temp.t1_num(d.frozen_baseline->>'potassium_t3') AS "potassium_t3",
pg_temp.t1_num(d.frozen_baseline->>'potassium_t4') AS "potassium_t4",
pg_temp.t1_num(d.frozen_baseline->>'chloride_t1') AS "chloride_t1",
pg_temp.t1_num(d.frozen_baseline->>'chloride_t2') AS "chloride_t2",
pg_temp.t1_num(d.frozen_baseline->>'chloride_t3') AS "chloride_t3",
pg_temp.t1_num(d.frozen_baseline->>'chloride_t4') AS "chloride_t4",
pg_temp.t1_num(d.frozen_baseline->>'bicarbonate_t1') AS "bicarbonate_t1",
pg_temp.t1_num(d.frozen_baseline->>'bicarbonate_t2') AS "bicarbonate_t2",
pg_temp.t1_num(d.frozen_baseline->>'bicarbonate_t3') AS "bicarbonate_t3",
pg_temp.t1_num(d.frozen_baseline->>'bicarbonate_t4') AS "bicarbonate_t4",
pg_temp.t1_num(d.frozen_baseline->>'age') AS "age",
d.frozen_baseline->>'sex' AS "sex",
d.frozen_baseline->>'race_ethnicity' AS "race_ethnicity",
pg_temp.t1_num(d.frozen_baseline->>'weight_kg') AS "weight_kg",
pg_temp.t1_num(d.frozen_baseline->>'height_cm') AS "height_cm",
pg_temp.t1_num(d.frozen_baseline->>'bmi') AS "bmi",
pg_temp.t1_num(d.frozen_baseline->>'sepsis_sofa_score') AS "sepsis_sofa_score",
pg_temp.t1_num(d.frozen_baseline->>'saps3') AS "saps3",
pg_temp.t1_num(d.frozen_baseline->>'historical_creatinine_median') AS "historical_creatinine_median",
pg_temp.t1_num(d.frozen_baseline->>'first_admission_creatinine') AS "first_admission_creatinine",
pg_temp.t1_num(d.frozen_baseline->>'mdrd_backcalc_creatinine') AS "mdrd_backcalc_creatinine",
d.frozen_baseline->>'baseline_creatinine_source' AS "baseline_creatinine_source",
pg_temp.t1_num(d.frozen_baseline->>'baseline_creatinine') AS "baseline_creatinine",
pg_temp.t1_num(d.frozen_baseline->>'baseline_egfr_mdrd') AS "baseline_egfr_mdrd",
pg_temp.t1_num(d.frozen_baseline->>'ckd_or_renal_disease') AS "ckd_or_renal_disease",
pg_temp.t1_num(d.frozen_baseline->>'diabetes') AS "diabetes",
pg_temp.t1_num(d.frozen_baseline->>'congestive_heart_failure') AS "congestive_heart_failure",
pg_temp.t1_num(d.frozen_baseline->>'liver_disease') AS "liver_disease",
pg_temp.t1_num(d.frozen_baseline->>'malignancy') AS "malignancy",
pg_temp.t1_num(d.frozen_baseline->>'mechanical_ventilation_0_24h') AS "mechanical_ventilation_0_24h",
pg_temp.t1_num(d.frozen_baseline->>'vasopressor_0_24h') AS "vasopressor_0_24h",
d.frozen_baseline->>'icd10_main' AS "icd10_main",
d.frozen_baseline->>'icd10_main_text' AS "icd10_main_text",
pg_temp.t1_num(d.frozen_baseline->>'make30_risk_start_offset_sec') AS "make30_risk_start_offset_sec",
pg_temp.t1_num(d.frozen_baseline->>'make30_censor_offset_sec') AS "make30_censor_offset_sec",
pg_temp.t1_num(d.frozen_baseline->>'hospital_discharge_offset_sec') AS "hospital_discharge_offset_sec",
pg_temp.t1_num(d.frozen_baseline->>'make30_followup_days_after_landmark') AS "make30_followup_days_after_landmark",
pg_temp.t1_num(d.frozen_baseline->>'make30_evaluable_window') AS "make30_evaluable_window",
pg_temp.t1_num(d.frozen_baseline->>'make30_death') AS "make30_death",
pg_temp.t1_num(d.frozen_baseline->>'prior_krt_before_landmark') AS "prior_krt_before_landmark",
pg_temp.t1_num(d.frozen_baseline->>'make30_new_krt') AS "make30_new_krt",
pg_temp.t1_num(d.frozen_baseline->>'first_new_krt_offset_sec') AS "first_new_krt_offset_sec",
pg_temp.t1_num(d.frozen_baseline->>'final_inpatient_creatinine') AS "final_inpatient_creatinine",
pg_temp.t1_num(d.frozen_baseline->>'final_inpatient_creatinine_offset_sec') AS "final_inpatient_creatinine_offset_sec",
pg_temp.t1_num(d.frozen_baseline->>'make30_persistent_kidney_dysfunction') AS "make30_persistent_kidney_dysfunction",
pg_temp.t1_num(d.frozen_baseline->>'make30') AS "make30",p.* FROM pg_temp.t1_current_data d LEFT JOIN t1_pivot_values p ON p.unit_id=d.unit_id$base$;
 ALTER TABLE table1_direct_20261016.sicdb_wide_v2 DROP COLUMN unit_id;
 ALTER TABLE table1_direct_20261016.sicdb_wide_v2 ADD PRIMARY KEY(cohort_unit_id);
 DROP TABLE t1_pivot_values;
 IF (SELECT count(*) FROM table1_direct_20261016.sicdb_wide_v2)<>3831 THEN RAISE EXCEPTION 'Wide table row count mismatch'; END IF;
END $pivot$;
COMMIT;
-- This is the ONLY result set. No JSON/array-valued columns are exported.
SELECT * FROM table1_direct_20261016.sicdb_wide_v2 ORDER BY cohort_unit_id;
