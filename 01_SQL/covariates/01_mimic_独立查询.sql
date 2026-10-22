-- STANDALONE: run only this file in the corresponding native database.
-- Requires privately populated local_inputs.mimic_baseline_seed; no patient records are distributed.
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
INSERT INTO t1_seed SELECT * FROM local_inputs.mimic_baseline_seed;
DO $seed_check$ BEGIN IF NOT EXISTS (SELECT 1 FROM t1_seed) THEN
RAISE EXCEPTION 'Populate local_inputs.mimic_baseline_seed with authorised local data first';
END IF; END $seed_check$;
ANALYZE t1_seed;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('identity', ARRAY['mimiciv_icu.icustays','mimiciv_hosp.admissions'], $query$
SELECT s.unit_id,jsonb_build_object('icu',to_jsonb(i),'hospital',to_jsonb(a)) payload
FROM t1_seed s JOIN mimiciv_icu.icustays i ON i.stay_id=s.unit_id AND i.hadm_id=s.hospital_id AND i.subject_id::text=s.patient_id
JOIN mimiciv_hosp.admissions a ON a.hadm_id=s.hospital_id
$query$, true); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('labs_24h', ARRAY['mimiciv_hosp.labevents','mimiciv_hosp.d_labitems'], $query$
SELECT s.unit_id,to_jsonb(l)||jsonb_build_object('dictionary',to_jsonb(d)) payload
FROM t1_seed s JOIN mimiciv_hosp.labevents l ON l.hadm_id=s.hospital_id
LEFT JOIN mimiciv_hosp.d_labitems d ON d.itemid=l.itemid
WHERE l.charttime>=s.onset_timestamp AND l.charttime<s.onset_timestamp+interval '24 hours'
$query$, true); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('vitals_gcs_fio2_24h', ARRAY['mimiciv_icu.chartevents','mimiciv_icu.d_items'], $query$
WITH x AS (SELECT s.unit_id,c.itemid::text||':'||COALESCE(d.label,'') variable,c.valueuom::text unit,c.valuenum::numeric value
FROM t1_seed s JOIN mimiciv_icu.chartevents c ON c.hadm_id=s.hospital_id
LEFT JOIN mimiciv_icu.d_items d ON d.itemid=c.itemid
WHERE c.itemid IN (220045,220050,220051,220052,220179,220180,220181,225312,220210,224688,224689,224690,223761,223762,220277,223835,220739,223900,223901,226755,227013) AND c.charttime>=s.onset_timestamp AND c.charttime<s.onset_timestamp+interval '24 hours' ), h AS (
 SELECT unit_id,variable,unit,value,count(*) AS n FROM x
 WHERE value IS NOT NULL GROUP BY unit_id,variable,unit,value
)
SELECT unit_id,jsonb_build_object('variable',variable,'unit',unit,
 'mean_raw',sum(value*n)/sum(n),'n',sum(n),'min_raw',min(value),'max_raw',max(value),
 'value_counts',jsonb_agg(jsonb_build_array(value,n))) AS payload
FROM h GROUP BY unit_id,variable,unit
$query$, true); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('diagnoses_hospital', ARRAY['mimiciv_hosp.diagnoses_icd','mimiciv_hosp.d_icd_diagnoses'], $query$
SELECT s.unit_id,to_jsonb(d)||jsonb_build_object('dictionary',to_jsonb(n)) payload
FROM t1_seed s JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id=s.hospital_id
LEFT JOIN mimiciv_hosp.d_icd_diagnoses n ON n.icd_code=d.icd_code AND n.icd_version=d.icd_version
$query$, true); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('sofa_24h', ARRAY['mimiciv_derived.sofa'], $query$
SELECT s.unit_id,to_jsonb(d) payload FROM t1_seed s JOIN mimiciv_derived.sofa d ON d.stay_id=s.unit_id WHERE d.endtime>=s.onset_timestamp AND d.endtime<s.onset_timestamp+interval '24 hours'
$query$, false); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('gcs_24h', ARRAY['mimiciv_derived.gcs'], $query$
SELECT s.unit_id,to_jsonb(d) payload FROM t1_seed s JOIN mimiciv_derived.gcs d ON d.stay_id=s.unit_id WHERE d.charttime>=s.onset_timestamp AND d.charttime<s.onset_timestamp+interval '24 hours'
$query$, false); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('blood_gas_24h', ARRAY['mimiciv_derived.bg'], $query$
SELECT s.unit_id,to_jsonb(d) payload FROM t1_seed s JOIN mimiciv_derived.bg d ON d.hadm_id=s.hospital_id WHERE d.charttime>=s.onset_timestamp AND d.charttime<s.onset_timestamp+interval '24 hours'
$query$, false); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('charlson_hospital', ARRAY['mimiciv_derived.charlson'], $query$
SELECT s.unit_id,to_jsonb(d) payload FROM t1_seed s JOIN mimiciv_derived.charlson d ON d.hadm_id=s.hospital_id
$query$, false); END $collect$;

DO $collect$ BEGIN PERFORM pg_temp.t1_collect('apsiii_admission', ARRAY['mimiciv_derived.apsiii'], $query$
SELECT s.unit_id,to_jsonb(d) payload FROM t1_seed s JOIN mimiciv_derived.apsiii d ON d.stay_id=s.unit_id
$query$, false); END $collect$;
DROP TABLE IF EXISTS pg_temp.t1_current_data;
CREATE TEMP TABLE t1_current_data AS SELECT s.unit_id,s.hospital_id,s.patient_id,
s.onset_timestamp AS sepsis_onset_timestamp,s.onset_offset AS sepsis_onset_offset,
'timestamp'::text AS onset_offset_unit,
s.assigned_class,s.frozen_baseline,COALESCE(p.data,'{}'::jsonb) AS collected_data
FROM t1_seed s LEFT JOIN (SELECT unit_id,jsonb_object_agg(source,records) data FROM t1_payload GROUP BY unit_id) p USING(unit_id);
CREATE UNIQUE INDEX ON t1_current_data(unit_id);
CREATE SCHEMA IF NOT EXISTS table1_direct_20261016;
DROP TABLE IF EXISTS table1_direct_20261016.mimic_source_audit;
CREATE TABLE table1_direct_20261016.mimic_source_audit AS SELECT * FROM t1_audit;

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
 IF (SELECT count(*) FROM pg_temp.t1_current_data) <> 9804 THEN RAISE EXCEPTION 'Unexpected collected cohort count'; END IF;
END $$;
INSERT INTO t1_features
SELECT d.unit_id,'lab',COALESCE(pg_temp.t1_lab(r#>>'{dictionary,label}',r->>'itemid','mimic'),'native_'||COALESCE(r#>>'{dictionary,label}',r->>'itemid','unknown')),r->>'itemid',r#>>'{dictionary,label}',r->>'valueuom',r#>>'{dictionary,fluid}',
 pg_temp.t1_num(r->>'valuenum'),1,'arithmetic_mean_of_stored_records_0_24h'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'labs_24h','[]'::jsonb)) r;
INSERT INTO t1_features
SELECT d.unit_id,'vital',split_part(r->>'variable',':',2),r->>'variable',r->>'variable',r->>'unit',NULL,
 pg_temp.t1_num(r->>'mean_raw'),pg_temp.t1_num(r->>'n'),'record_count_weighted_mean_0_24h'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'vitals_gcs_fio2_24h','[]'::jsonb)) r;
UPDATE t1_features SET variable=COALESCE(NULLIF(variable,''),source_key) WHERE domain='vital';
INSERT INTO t1_features
SELECT d.unit_id,'sofa_24h',e.key,e.key,e.key,'score_or_source_unit',NULL,pg_temp.t1_num(e.value),1,'arithmetic_mean_of_stored_scores_0_24h'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'sofa_24h','[]'::jsonb)) r
CROSS JOIN LATERAL jsonb_each_text(r) e WHERE e.key IN ('sofa_24hours','respiration_24hours','coagulation_24hours','liver_24hours','cardiovascular_24hours','cns_24hours','renal_24hours');
INSERT INTO t1_features
SELECT d.unit_id,'gcs_24h',e.key,e.key,e.key,'score_or_source_unit',NULL,pg_temp.t1_num(e.value),1,'arithmetic_mean_of_stored_scores_0_24h'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'gcs_24h','[]'::jsonb)) r
CROSS JOIN LATERAL jsonb_each_text(r) e WHERE e.key IN ('gcs','gcs_motor','gcs_verbal','gcs_eyes');
INSERT INTO t1_features
SELECT d.unit_id,'blood_gas_24h',e.key,e.key,e.key,'score_or_source_unit',NULL,pg_temp.t1_num(e.value),1,'arithmetic_mean_of_stored_scores_0_24h'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'blood_gas_24h','[]'::jsonb)) r
CROSS JOIN LATERAL jsonb_each_text(r) e WHERE e.key IN ('pao2fio2ratio');
INSERT INTO t1_features
SELECT d.unit_id,'apsiii_admission',e.key,e.key,e.key,'score_or_source_unit',NULL,pg_temp.t1_num(e.value),1,'unique_value_else_null'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'apsiii_admission','[]'::jsonb)) r
CROSS JOIN LATERAL jsonb_each_text(r) e WHERE e.key IN ('apsiii');
INSERT INTO t1_features
SELECT d.unit_id,'charlson_hospital',e.key,e.key,e.key,'score_or_source_unit',NULL,pg_temp.t1_num(e.value),1,'unique_value_else_null'
FROM pg_temp.t1_current_data d CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'charlson_hospital','[]'::jsonb)) r
CROSS JOIN LATERAL jsonb_each_text(r) e WHERE e.key IN ('myocardial_infarct','congestive_heart_failure','peripheral_vascular_disease','cerebrovascular_disease','dementia','chronic_pulmonary_disease','rheumatic_disease','peptic_ulcer_disease','mild_liver_disease','diabetes_without_cc','diabetes_with_cc','paraplegia','renal_disease','malignant_cancer','severe_liver_disease','metastatic_solid_tumor','aids','charlson_comorbidity_index');
WITH diagnoses AS (
 SELECT d.unit_id,r,replace(upper(r->>'icd_code'),'.','') code FROM pg_temp.t1_current_data d
 CROSS JOIN LATERAL jsonb_array_elements(COALESCE(d.collected_data->'diagnoses_hospital','[]'::jsonb)) r
), coded AS (SELECT * FROM diagnoses WHERE NULLIF(code,'') IS NOT NULL)
INSERT INTO t1_features
SELECT unit_id,'diagnosis','hypertension_icd_recorded','hypertension_icd_recorded','hypertension ICD code recorded','0=not coded;1=coded',NULL,
 max(CASE WHEN ((r->>'icd_version'='9' AND code ~ '^40[1-5]') OR (r->>'icd_version'='10' AND code ~ '^I1[0-5]')) THEN 1 ELSE 0 END),1,'any_matching_hospital_icd_code' FROM coded GROUP BY unit_id
UNION ALL
SELECT unit_id,'diagnosis','ischemic_heart_disease_icd_recorded','ischemic_heart_disease_icd_recorded','ischemic heart disease ICD code recorded','0=not coded;1=coded',NULL,
 max(CASE WHEN ((r->>'icd_version'='9' AND code ~ '^41[0-4]') OR (r->>'icd_version'='10' AND code ~ '^I2[0-5]')) THEN 1 ELSE 0 END),1,'any_matching_hospital_icd_code' FROM coded GROUP BY unit_id;

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
DROP TABLE IF EXISTS table1_direct_20261016.mimic_field_map_v2;
CREATE TABLE table1_direct_20261016.mimic_field_map_v2 AS
SELECT column_name,domain,variable,source_key,source_label,native_unit,specimen,method,
 variable IS NOT NULL AS included_in_wide,
 count(*) AS patients_with_source_records,count(feature_value) AS patients_with_numeric_value,
 sum(n_measurements) AS numeric_measurements,count(*) FILTER(WHERE conflicting_values) AS patients_with_conflicting_values
FROM t1_agg GROUP BY column_name,domain,variable,source_key,source_label,native_unit,specimen,method;
DO $$ BEGIN
 IF EXISTS(SELECT column_name FROM table1_direct_20261016.mimic_field_map_v2 GROUP BY column_name HAVING count(*)>1)
 THEN RAISE EXCEPTION 'Feature column naming collision; report before exporting'; END IF;
END $$;
DO $pivot$
DECLARE feature_select text; k integer;
BEGIN
 SELECT count(*) INTO k FROM table1_direct_20261016.mimic_field_map_v2 WHERE included_in_wide;
 IF k+74>1500 THEN RAISE EXCEPTION 'Too many columns (%). Stop and report; do not silently truncate.',k; END IF;
 SELECT string_agg(format('max(a.feature_value) FILTER (WHERE a.column_name=%L) AS %I',column_name,column_name),',' ORDER BY domain,variable,column_name)
 INTO feature_select FROM table1_direct_20261016.mimic_field_map_v2 WHERE included_in_wide;
 DROP TABLE IF EXISTS table1_direct_20261016.mimic_wide_v2;
 IF feature_select IS NULL THEN RAISE EXCEPTION 'No mapped features; report before export'; END IF;
 EXECUTE 'CREATE TEMP TABLE t1_pivot_values AS SELECT unit_id,'||feature_select||' FROM t1_agg a WHERE variable IS NOT NULL GROUP BY unit_id';
 EXECUTE $base$CREATE TABLE table1_direct_20261016.mimic_wide_v2 AS
 SELECT d.unit_id AS cohort_unit_id,
d.assigned_class,
d.sepsis_onset_timestamp AS onset_timestamp_iso,
d.sepsis_onset_offset AS onset_offset,
d.onset_offset_unit,
pg_temp.t1_num(d.frozen_baseline->>'subject_id') AS "subject_id",
pg_temp.t1_num(d.frozen_baseline->>'hadm_id') AS "hadm_id",
pg_temp.t1_num(d.frozen_baseline->>'stay_id') AS "stay_id",
d.frozen_baseline->>'sepsis_time' AS "sepsis_time",
d.frozen_baseline->>'landmark_time' AS "landmark_time",
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
pg_temp.t1_num(d.frozen_baseline->>'historical_creatinine_median') AS "historical_creatinine_median",
pg_temp.t1_num(d.frozen_baseline->>'first_admission_creatinine') AS "first_admission_creatinine",
pg_temp.t1_num(d.frozen_baseline->>'mdrd_backcalc_creatinine') AS "mdrd_backcalc_creatinine",
d.frozen_baseline->>'baseline_creatinine_source' AS "baseline_creatinine_source",
pg_temp.t1_num(d.frozen_baseline->>'baseline_creatinine') AS "baseline_creatinine",
pg_temp.t1_num(d.frozen_baseline->>'baseline_egfr_mdrd') AS "baseline_egfr_mdrd",
pg_temp.t1_num(d.frozen_baseline->>'ckd_or_renal_disease') AS "ckd_or_renal_disease",
pg_temp.t1_num(d.frozen_baseline->>'diabetes_without_complication') AS "diabetes_without_complication",
pg_temp.t1_num(d.frozen_baseline->>'diabetes_with_complication') AS "diabetes_with_complication",
pg_temp.t1_num(d.frozen_baseline->>'diabetes') AS "diabetes",
pg_temp.t1_num(d.frozen_baseline->>'congestive_heart_failure') AS "congestive_heart_failure",
pg_temp.t1_num(d.frozen_baseline->>'liver_disease') AS "liver_disease",
pg_temp.t1_num(d.frozen_baseline->>'malignancy') AS "malignancy",
pg_temp.t1_num(d.frozen_baseline->>'charlson_comorbidity_index') AS "charlson_comorbidity_index",
pg_temp.t1_num(d.frozen_baseline->>'mechanical_ventilation_0_24h') AS "mechanical_ventilation_0_24h",
pg_temp.t1_num(d.frozen_baseline->>'vasopressor_0_24h') AS "vasopressor_0_24h",
d.frozen_baseline->>'make30_risk_start' AS "make30_risk_start",
d.frozen_baseline->>'make30_censor_time' AS "make30_censor_time",
pg_temp.t1_num(d.frozen_baseline->>'make30_followup_days_after_landmark') AS "make30_followup_days_after_landmark",
pg_temp.t1_num(d.frozen_baseline->>'make30_evaluable_window') AS "make30_evaluable_window",
pg_temp.t1_num(d.frozen_baseline->>'make30_death') AS "make30_death",
d.frozen_baseline->>'death_time' AS "death_time",
pg_temp.t1_num(d.frozen_baseline->>'prior_krt_before_landmark') AS "prior_krt_before_landmark",
pg_temp.t1_num(d.frozen_baseline->>'make30_new_krt') AS "make30_new_krt",
d.frozen_baseline->>'first_new_krt_time' AS "first_new_krt_time",
pg_temp.t1_num(d.frozen_baseline->>'final_inpatient_creatinine') AS "final_inpatient_creatinine",
d.frozen_baseline->>'final_inpatient_creatinine_time' AS "final_inpatient_creatinine_time",
pg_temp.t1_num(d.frozen_baseline->>'make30_persistent_kidney_dysfunction') AS "make30_persistent_kidney_dysfunction",
pg_temp.t1_num(d.frozen_baseline->>'make30') AS "make30",p.* FROM pg_temp.t1_current_data d LEFT JOIN t1_pivot_values p ON p.unit_id=d.unit_id$base$;
 ALTER TABLE table1_direct_20261016.mimic_wide_v2 DROP COLUMN unit_id;
 ALTER TABLE table1_direct_20261016.mimic_wide_v2 ADD PRIMARY KEY(cohort_unit_id);
 DROP TABLE t1_pivot_values;
 IF (SELECT count(*) FROM table1_direct_20261016.mimic_wide_v2)<>9804 THEN RAISE EXCEPTION 'Wide table row count mismatch'; END IF;
END $pivot$;
COMMIT;
-- This is the ONLY result set. No JSON/array-valued columns are exported.
SELECT * FROM table1_direct_20261016.mimic_wide_v2 ORDER BY cohort_unit_id;
