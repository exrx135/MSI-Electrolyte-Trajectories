-- Standalone PostgreSQL; native eICU tables only. Final 5839 patients.
-- Stable comorbidity ascertainment: same hospital admission, all linked ICU episodes.
CREATE TEMP TABLE IF NOT EXISTS cci_session(dummy integer);
DROP TABLE IF EXISTS pg_temp.cci_seed;
CREATE TEMP TABLE cci_seed(unit_id bigint PRIMARY KEY,hospital_id bigint,patient_id text,age numeric);
-- Public-code edition: patient rows are loaded from a private local table.
INSERT INTO cci_seed SELECT * FROM local_inputs.eicu_cci_seed;
DO $seed_check$ BEGIN IF NOT EXISTS (SELECT 1 FROM cci_seed) THEN
RAISE EXCEPTION 'Populate local_inputs.eicu_cci_seed with authorised local data first';
END IF; END $seed_check$;
DROP TABLE IF EXISTS pg_temp.cci_links;
CREATE TEMP TABLE cci_links AS
SELECT s.unit_id,p.patientunitstayid FROM cci_seed s
JOIN public.patient p ON p.patienthealthsystemstayid=s.hospital_id AND p.uniquepid::text=s.patient_id;
CREATE INDEX ON cci_links(patientunitstayid);
WITH dx AS (
 SELECT l.unit_id,count(*) AS diagnosis_rows,
 string_agg(DISTINCT NULLIF(trim(d.icd9code),''),' | ') AS diagnosis_codes,
 string_agg(DISTINCT NULLIF(trim(d.diagnosisstring),''),' | ') AS diagnosis_text
 FROM cci_links l JOIN public.diagnosis d USING(patientunitstayid) GROUP BY l.unit_id
), hx AS (
 SELECT l.unit_id,count(*) AS history_rows,
 string_agg(DISTINCT concat_ws(' :: ',h.pasthistorypath,h.pasthistoryvalue,h.pasthistoryvaluetext),' | ') AS history_text
 FROM cci_links l JOIN public.pasthistory h USING(patientunitstayid) GROUP BY l.unit_id
)
SELECT s.unit_id AS patientunitstayid,s.hospital_id AS patienthealthsystemstayid,s.age,
 COALESCE(dx.diagnosis_rows,0) AS diagnosis_rows,dx.diagnosis_codes,dx.diagnosis_text,
 COALESCE(hx.history_rows,0) AS history_rows,hx.history_text
FROM cci_seed s LEFT JOIN dx USING(unit_id) LEFT JOIN hx USING(unit_id) ORDER BY s.unit_id;
