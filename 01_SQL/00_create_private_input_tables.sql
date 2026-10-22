-- Empty local input tables. Populate ONLY in an authorised, private database.
CREATE SCHEMA IF NOT EXISTS local_inputs;

CREATE TABLE IF NOT EXISTS local_inputs.mimic_baseline_seed (
 unit_id bigint PRIMARY KEY, hospital_id bigint, patient_id text,
 onset_timestamp timestamp, onset_offset numeric,
 assigned_class integer, frozen_baseline jsonb NOT NULL
);

CREATE TABLE IF NOT EXISTS local_inputs.eicu_baseline_seed (
 unit_id bigint PRIMARY KEY, hospital_id bigint, patient_id text,
 onset_timestamp timestamp, onset_offset numeric,
 assigned_class integer, frozen_baseline jsonb NOT NULL
);

CREATE TABLE IF NOT EXISTS local_inputs.sicdb_baseline_seed (
 unit_id bigint PRIMARY KEY, hospital_id bigint, patient_id text,
 onset_timestamp timestamp, onset_offset numeric,
 assigned_class integer, frozen_baseline jsonb NOT NULL
);

CREATE TABLE IF NOT EXISTS local_inputs.sicdb_vitals_seed (unit_id bigint PRIMARY KEY,onset_offset numeric,assigned_class integer);

CREATE TABLE IF NOT EXISTS local_inputs.eicu_cci_seed (unit_id bigint PRIMARY KEY,hospital_id bigint,patient_id text,age numeric);

CREATE TABLE IF NOT EXISTS local_inputs.eicu_survival_review_ids (patientunitstayid bigint PRIMARY KEY);
