/* FINAL ANALYSIS EXTENSION: preserves frozen trajectory extraction and adds ORAKLE MAKE30 + covariates.
   Outcome events are counted only after Sepsis t0+24h landmark.
   Generated 2026-08-27. */

/* ========================================================================
   MSI STUDY — MIMIC-IV v3.1
   FULL-COHORT FORMAL RESEARCH PIPELINE — MIXED-FREQUENCY CORE-4 VERSION: MSI 2-H + ELECTROLYTES 6-H (CUMULATIVE, DIRECTLY RUNNABLE)

   PURPOSE
   -------
   This is NOT a database-only audit.
   This script runs the real study pipeline on a full patient cohort:

     1) full cohort of ICU patients
     2) ORAKLE-style Sepsis identification
     3) first qualifying Sepsis ICU stay per patient
     4) +24 h landmark requirement
     5) exclude active/typed KRT during Sepsis 0–24 h
     6) extract MSI in twelve 2-h bins and electrolytes in four 6-h bins
        (all 7 are still extracted so the reduced model can be compared fairly)
     7) calculate patient-level longitudinal coverage
     8) define the PRIMARY reduced trajectory:
          MSI + Na + K + Cl + HCO3
     9) retain the old 7-electrolyte eligibility only as a comparison
    10) create the actual Core-4 GBMTM input table
    11) output cohort flow + coverage + bin-level missingness + final-cohort value QC

   IMPORTANT
   ---------
   - This is the formal full-cohort extraction and analysis pipeline.
   - All downstream logic is the intended formal research logic.
   - Exposure window is [sepsis_time, sepsis_time + 24 h).
   - Landmark is exactly sepsis_time + 24 h.
   - MSI is formed by:
       a) hourly median HR
       b) hourly invasive MAP and NIBP summaries
       c) MAP QC: both invasive MAP and NIBP <=20 mmHg are excluded before hourly summaries;
          if both valid hourly sources exist and |invasive-NIBP| >=30 mmHg, use NIBP;
          otherwise valid invasive MAP is preferred
       d) hourly MSI = HR / selected MAP
       e) each 2-h MSI = maximum valid hourly MSI in that 2-h bin
   - Electrolytes:
       median value in each 6-h bin
   - Missing time-bin values remain NULL; no LOCF/NOCB/interpolation.

   PRIMARY CORE-4 TRAJECTORY ELIGIBILITY
   -------------------------------------
   MSI: >=9 distinct 2-h bins (formal 75% coverage threshold; 6/12 and 12/12 also reported)
   AND
   Na/K/Cl/HCO3: each >=2 distinct 6-h bins

   Electrolyte hard-QC before 6-h aggregation:
   - sodium >300 mmol/L is discarded; sodium =300 is retained.
   - chloride <10 mmol/L is discarded; chloride =10 is retained.
   These rules are applied to raw measurements before the 6-h median.

   COMPARATOR (OLD ALL-7 RULE)
   ---------------------------
   MSI: >=9 distinct 2-h bins (formal 75% coverage threshold; 6/12 and 12/12 also reported)
   AND
   Na/K/Cl/HCO3/Ca/Mg/P: each >=2 distinct 6-h bins

   NOTE
   ----
   Ca/Mg/P are still extracted in this script only so that:
   - their missingness remains visible;
   - Core-4 versus All-7 cohort retention can be compared on the same patients.
   They are NOT included in the primary Core-4 GBMTM input table.

   MIMIC-IV vital itemids used:
     HR                  220045
     invasive MAP        220052, 225312
     non-invasive MAP    220181

   MIMIC-IV lab itemids used:
     Magnesium           50960
     Phosphate           50970
   ======================================================================== */


/* ========================================================================
   0. SCHEMA AND CLEANUP
   ======================================================================== */

CREATE SCHEMA IF NOT EXISTS study_make30_final;

DROP TABLE IF EXISTS study_make30_final.mimic_analysis_master;
DROP TABLE IF EXISTS study_make30_final.mimic_make30;
DROP TABLE IF EXISTS study_make30_final.mimic_covariates;
DROP TABLE IF EXISTS study_make30_final.mimic_baseline_kidney;
DROP TABLE IF EXISTS study_make30_final.mimic_qc_summary;


DROP TABLE IF EXISTS study_make30_final.mimic_msi_gt4_trace;
DROP TABLE IF EXISTS study_make30_final.mimic_core4_multitrajectory_wide;
DROP TABLE IF EXISTS study_make30_final.mimic_core4_electrolyte_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.mimic_core4_msi_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.mimic_core4_trajectory_cohort;
DROP TABLE IF EXISTS study_make30_final.mimic_trajectory_coverage;
DROP TABLE IF EXISTS study_make30_final.mimic_electrolyte_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.mimic_msi_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.mimic_electrolyte_6h;
DROP TABLE IF EXISTS study_make30_final.mimic_msi_2h;
DROP TABLE IF EXISTS study_make30_final.mimic_msi_hourly;
DROP TABLE IF EXISTS study_make30_final.mimic_trajectory_pool;
DROP TABLE IF EXISTS study_make30_final.mimic_landmark_cohort;
DROP TABLE IF EXISTS study_make30_final.mimic_sepsis_index;
DROP TABLE IF EXISTS study_make30_final.mimic_sepsis_adult;
DROP TABLE IF EXISTS study_make30_final.mimic_sepsis_first_stay;
DROP TABLE IF EXISTS study_make30_final.mimic_sepsis_pairs;
DROP TABLE IF EXISTS study_make30_final.mimic_sofa_ge2;
DROP TABLE IF EXISTS study_make30_final.mimic_icustays;
DROP TABLE IF EXISTS study_make30_final.mimic_patients;


/* ========================================================================
   1. FULL ICU PATIENT COHORT
      Unit = subject_id
      All ICU stays of all ICU patients are retained.
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_patients AS

SELECT DISTINCT
    subject_id

FROM mimiciv_icu.icustays

WHERE subject_id IS NOT NULL;


CREATE UNIQUE INDEX idx_mimic_patient
ON study_make30_final.mimic_patients(subject_id);


/* ========================================================================
   2. ALL ICU STAYS FOR THE FULL PATIENT COHORT
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_icustays AS

SELECT
    i.subject_id,
    i.hadm_id,
    i.stay_id,
    i.intime AS icu_intime,
    i.outtime AS icu_outtime,
    i.los AS icu_los

FROM mimiciv_icu.icustays i

INNER JOIN study_make30_final.mimic_patients p
    ON i.subject_id = p.subject_id;


CREATE UNIQUE INDEX idx_mimic_icustays_stay
ON study_make30_final.mimic_icustays(stay_id);

CREATE INDEX idx_mimic_icustays_hadm
ON study_make30_final.mimic_icustays(hadm_id);

CREATE INDEX idx_mimic_icustays_subject
ON study_make30_final.mimic_icustays(subject_id);


/* ========================================================================
   3. ORAKLE-STYLE SOFA COMPONENT
      Qualifying SOFA >= 2.
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_sofa_ge2 AS

SELECT
    s.stay_id,
    s.starttime,
    s.endtime AS sofa_time,
    s.sofa_24hours AS sofa_score,

    s.respiration_24hours AS sofa_respiration,
    s.coagulation_24hours AS sofa_coagulation,
    s.liver_24hours AS sofa_liver,
    s.cardiovascular_24hours AS sofa_cardiovascular,
    s.cns_24hours AS sofa_cns,
    s.renal_24hours AS sofa_renal

FROM mimiciv_derived.sofa s

INNER JOIN study_make30_final.mimic_icustays x
    ON s.stay_id = x.stay_id

WHERE s.sofa_24hours >= 2;


CREATE INDEX idx_mimic_sofa_stay_time
ON study_make30_final.mimic_sofa_ge2(stay_id, sofa_time);


/* ========================================================================
   4. ORAKLE-STYLE SUSPECTED INFECTION + SOFA MATCHING

      qualifying sofa_time:
        suspected infection -24 h
        through
        suspected infection +12 h

      sepsis_time:
        earlier of suspected infection time and qualifying SOFA time
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_sepsis_pairs AS

SELECT
    soi.subject_id,
    soi.hadm_id,
    soi.stay_id,

    soi.antibiotic,
    soi.antibiotic_time,
    soi.culture_time,
    soi.suspected_infection_time,

    s.sofa_time,
    s.sofa_score,

    s.sofa_respiration,
    s.sofa_coagulation,
    s.sofa_liver,
    s.sofa_cardiovascular,
    s.sofa_cns,
    s.sofa_renal,

    LEAST(
        soi.suspected_infection_time,
        s.sofa_time
    ) AS sepsis_time

FROM mimiciv_derived.suspicion_of_infection soi

INNER JOIN study_make30_final.mimic_icustays x
    ON soi.stay_id = x.stay_id

INNER JOIN study_make30_final.mimic_sofa_ge2 s
    ON soi.stay_id = s.stay_id

   AND s.sofa_time >=
       soi.suspected_infection_time - INTERVAL '24 hours'

   AND s.sofa_time <=
       soi.suspected_infection_time + INTERVAL '12 hours'

WHERE soi.suspected_infection = 1
  AND soi.suspected_infection_time IS NOT NULL;


CREATE INDEX idx_mimic_pairs_stay
ON study_make30_final.mimic_sepsis_pairs(stay_id);

CREATE INDEX idx_mimic_pairs_subject
ON study_make30_final.mimic_sepsis_pairs(subject_id);


/* ========================================================================
   5. FIRST QUALIFYING SEPSIS EVENT PER ICU STAY
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_sepsis_first_stay AS

SELECT
    subject_id,
    hadm_id,
    stay_id,

    antibiotic,
    antibiotic_time,
    culture_time,
    suspected_infection_time,

    sofa_time,
    sofa_score,

    sofa_respiration,
    sofa_coagulation,
    sofa_liver,
    sofa_cardiovascular,
    sofa_cns,
    sofa_renal,

    sepsis_time

FROM
(
    SELECT
        p.*,

        ROW_NUMBER() OVER
        (
            PARTITION BY p.stay_id

            ORDER BY
                p.sepsis_time,
                p.suspected_infection_time,
                p.sofa_time,
                p.antibiotic_time,
                p.culture_time
        ) AS rn_stay

    FROM study_make30_final.mimic_sepsis_pairs p
) q

WHERE rn_stay = 1;


/* ========================================================================
   6. ADD DEMOGRAPHICS / HOSPITAL / ICU INFORMATION
      Adult only.
      Sepsis onset must occur inside ICU.
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_sepsis_adult AS

SELECT
    s.subject_id,
    s.hadm_id,
    s.stay_id,

    d.gender AS sex,
    d.admission_age AS age,
    d.race,

    d.admittime,
    d.dischtime,

    d.icu_intime,
    d.icu_outtime,
    d.los_icu,

    s.antibiotic,
    s.antibiotic_time,
    s.culture_time,
    s.suspected_infection_time,

    s.sofa_time,
    s.sofa_score,

    s.sofa_respiration,
    s.sofa_coagulation,
    s.sofa_liver,
    s.sofa_cardiovascular,
    s.sofa_cns,
    s.sofa_renal,

    s.sepsis_time,

    s.sepsis_time + INTERVAL '24 hours'
        AS landmark_time,

    s.sepsis_time + INTERVAL '30 days'
        AS sepsis_day30,

    d.hospital_expire_flag

FROM study_make30_final.mimic_sepsis_first_stay s

INNER JOIN mimiciv_derived.icustay_detail d
    ON s.stay_id = d.stay_id

WHERE d.admission_age >= 18

  AND s.sepsis_time >= d.icu_intime

  AND s.sepsis_time < d.icu_outtime;


/* ========================================================================
   7. FIRST QUALIFYING SEPSIS ICU STAY PER PATIENT
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_sepsis_index AS

SELECT
    subject_id,
    hadm_id,
    stay_id,

    sex,
    age,
    race,

    admittime,
    dischtime,

    icu_intime,
    icu_outtime,
    los_icu,

    antibiotic,
    antibiotic_time,
    culture_time,
    suspected_infection_time,

    sofa_time,
    sofa_score,

    sofa_respiration,
    sofa_coagulation,
    sofa_liver,
    sofa_cardiovascular,
    sofa_cns,
    sofa_renal,

    sepsis_time,
    landmark_time,
    sepsis_day30,

    hospital_expire_flag

FROM
(
    SELECT
        a.*,

        ROW_NUMBER() OVER
        (
            PARTITION BY a.subject_id

            ORDER BY
                a.sepsis_time,
                a.icu_intime,
                a.stay_id
        ) AS rn_patient

    FROM study_make30_final.mimic_sepsis_adult a
) q

WHERE rn_patient = 1;


CREATE UNIQUE INDEX idx_mimic_sepsis_index_stay
ON study_make30_final.mimic_sepsis_index(stay_id);

CREATE UNIQUE INDEX idx_mimic_sepsis_index_subject
ON study_make30_final.mimic_sepsis_index(subject_id);


/* ========================================================================
   8. +24 H LANDMARK COHORT
      Must still be in ICU and still hospitalized at landmark.
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_landmark_cohort AS

SELECT
    s.*,

    CASE
        WHEN s.icu_outtime >= s.landmark_time
         AND s.dischtime >= s.landmark_time
        THEN 1
        ELSE 0
    END AS reached_24h_landmark

FROM study_make30_final.mimic_sepsis_index s;


CREATE UNIQUE INDEX idx_mimic_landmark_stay
ON study_make30_final.mimic_landmark_cohort(stay_id);


/* ========================================================================
   9. KRT DURING PHENOTYPE WINDOW [t0, t0+24 h)

      IMPORTANT:
      dialysis_present alone is NOT used for exclusion.
      Primary contamination signal:
        dialysis_active = 1
        OR a non-empty dialysis_type exists.

      At this stage this is deliberately a phenotype-window KRT exclusion.
      Formal post-landmark "new KRT" outcome will be constructed separately.
   ======================================================================== */

DROP TABLE IF EXISTS tmp_early_krt;

CREATE TEMP TABLE tmp_early_krt AS

SELECT
    s.stay_id,

    MAX(
        CASE
            WHEN r.dialysis_present = 1
            THEN 1
            ELSE 0
        END
    ) AS any_dialysis_signal_0_24h,

    MAX(
        CASE
            WHEN r.dialysis_active = 1
              OR NULLIF(BTRIM(r.dialysis_type), '') IS NOT NULL
            THEN 1
            ELSE 0
        END
    ) AS active_or_typed_krt_0_24h

FROM study_make30_final.mimic_landmark_cohort s

LEFT JOIN mimiciv_derived.rrt r
    ON r.stay_id = s.stay_id
   AND r.charttime >= s.sepsis_time
   AND r.charttime <  s.landmark_time

WHERE s.reached_24h_landmark = 1

GROUP BY s.stay_id;


CREATE INDEX idx_tmp_early_krt_stay
ON tmp_early_krt(stay_id);


/* ========================================================================
   10. TRAJECTORY POOL AFTER LANDMARK + KRT EXCLUSION
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_trajectory_pool AS

SELECT
    s.*,

    COALESCE(k.any_dialysis_signal_0_24h, 0)
        AS any_dialysis_signal_0_24h,

    COALESCE(k.active_or_typed_krt_0_24h, 0)
        AS active_or_typed_krt_0_24h

FROM study_make30_final.mimic_landmark_cohort s

LEFT JOIN tmp_early_krt k
    ON s.stay_id = k.stay_id

WHERE s.reached_24h_landmark = 1

  AND COALESCE(k.active_or_typed_krt_0_24h, 0) = 0;


CREATE UNIQUE INDEX idx_mimic_trajcand_stay
ON study_make30_final.mimic_trajectory_pool(stay_id);

CREATE INDEX idx_mimic_trajcand_hadm
ON study_make30_final.mimic_trajectory_pool(hadm_id);


/* ========================================================================
   11. RAW VITALS -> HOURLY HR / MAP

      HR:
        itemid 220045

      invasive MAP:
        220052
        225312

      non-invasive MAP:
        220181

      MIMIC MSI-input plausibility bounds:
        HR  >0 and <300
        MAP >20 and <300 (both invasive MAP and NIBP)

      MAP QC / source selection:
        both invasive MAP and NIBP <=20 mmHg are excluded before hourly summaries;
        if both valid hourly MAP summaries differ by >=30 mmHg, use NIBP;
        otherwise valid invasive MAP has priority
   ======================================================================== */

DROP TABLE IF EXISTS tmp_vital_hourly;

CREATE TEMP TABLE tmp_vital_hourly AS

SELECT
    b.stay_id,

    FLOOR(
        EXTRACT(
            EPOCH FROM (ce.charttime - b.sepsis_time)
        ) / 3600.0
    )::integer AS hour_bin,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY ce.valuenum)
    FILTER
    (
        WHERE ce.itemid = 220045
          AND ce.valuenum > 0
          AND ce.valuenum < 300
    ) AS hr_median,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY ce.valuenum)
    FILTER
    (
        WHERE ce.itemid IN (220052, 225312)
          AND ce.valuenum > 20
          AND ce.valuenum < 300
    ) AS map_invasive_median,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY ce.valuenum)
    FILTER
    (
        WHERE ce.itemid = 220181
          AND ce.valuenum > 20
          AND ce.valuenum < 300
    ) AS map_noninvasive_median

FROM study_make30_final.mimic_trajectory_pool b

INNER JOIN mimiciv_icu.chartevents ce
    ON ce.stay_id = b.stay_id
   AND ce.charttime >= b.sepsis_time
   AND ce.charttime <  b.landmark_time

WHERE ce.itemid IN
(
    220045,
    220052,
    225312,
    220181
)

GROUP BY
    b.stay_id,

    FLOOR(
        EXTRACT(
            EPOCH FROM (ce.charttime - b.sepsis_time)
        ) / 3600.0
    )::integer;


CREATE INDEX idx_tmp_vital_hourly
ON tmp_vital_hourly(stay_id, hour_bin);


/* ========================================================================
   12. VALID HOURLY MSI — UNIFIED MAP QC (FINAL <=20 RULE)

   HR:
      hourly median primary monitored HR.

   MAP preprocessing and selection order:
      1) BEFORE the hourly MAP summaries are formed, discard every
         invasive MAP value <=20 mmHg and every NIBP MAP value <=20 mmHg;
      2) derive hourly median invasive MAP and hourly median NIBP from the
         remaining valid source-specific measurements;
      3) if only one valid hourly MAP source exists, use that source;
      4) if both valid hourly sources exist:
           - if |invasive MAP - NIBP| >=30 mmHg, use NIBP;
           - otherwise use invasive MAP;
      5) if neither valid source exists, MAP is missing.

   MSI_hour = hourly HR summary / selected hourly MAP summary.
   No minimum raw-measurement count is imposed within an hour.

   NOTE:
      The <=20-mmHg exclusion is applied symmetrically to invasive MAP and
      NIBP, before source-priority and discordance rules are evaluated.
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_msi_hourly AS

WITH selected AS
(
    SELECT
        stay_id,
        hour_bin,
        hr_median,
        map_invasive_median,
        map_noninvasive_median,

        CASE
            WHEN map_invasive_median IS NOT NULL
             AND map_noninvasive_median IS NOT NULL
            THEN ABS(map_invasive_median - map_noninvasive_median)
            ELSE NULL
        END::double precision AS map_abs_diff,

        CASE
            WHEN map_invasive_median IS NOT NULL
             AND map_noninvasive_median IS NOT NULL
             AND ABS(map_invasive_median - map_noninvasive_median) >= 30
            THEN map_noninvasive_median

            WHEN map_invasive_median IS NOT NULL
            THEN map_invasive_median

            WHEN map_noninvasive_median IS NOT NULL
            THEN map_noninvasive_median

            ELSE NULL
        END::double precision AS map_used,

        CASE
            WHEN map_invasive_median IS NOT NULL
             AND map_noninvasive_median IS NOT NULL
             AND ABS(map_invasive_median - map_noninvasive_median) >= 30
            THEN 'noninvasive_fallback_discordant_ge30'

            WHEN map_invasive_median IS NOT NULL
            THEN 'invasive'

            WHEN map_noninvasive_median IS NOT NULL
            THEN 'noninvasive'

            ELSE NULL
        END::text AS map_source

    FROM tmp_vital_hourly

    WHERE hour_bin BETWEEN 0 AND 23
)

SELECT
    stay_id,
    hour_bin,
    hr_median,
    map_invasive_median,
    map_noninvasive_median,
    map_abs_diff,
    map_used,
    map_source,

    CASE
        WHEN hr_median IS NOT NULL
         AND map_used IS NOT NULL
         AND map_used > 0
        THEN hr_median / map_used
        ELSE NULL
    END::double precision AS msi

FROM selected;

CREATE INDEX idx_mimic_msi_hourly
ON study_make30_final.mimic_msi_hourly(stay_id, hour_bin);


/* ========================================================================
   13. HOURLY MSI -> TWELVE 2-H MSI WINDOWS

   Current feasibility definition:
      MSI_2h = MAX(valid hourly MSI) within each 2-h window.

   The MAX definition is deliberately retained from the prior pipeline so
   that only the temporal resolution changes in this test. Hourly MSI >4 is
   traced explicitly downstream for source-level QC.
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_msi_2h AS

SELECT
    stay_id,
    FLOOR(hour_bin / 2.0)::integer AS time_bin,
    MAX(msi)::double precision AS msi,
    COUNT(msi)::integer AS valid_hourly_msi_count

FROM study_make30_final.mimic_msi_hourly

WHERE msi IS NOT NULL

GROUP BY
    stay_id,
    FLOOR(hour_bin / 2.0)::integer;

CREATE UNIQUE INDEX idx_mimic_msi_2h
ON study_make30_final.mimic_msi_2h(stay_id, time_bin);

/* ========================================================================
   14. ELECTROLYTE RAW LONG TABLE

       Na / K / Cl / HCO3 / total Ca:
         mimiciv_derived.chemistry

       Mg:
         labevents itemid 50960

       Phosphorus/phosphate:
         labevents itemid 50970
   ======================================================================== */

DROP TABLE IF EXISTS tmp_electrolyte_raw;

CREATE TEMP TABLE tmp_electrolyte_raw
(
    stay_id          bigint,
    charttime        timestamp,
    variable_name    text,
    value_num        double precision
);


/* 14a. Na / K / Cl / HCO3 / total Ca */

INSERT INTO tmp_electrolyte_raw
(
    stay_id,
    charttime,
    variable_name,
    value_num
)

SELECT
    b.stay_id,
    c.charttime,
    x.variable_name,
    x.value_num

FROM study_make30_final.mimic_trajectory_pool b

INNER JOIN mimiciv_derived.chemistry c
    ON c.hadm_id = b.hadm_id
   AND c.charttime >= b.sepsis_time
   AND c.charttime <  b.landmark_time

CROSS JOIN LATERAL
(
    VALUES
        ('sodium'::text,      c.sodium::double precision),
        ('potassium'::text,   c.potassium::double precision),
        ('chloride'::text,    c.chloride::double precision),
        ('bicarbonate'::text, c.bicarbonate::double precision),
        ('calcium'::text,     c.calcium::double precision)
) AS x(variable_name, value_num)

WHERE x.value_num IS NOT NULL
  AND NOT (x.variable_name = 'sodium'   AND x.value_num > 300)
  AND NOT (x.variable_name = 'chloride' AND x.value_num < 10);


/* 14b. Magnesium + phosphorus */

INSERT INTO tmp_electrolyte_raw
(
    stay_id,
    charttime,
    variable_name,
    value_num
)

SELECT
    b.stay_id,
    l.charttime,

    CASE
        WHEN l.itemid = 50960
        THEN 'magnesium'

        WHEN l.itemid = 50970
        THEN 'phosphorus'
    END AS variable_name,

    l.valuenum::double precision

FROM study_make30_final.mimic_trajectory_pool b

INNER JOIN mimiciv_hosp.labevents l
    ON l.hadm_id = b.hadm_id
   AND l.charttime >= b.sepsis_time
   AND l.charttime <  b.landmark_time

WHERE l.itemid IN (50960, 50970)
  AND l.valuenum IS NOT NULL;


CREATE INDEX idx_tmp_electrolyte_raw
ON tmp_electrolyte_raw(stay_id, variable_name, charttime);


/* ========================================================================
   15. ELECTROLYTES -> FOUR 6-H BINS
       Value per bin = median of all available values in the bin
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_electrolyte_6h AS

SELECT
    e.stay_id,

    FLOOR(
        EXTRACT(
            EPOCH FROM (e.charttime - b.sepsis_time)
        ) / 21600.0
    )::integer AS time_bin,

    e.variable_name,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY e.value_num)
        AS median_value,

    COUNT(*) AS raw_measurement_count

FROM tmp_electrolyte_raw e

INNER JOIN study_make30_final.mimic_trajectory_pool b
    ON e.stay_id = b.stay_id

WHERE e.charttime >= b.sepsis_time
  AND e.charttime <  b.landmark_time

GROUP BY
    e.stay_id,

    FLOOR(
        EXTRACT(
            EPOCH FROM (e.charttime - b.sepsis_time)
        ) / 21600.0
    )::integer,

    e.variable_name;


CREATE UNIQUE INDEX idx_mimic_electrolyte_6h
ON study_make30_final.mimic_electrolyte_6h
(stay_id, time_bin, variable_name);



/* ========================================================================
   16. MSI 2-H TRAJECTORY LONG TABLE — 12 ROWS PER TRAJECTORY-POOL PATIENT
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_msi_trajectory_long AS

WITH patient_bins AS
(
    SELECT
        b.subject_id,
        b.hadm_id,
        b.stay_id,
        b.sepsis_time,
        b.landmark_time,
        gs.time_bin
    FROM study_make30_final.mimic_trajectory_pool b
    CROSS JOIN generate_series(0,11) AS gs(time_bin)
)
SELECT
    p.subject_id,
    p.hadm_id,
    p.stay_id,
    p.sepsis_time,
    p.landmark_time,
    p.time_bin,
    ('M' || (p.time_bin + 1)::text)::text AS time_label,
    (p.time_bin * 2)::integer AS window_start_hour,
    ((p.time_bin + 1) * 2)::integer AS window_end_hour,
    (p.time_bin * 2 + 1)::integer AS time_midpoint_hour,
    m.msi,
    m.valid_hourly_msi_count
FROM patient_bins p
LEFT JOIN study_make30_final.mimic_msi_2h m
    ON p.stay_id = m.stay_id
   AND p.time_bin = m.time_bin
ORDER BY p.stay_id, p.time_bin;

CREATE UNIQUE INDEX idx_mimic_msi_trajectory_long
ON study_make30_final.mimic_msi_trajectory_long(stay_id,time_bin);

/* ========================================================================
   17. ELECTROLYTE 6-H TRAJECTORY LONG TABLE — 4 ROWS PER TRAJECTORY-POOL PATIENT
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_electrolyte_trajectory_long AS

WITH patient_bins AS
(
    SELECT
        b.subject_id,
        b.hadm_id,
        b.stay_id,
        b.sepsis_time,
        b.landmark_time,
        gs.time_bin
    FROM study_make30_final.mimic_trajectory_pool b
    CROSS JOIN generate_series(0,3) AS gs(time_bin)
),
electrolyte_wide AS
(
    SELECT
        stay_id,
        time_bin,
        MAX(median_value) FILTER (WHERE variable_name='sodium') AS sodium,
        MAX(median_value) FILTER (WHERE variable_name='potassium') AS potassium,
        MAX(median_value) FILTER (WHERE variable_name='chloride') AS chloride,
        MAX(median_value) FILTER (WHERE variable_name='bicarbonate') AS bicarbonate,
        MAX(median_value) FILTER (WHERE variable_name='calcium') AS calcium,
        MAX(median_value) FILTER (WHERE variable_name='magnesium') AS magnesium,
        MAX(median_value) FILTER (WHERE variable_name='phosphorus') AS phosphorus
    FROM study_make30_final.mimic_electrolyte_6h
    GROUP BY stay_id,time_bin
)
SELECT
    p.subject_id,
    p.hadm_id,
    p.stay_id,
    p.sepsis_time,
    p.landmark_time,
    p.time_bin,
    ('T' || (p.time_bin + 1)::text)::text AS time_label,
    (p.time_bin * 6)::integer AS window_start_hour,
    ((p.time_bin + 1) * 6)::integer AS window_end_hour,
    (p.time_bin * 6 + 3)::integer AS time_midpoint_hour,
    e.sodium,
    e.potassium,
    e.chloride,
    e.bicarbonate,
    e.calcium,
    e.magnesium,
    e.phosphorus
FROM patient_bins p
LEFT JOIN electrolyte_wide e
    ON p.stay_id=e.stay_id AND p.time_bin=e.time_bin
ORDER BY p.stay_id,p.time_bin;

CREATE UNIQUE INDEX idx_mimic_electrolyte_trajectory_long
ON study_make30_final.mimic_electrolyte_trajectory_long(stay_id,time_bin);

/* ========================================================================
   18. PATIENT-LEVEL COVERAGE + FINAL PRIMARY ELIGIBILITY

   MSI primary threshold for this full-cohort analysis:
      >=9/12 2-h bins = 75% coverage, preserving the prior 3/4 requirement.

   Sensitivity counts for >=6/12 and 12/12 are emitted downstream.
   Electrolyte requirement remains unchanged: each >=2/4 6-h bins.
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_trajectory_coverage AS
WITH m AS
(
    SELECT subject_id,hadm_id,stay_id,COUNT(msi)::integer AS msi_bins
    FROM study_make30_final.mimic_msi_trajectory_long
    GROUP BY subject_id,hadm_id,stay_id
),
e AS
(
    SELECT subject_id,hadm_id,stay_id,
           COUNT(sodium)::integer AS sodium_bins,
           COUNT(potassium)::integer AS potassium_bins,
           COUNT(chloride)::integer AS chloride_bins,
           COUNT(bicarbonate)::integer AS bicarbonate_bins,
           COUNT(calcium)::integer AS calcium_bins,
           COUNT(magnesium)::integer AS magnesium_bins,
           COUNT(phosphorus)::integer AS phosphorus_bins
    FROM study_make30_final.mimic_electrolyte_trajectory_long
    GROUP BY subject_id,hadm_id,stay_id
)
SELECT
    m.subject_id,m.hadm_id,m.stay_id,
    m.msi_bins,
    e.sodium_bins,e.potassium_bins,e.chloride_bins,e.bicarbonate_bins,
    e.calcium_bins,e.magnesium_bins,e.phosphorus_bins,
    CASE WHEN m.msi_bins>=6 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END AS core4_eligible_msi6,
    CASE WHEN m.msi_bins>=9 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END AS core4_eligible_msi9,
    CASE WHEN m.msi_bins=12 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END AS core4_eligible_msi12,
    CASE WHEN m.msi_bins>=9 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END AS core4_trajectory_eligible,
    CASE WHEN m.msi_bins>=9 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 AND e.calcium_bins>=2 AND e.magnesium_bins>=2 AND e.phosphorus_bins>=2 THEN 1 ELSE 0 END AS all7_trajectory_eligible
FROM m INNER JOIN e USING(subject_id,hadm_id,stay_id);

CREATE UNIQUE INDEX idx_mimic_trajectory_coverage
ON study_make30_final.mimic_trajectory_coverage(stay_id);

/* ========================================================================
   19. PRIMARY CORE-4 PATIENT COHORT
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_core4_trajectory_cohort AS
SELECT c.*,cov.msi_bins,cov.sodium_bins,cov.potassium_bins,cov.chloride_bins,cov.bicarbonate_bins,
       cov.calcium_bins,cov.magnesium_bins,cov.phosphorus_bins,
       cov.core4_eligible_msi6,cov.core4_eligible_msi9,cov.core4_eligible_msi12,
       cov.core4_trajectory_eligible,cov.all7_trajectory_eligible
FROM study_make30_final.mimic_trajectory_pool c
INNER JOIN study_make30_final.mimic_trajectory_coverage cov ON c.stay_id=cov.stay_id
WHERE cov.core4_trajectory_eligible=1;

CREATE UNIQUE INDEX idx_mimic_core4_cohort_stay ON study_make30_final.mimic_core4_trajectory_cohort(stay_id);
CREATE UNIQUE INDEX idx_mimic_core4_cohort_subject ON study_make30_final.mimic_core4_trajectory_cohort(subject_id);

/* ========================================================================
   20. PRIMARY MIXED-FREQUENCY GBMTM INPUT TABLES
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_core4_msi_trajectory_long AS
SELECT t.* FROM study_make30_final.mimic_msi_trajectory_long t
INNER JOIN study_make30_final.mimic_trajectory_coverage c USING(stay_id)
WHERE c.core4_trajectory_eligible=1
ORDER BY t.stay_id,t.time_bin;
CREATE UNIQUE INDEX idx_mimic_core4_msi_long ON study_make30_final.mimic_core4_msi_trajectory_long(stay_id,time_bin);

CREATE TABLE study_make30_final.mimic_core4_electrolyte_trajectory_long AS
SELECT t.* FROM study_make30_final.mimic_electrolyte_trajectory_long t
INNER JOIN study_make30_final.mimic_trajectory_coverage c USING(stay_id)
WHERE c.core4_trajectory_eligible=1
ORDER BY t.stay_id,t.time_bin;
CREATE UNIQUE INDEX idx_mimic_core4_electrolyte_long ON study_make30_final.mimic_core4_electrolyte_trajectory_long(stay_id,time_bin);

CREATE TABLE study_make30_final.mimic_core4_multitrajectory_wide AS
WITH mw AS
(
    SELECT stay_id,
        MAX(msi) FILTER (WHERE time_bin=0) AS msi_m1, MAX(msi) FILTER (WHERE time_bin=1) AS msi_m2,
        MAX(msi) FILTER (WHERE time_bin=2) AS msi_m3, MAX(msi) FILTER (WHERE time_bin=3) AS msi_m4,
        MAX(msi) FILTER (WHERE time_bin=4) AS msi_m5, MAX(msi) FILTER (WHERE time_bin=5) AS msi_m6,
        MAX(msi) FILTER (WHERE time_bin=6) AS msi_m7, MAX(msi) FILTER (WHERE time_bin=7) AS msi_m8,
        MAX(msi) FILTER (WHERE time_bin=8) AS msi_m9, MAX(msi) FILTER (WHERE time_bin=9) AS msi_m10,
        MAX(msi) FILTER (WHERE time_bin=10) AS msi_m11, MAX(msi) FILTER (WHERE time_bin=11) AS msi_m12
    FROM study_make30_final.mimic_core4_msi_trajectory_long GROUP BY stay_id
),
ew AS
(
    SELECT stay_id,
        MAX(sodium) FILTER (WHERE time_bin=0) AS sodium_t1, MAX(sodium) FILTER (WHERE time_bin=1) AS sodium_t2, MAX(sodium) FILTER (WHERE time_bin=2) AS sodium_t3, MAX(sodium) FILTER (WHERE time_bin=3) AS sodium_t4,
        MAX(potassium) FILTER (WHERE time_bin=0) AS potassium_t1, MAX(potassium) FILTER (WHERE time_bin=1) AS potassium_t2, MAX(potassium) FILTER (WHERE time_bin=2) AS potassium_t3, MAX(potassium) FILTER (WHERE time_bin=3) AS potassium_t4,
        MAX(chloride) FILTER (WHERE time_bin=0) AS chloride_t1, MAX(chloride) FILTER (WHERE time_bin=1) AS chloride_t2, MAX(chloride) FILTER (WHERE time_bin=2) AS chloride_t3, MAX(chloride) FILTER (WHERE time_bin=3) AS chloride_t4,
        MAX(bicarbonate) FILTER (WHERE time_bin=0) AS bicarbonate_t1, MAX(bicarbonate) FILTER (WHERE time_bin=1) AS bicarbonate_t2, MAX(bicarbonate) FILTER (WHERE time_bin=2) AS bicarbonate_t3, MAX(bicarbonate) FILTER (WHERE time_bin=3) AS bicarbonate_t4
    FROM study_make30_final.mimic_core4_electrolyte_trajectory_long GROUP BY stay_id
)
SELECT c.subject_id,c.hadm_id,c.stay_id,c.sepsis_time,c.landmark_time,
       mw.msi_m1,mw.msi_m2,mw.msi_m3,mw.msi_m4,mw.msi_m5,mw.msi_m6,mw.msi_m7,mw.msi_m8,mw.msi_m9,mw.msi_m10,mw.msi_m11,mw.msi_m12,
       ew.sodium_t1,ew.sodium_t2,ew.sodium_t3,ew.sodium_t4,
       ew.potassium_t1,ew.potassium_t2,ew.potassium_t3,ew.potassium_t4,
       ew.chloride_t1,ew.chloride_t2,ew.chloride_t3,ew.chloride_t4,
       ew.bicarbonate_t1,ew.bicarbonate_t2,ew.bicarbonate_t3,ew.bicarbonate_t4
FROM study_make30_final.mimic_core4_trajectory_cohort c
LEFT JOIN mw USING(stay_id)
LEFT JOIN ew USING(stay_id);
CREATE UNIQUE INDEX idx_mimic_core4_multitrajectory_wide ON study_make30_final.mimic_core4_multitrajectory_wide(stay_id);

/* ========================================================================
   21. HOURLY MSI >4 SOURCE TRACE
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_msi_gt4_trace AS
SELECT c.subject_id,c.hadm_id,h.stay_id,h.hour_bin,
       FLOOR(h.hour_bin/2.0)::integer AS msi_2h_bin,
       ('M'||(FLOOR(h.hour_bin/2.0)::integer+1)::text)::text AS msi_2h_label,
       h.hr_median,
       h.map_invasive_median,
       h.map_noninvasive_median,
       h.map_abs_diff,
       h.map_used,h.map_source,h.msi AS hourly_msi,
       m.msi AS msi_2h_max
FROM study_make30_final.mimic_msi_hourly h
INNER JOIN study_make30_final.mimic_trajectory_pool c ON h.stay_id=c.stay_id
LEFT JOIN study_make30_final.mimic_msi_2h m ON h.stay_id=m.stay_id AND FLOOR(h.hour_bin/2.0)::integer=m.time_bin
WHERE h.msi>4
ORDER BY h.msi DESC,h.stay_id,h.hour_bin;

CREATE INDEX idx_mimic_msi_gt4_trace ON study_make30_final.mimic_msi_gt4_trace(stay_id,hour_bin);



/* ========================================================================
   22. ORAKLE MAKE30 + BASELINE / ADJUSTMENT COVARIATES

   FINAL STUDY CONVENTION
   ----------------------
   Trajectory exposure remains frozen at Sepsis t0 through t0+24 h.
   Outcome risk starts strictly AFTER the t0+24 h landmark.

   MAKE30 components follow ORAKLE:
     1) death;
     2) initiation of new dialysis/KRT;
     3) final inpatient serum creatinine >=200% of baseline creatinine.

   ORAKLE censoring is retained:
     earlier of hospital discharge or 30 days after INDEX ICU admission.

   Baseline creatinine follows ORAKLE:
     - median serum creatinine during the 12 months before hospital admission;
     - if unavailable, MDRD back-calculation assuming eGFR=75;
     - reference baseline = lower of that value and first admission creatinine.

   IMPORTANT:
     the frozen trajectory cohort itself is NOT changed here.
   ======================================================================== */

CREATE TABLE study_make30_final.mimic_baseline_kidney AS
WITH hist AS
(
    SELECT
        c.stay_id,
        PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY l.valuenum)::double precision AS historical_creatinine_median
    FROM study_make30_final.mimic_core4_trajectory_cohort c
    INNER JOIN mimiciv_hosp.labevents l
        ON l.subject_id = c.subject_id
       AND l.itemid = 50912
       AND l.charttime >= c.admittime - INTERVAL '365 days'
       AND l.charttime <  c.admittime
       AND l.valuenum > 0.1
       AND l.valuenum < 20
    GROUP BY c.stay_id
),
first_admit AS
(
    SELECT DISTINCT ON (c.stay_id)
        c.stay_id,
        l.charttime AS first_admission_creatinine_time,
        l.valuenum::double precision AS first_admission_creatinine
    FROM study_make30_final.mimic_core4_trajectory_cohort c
    INNER JOIN mimiciv_hosp.labevents l
        ON l.subject_id = c.subject_id
       AND l.itemid = 50912
       AND l.charttime >= c.admittime
       AND l.charttime <  c.landmark_time
       AND l.valuenum > 0.1
       AND l.valuenum < 20
    ORDER BY c.stay_id, l.charttime, l.labevent_id
),
calc AS
(
    SELECT
        c.subject_id,c.hadm_id,c.stay_id,
        c.age,c.sex,c.race,
        h.historical_creatinine_median,
        f.first_admission_creatinine_time,
        f.first_admission_creatinine,
        POWER(
            175.0
            * POWER(GREATEST(c.age::double precision,18.0),-0.203)
            * CASE WHEN UPPER(COALESCE(c.sex,''))='F' THEN 0.742 ELSE 1.0 END
            * CASE WHEN UPPER(COALESCE(c.race,'')) LIKE '%BLACK%' THEN 1.212 ELSE 1.0 END
            / 75.0,
            1.0/1.154
        )::double precision AS mdrd_backcalc_creatinine
    FROM study_make30_final.mimic_core4_trajectory_cohort c
    LEFT JOIN hist h USING(stay_id)
    LEFT JOIN first_admit f USING(stay_id)
),
base AS
(
    SELECT
        x.*,
        LEAST(
            COALESCE(x.historical_creatinine_median,x.mdrd_backcalc_creatinine),
            COALESCE(x.first_admission_creatinine,x.historical_creatinine_median,x.mdrd_backcalc_creatinine)
        )::double precision AS baseline_creatinine
    FROM calc x
)
SELECT
    b.*,
    CASE
        WHEN b.historical_creatinine_median IS NOT NULL
             AND b.baseline_creatinine = b.historical_creatinine_median
            THEN '12mo historical median'
        WHEN b.first_admission_creatinine IS NOT NULL
             AND b.baseline_creatinine = b.first_admission_creatinine
            THEN 'first admission creatinine'
        ELSE 'MDRD back-calculation eGFR=75'
    END::text AS baseline_creatinine_source,
    CASE WHEN b.baseline_creatinine IS NOT NULL AND b.baseline_creatinine>0 THEN
        175.0 * POWER(b.baseline_creatinine,-1.154)
        * POWER(GREATEST(b.age::double precision,18.0),-0.203)
        * CASE WHEN UPPER(COALESCE(b.sex,''))='F' THEN 0.742 ELSE 1.0 END
        * CASE WHEN UPPER(COALESCE(b.race,'')) LIKE '%BLACK%' THEN 1.212 ELSE 1.0 END
    END::double precision AS baseline_egfr_mdrd
FROM base b;

CREATE UNIQUE INDEX idx_mimic_baseline_kidney_stay
ON study_make30_final.mimic_baseline_kidney(stay_id);


CREATE TABLE study_make30_final.mimic_covariates AS
WITH wt AS
(
    SELECT DISTINCT ON (c.stay_id)
        c.stay_id,
        ce.valuenum::double precision AS weight_kg
    FROM study_make30_final.mimic_core4_trajectory_cohort c
    INNER JOIN mimiciv_icu.chartevents ce
        ON ce.stay_id=c.stay_id
       AND ce.itemid IN (226512,224639)
       AND ce.charttime < c.landmark_time
       AND ce.valuenum BETWEEN 25 AND 300
    ORDER BY c.stay_id,
             CASE WHEN ce.itemid=226512 THEN 0 ELSE 1 END,
             ce.charttime
),
ht0 AS
(
    SELECT
        c.stay_id,
        ce.charttime,
        ce.itemid,
        CASE
            WHEN ce.itemid=226730 AND ce.valuenum BETWEEN 100 AND 250 THEN ce.valuenum
            WHEN ce.itemid=226707 AND ce.valuenum BETWEEN 39.37 AND 98.43 THEN ce.valuenum*2.54
        END::double precision AS height_cm
    FROM study_make30_final.mimic_core4_trajectory_cohort c
    INNER JOIN mimiciv_icu.chartevents ce
        ON ce.stay_id=c.stay_id
       AND ce.itemid IN (226730,226707)
       AND ce.charttime < c.landmark_time
       AND ce.valuenum IS NOT NULL
),
ht AS
(
    SELECT DISTINCT ON (stay_id)
        stay_id,height_cm
    FROM ht0
    WHERE height_cm BETWEEN 100 AND 250
    ORDER BY stay_id,CASE WHEN itemid=226730 THEN 0 ELSE 1 END,charttime
),
vent AS
(
    SELECT c.stay_id,1::integer AS mechanical_ventilation_0_24h
    FROM study_make30_final.mimic_core4_trajectory_cohort c
    WHERE EXISTS
    (
        SELECT 1 FROM mimiciv_derived.ventilation v
        WHERE v.stay_id=c.stay_id
          AND v.starttime < c.landmark_time
          AND v.endtime > c.sepsis_time
          AND v.ventilation_status IN ('InvasiveVent','Tracheostomy')
    )
),
vaso AS
(
    /* Raw MIMIC-IV inputevents implementation avoids dependence on an
       optional locally materialized derived vasoactive_agent table.
       ItemIDs: norepinephrine 221906; epinephrine 221289;
       phenylephrine 221749; vasopressin 222315; dopamine 221662. */
    SELECT c.stay_id,1::integer AS vasopressor_0_24h
    FROM study_make30_final.mimic_core4_trajectory_cohort c
    WHERE EXISTS
    (
        SELECT 1
        FROM mimiciv_icu.inputevents ie
        WHERE ie.stay_id=c.stay_id
          AND ie.itemid IN (221906,221289,221749,222315,221662)
          AND ie.starttime < c.landmark_time
          AND COALESCE(ie.endtime,ie.starttime) > c.sepsis_time
          AND (COALESCE(ie.rate,0)>0 OR COALESCE(ie.amount,0)>0)
          AND COALESCE(ie.statusdescription,'') <> 'Rewritten'
    )
)
SELECT
    c.subject_id,c.hadm_id,c.stay_id,
    c.age::double precision AS age,
    c.sex,
    c.race AS race_ethnicity,
    w.weight_kg,
    h.height_cm,
    CASE WHEN w.weight_kg IS NOT NULL AND h.height_cm IS NOT NULL AND h.height_cm>0
         THEN w.weight_kg/POWER(h.height_cm/100.0,2) END::double precision AS bmi,
    c.sofa_score::double precision AS sepsis_sofa_score,
    bk.historical_creatinine_median,
    bk.first_admission_creatinine,
    bk.mdrd_backcalc_creatinine,
    bk.baseline_creatinine_source,
    bk.baseline_creatinine,
    bk.baseline_egfr_mdrd,
    COALESCE(ch.renal_disease,0)::integer AS ckd_or_renal_disease,
    COALESCE(ch.diabetes_without_cc,0)::integer AS diabetes_without_complication,
    COALESCE(ch.diabetes_with_cc,0)::integer AS diabetes_with_complication,
    GREATEST(COALESCE(ch.diabetes_without_cc,0),COALESCE(ch.diabetes_with_cc,0))::integer AS diabetes,
    COALESCE(ch.congestive_heart_failure,0)::integer AS congestive_heart_failure,
    GREATEST(COALESCE(ch.mild_liver_disease,0),COALESCE(ch.severe_liver_disease,0))::integer AS liver_disease,
    GREATEST(COALESCE(ch.malignant_cancer,0),COALESCE(ch.metastatic_solid_tumor,0))::integer AS malignancy,
    ch.charlson_comorbidity_index::double precision AS charlson_comorbidity_index,
    COALESCE(v.mechanical_ventilation_0_24h,0)::integer AS mechanical_ventilation_0_24h,
    COALESCE(p.vasopressor_0_24h,0)::integer AS vasopressor_0_24h
FROM study_make30_final.mimic_core4_trajectory_cohort c
LEFT JOIN wt w USING(stay_id)
LEFT JOIN ht h USING(stay_id)
LEFT JOIN study_make30_final.mimic_baseline_kidney bk USING(stay_id)
LEFT JOIN mimiciv_derived.charlson ch ON ch.hadm_id=c.hadm_id
LEFT JOIN vent v USING(stay_id)
LEFT JOIN vaso p USING(stay_id);

CREATE UNIQUE INDEX idx_mimic_covariates_stay
ON study_make30_final.mimic_covariates(stay_id);


CREATE TABLE study_make30_final.mimic_make30 AS
WITH fup AS
(
    SELECT
        c.subject_id,c.hadm_id,c.stay_id,c.sepsis_time,c.landmark_time,c.icu_intime,c.dischtime,
        LEAST(
            COALESCE(c.dischtime,c.icu_intime+INTERVAL '30 days'),
            c.icu_intime+INTERVAL '30 days'
        ) AS make30_censor_time
    FROM study_make30_final.mimic_core4_trajectory_cohort c
),
death AS
(
    SELECT
        f.stay_id,
        CASE WHEN a.hospital_expire_flag=1
                   AND COALESCE(a.deathtime,a.dischtime) > f.landmark_time
                   AND COALESCE(a.deathtime,a.dischtime) <= f.make30_censor_time
             THEN 1 ELSE 0 END::integer AS make30_death,
        CASE WHEN a.hospital_expire_flag=1 THEN COALESCE(a.deathtime,a.dischtime) END AS death_time
    FROM fup f
    INNER JOIN mimiciv_hosp.admissions a ON a.hadm_id=f.hadm_id
),
prior_krt AS
(
    SELECT f.stay_id,1::integer AS prior_krt_before_landmark
    FROM fup f
    WHERE EXISTS
    (
        SELECT 1
        FROM mimiciv_icu.icustays ie
        INNER JOIN mimiciv_derived.rrt r ON r.stay_id=ie.stay_id
        WHERE ie.hadm_id=f.hadm_id
          AND r.charttime <= f.landmark_time
          AND (r.dialysis_active=1 OR NULLIF(BTRIM(r.dialysis_type),'') IS NOT NULL)
    )
),
krt AS
(
    SELECT
        f.stay_id,
        MIN(r.charttime) AS first_new_krt_time
    FROM fup f
    LEFT JOIN prior_krt pk USING(stay_id)
    INNER JOIN mimiciv_icu.icustays ie ON ie.hadm_id=f.hadm_id
    INNER JOIN mimiciv_derived.rrt r
        ON r.stay_id=ie.stay_id
       AND r.charttime > f.landmark_time
       AND r.charttime <= f.make30_censor_time
       AND (r.dialysis_active=1 OR NULLIF(BTRIM(r.dialysis_type),'') IS NOT NULL)
    WHERE pk.stay_id IS NULL
    GROUP BY f.stay_id
),
final_cr AS
(
    SELECT DISTINCT ON (f.stay_id)
        f.stay_id,
        l.charttime AS final_inpatient_creatinine_time,
        l.valuenum::double precision AS final_inpatient_creatinine
    FROM fup f
    INNER JOIN mimiciv_hosp.labevents l
        ON l.subject_id=f.subject_id
       AND l.itemid=50912
       AND l.charttime > f.landmark_time
       AND l.charttime <= f.make30_censor_time
       AND l.valuenum > 0.1
       AND l.valuenum < 20
    ORDER BY f.stay_id,l.charttime DESC,l.labevent_id DESC
)
SELECT
    f.subject_id,f.hadm_id,f.stay_id,
    f.landmark_time AS make30_risk_start,
    f.make30_censor_time,
    EXTRACT(EPOCH FROM (f.make30_censor_time-f.landmark_time))/86400.0 AS make30_followup_days_after_landmark,
    CASE WHEN f.make30_censor_time>f.landmark_time THEN 1 ELSE 0 END::integer AS make30_evaluable_window,
    d.make30_death,
    d.death_time,
    COALESCE(pk.prior_krt_before_landmark,0)::integer AS prior_krt_before_landmark,
    CASE WHEN k.first_new_krt_time IS NOT NULL THEN 1 ELSE 0 END::integer AS make30_new_krt,
    k.first_new_krt_time,
    fc.final_inpatient_creatinine,
    fc.final_inpatient_creatinine_time,
    bk.baseline_creatinine,
    CASE WHEN fc.final_inpatient_creatinine IS NOT NULL AND bk.baseline_creatinine IS NOT NULL
                   AND fc.final_inpatient_creatinine >= 2.0*bk.baseline_creatinine
         THEN 1
         WHEN fc.final_inpatient_creatinine IS NOT NULL AND bk.baseline_creatinine IS NOT NULL
         THEN 0
         ELSE NULL END::integer AS make30_persistent_kidney_dysfunction,
    CASE
        WHEN f.make30_censor_time<=f.landmark_time THEN NULL
        WHEN COALESCE(d.make30_death,0)=1 OR k.first_new_krt_time IS NOT NULL THEN 1
        WHEN fc.final_inpatient_creatinine IS NOT NULL AND bk.baseline_creatinine IS NOT NULL
             AND fc.final_inpatient_creatinine >= 2.0*bk.baseline_creatinine THEN 1
        WHEN fc.final_inpatient_creatinine IS NOT NULL AND bk.baseline_creatinine IS NOT NULL THEN 0
        ELSE NULL
    END::integer AS make30
FROM fup f
LEFT JOIN death d USING(stay_id)
LEFT JOIN prior_krt pk USING(stay_id)
LEFT JOIN krt k USING(stay_id)
LEFT JOIN final_cr fc USING(stay_id)
LEFT JOIN study_make30_final.mimic_baseline_kidney bk USING(stay_id);

CREATE UNIQUE INDEX idx_mimic_make30_stay
ON study_make30_final.mimic_make30(stay_id);


CREATE TABLE study_make30_final.mimic_analysis_master AS
SELECT
    w.*,
    cv.age,cv.sex,cv.race_ethnicity,cv.weight_kg,cv.height_cm,cv.bmi,
    cv.sepsis_sofa_score,cv.historical_creatinine_median,cv.first_admission_creatinine,cv.mdrd_backcalc_creatinine,cv.baseline_creatinine_source,
    cv.baseline_creatinine,cv.baseline_egfr_mdrd,
    cv.ckd_or_renal_disease,cv.diabetes_without_complication,cv.diabetes_with_complication,cv.diabetes,
    cv.congestive_heart_failure,cv.liver_disease,cv.malignancy,cv.charlson_comorbidity_index,
    cv.mechanical_ventilation_0_24h,cv.vasopressor_0_24h,
    o.make30_risk_start,o.make30_censor_time,o.make30_followup_days_after_landmark,o.make30_evaluable_window,
    o.make30_death,o.death_time,o.prior_krt_before_landmark,o.make30_new_krt,o.first_new_krt_time,
    o.final_inpatient_creatinine,o.final_inpatient_creatinine_time,o.make30_persistent_kidney_dysfunction,o.make30
FROM study_make30_final.mimic_core4_multitrajectory_wide w
LEFT JOIN study_make30_final.mimic_covariates cv USING(stay_id)
LEFT JOIN study_make30_final.mimic_make30 o USING(stay_id);

CREATE UNIQUE INDEX idx_mimic_analysis_master_stay
ON study_make30_final.mimic_analysis_master(stay_id);


/* ========================================================================
   23. FINAL QC + ANALYSIS RESULT
   ======================================================================== */
CREATE TABLE study_make30_final.mimic_qc_summary AS
WITH
implementation_info AS
(
    SELECT '00_IMPLEMENTATION_INFO'::text AS section,0::integer AS sort1,1::integer AS sort2,
           'FINAL analysis pipeline: trajectories + landmark-start MAKE30 + covariates; MSI: 12 x 2-h windows; Na/K/Cl/HCO3: 4 x 6-h windows; hourly MAP QC = both invasive MAP and NIBP <=20 excluded before hourly summaries; if both valid sources exist and |invasive-NIBP| >=30 use NIBP; otherwise valid invasive preferred; MSI_2h = maximum hourly MSI; primary completeness MSI>=9/12 and each electrolyte>=2/4'::text AS item,
           NULL::text AS time_bin,NULL::bigint AS n_total,NULL::bigint AS n_observed,NULL::bigint AS n_missing,NULL::numeric AS observed_pct,NULL::numeric AS missing_pct,
           NULL::bigint AS n_meet_requirement,NULL::bigint AS n_all_bins,
           NULL::numeric AS mean_value,NULL::numeric AS sd_value,NULL::numeric AS min_value,NULL::numeric AS p01,NULL::numeric AS p05,NULL::numeric AS p25,NULL::numeric AS median_value,NULL::numeric AS p75,NULL::numeric AS p95,NULL::numeric AS p99,NULL::numeric AS max_value,NULL::text AS detail
),
cohort_flow AS
(
    SELECT '01_COHORT_FLOW'::text section,1 sort1,1 sort2,'all ICU patients'::text item,NULL::text time_bin,COUNT(*)::bigint n_total,COUNT(*)::bigint n_observed,NULL::bigint n_missing,100::numeric observed_pct,0::numeric missing_pct,NULL::bigint n_meet_requirement,NULL::bigint n_all_bins,NULL::numeric mean_value,NULL::numeric sd_value,NULL::numeric min_value,NULL::numeric p01,NULL::numeric p05,NULL::numeric p25,NULL::numeric median_value,NULL::numeric p75,NULL::numeric p95,NULL::numeric p99,NULL::numeric max_value,NULL::text detail FROM study_make30_final.mimic_patients
    UNION ALL SELECT '01_COHORT_FLOW',1,2,'adult patients with ORAKLE-style Sepsis onset in ICU',NULL,(SELECT COUNT(*) FROM study_make30_final.mimic_patients),COUNT(DISTINCT subject_id),NULL,ROUND(100.0*COUNT(DISTINCT subject_id)/NULLIF((SELECT COUNT(*) FROM study_make30_final.mimic_patients),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.mimic_sepsis_adult
    UNION ALL SELECT '01_COHORT_FLOW',1,3,'first qualifying Sepsis ICU stay per patient',NULL,(SELECT COUNT(*) FROM study_make30_final.mimic_patients),COUNT(*),NULL,ROUND(100.0*COUNT(*)/NULLIF((SELECT COUNT(*) FROM study_make30_final.mimic_patients),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.mimic_sepsis_index
    UNION ALL SELECT '01_COHORT_FLOW',1,4,'reached Sepsis +24 h landmark in ICU and hospital',NULL,COUNT(*),COUNT(*) FILTER(WHERE reached_24h_landmark=1),COUNT(*) FILTER(WHERE reached_24h_landmark=0),ROUND(100.0*COUNT(*) FILTER(WHERE reached_24h_landmark=1)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE reached_24h_landmark=0)/NULLIF(COUNT(*),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.mimic_landmark_cohort
    UNION ALL SELECT '01_COHORT_FLOW',1,5,'trajectory pool after 24 h landmark + KRT exclusion',NULL,(SELECT COUNT(*) FROM study_make30_final.mimic_landmark_cohort WHERE reached_24h_landmark=1),COUNT(*),NULL,ROUND(100.0*COUNT(*)/NULLIF((SELECT COUNT(*) FROM study_make30_final.mimic_landmark_cohort WHERE reached_24h_landmark=1),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.mimic_trajectory_pool
    UNION ALL SELECT '01_COHORT_FLOW',1,6,'PRIMARY mixed-frequency Core-4 eligible: MSI>=9/12 + Na/K/Cl/HCO3 each>=2/4',NULL,COUNT(*),COUNT(*) FILTER(WHERE core4_trajectory_eligible=1),COUNT(*) FILTER(WHERE core4_trajectory_eligible=0),ROUND(100.0*COUNT(*) FILTER(WHERE core4_trajectory_eligible=1)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE core4_trajectory_eligible=0)/NULLIF(COUNT(*),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.mimic_trajectory_coverage
),
variable_list AS
(
    SELECT * FROM (VALUES (1,'MSI'::text,9,12),(2,'sodium'::text,2,4),(3,'potassium'::text,2,4),(4,'chloride'::text,2,4),(5,'bicarbonate'::text,2,4)) v(variable_order,variable_name,required_bins,max_bins)
),
patient_coverage AS
(
    SELECT '02_PATIENT_LEVEL_COVERAGE'::text section,2 sort1,v.variable_order sort2,v.variable_name item,NULL::text time_bin,COUNT(*)::bigint n_total,
           COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END>=1)::bigint n_observed,
           COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END=0)::bigint n_missing,
           ROUND(100.0*COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END>=1)/NULLIF(COUNT(*),0),2) observed_pct,
           ROUND(100.0*COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END=0)/NULLIF(COUNT(*),0),2) missing_pct,
           COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END>=v.required_bins)::bigint n_meet_requirement,
           COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END=v.max_bins)::bigint n_all_bins,
           NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,
           ('required_bins='||v.required_bins||'; max_bins='||v.max_bins)::text detail
    FROM study_make30_final.mimic_trajectory_coverage c CROSS JOIN variable_list v GROUP BY v.variable_order,v.variable_name,v.required_bins,v.max_bins
),
all_series AS
(
    SELECT stay_id,1 variable_order,'MSI'::text variable_name,time_bin,time_label,msi::double precision value_num FROM study_make30_final.mimic_msi_trajectory_long
    UNION ALL SELECT stay_id,2,'sodium',time_bin,time_label,sodium::double precision FROM study_make30_final.mimic_electrolyte_trajectory_long
    UNION ALL SELECT stay_id,3,'potassium',time_bin,time_label,potassium::double precision FROM study_make30_final.mimic_electrolyte_trajectory_long
    UNION ALL SELECT stay_id,4,'chloride',time_bin,time_label,chloride::double precision FROM study_make30_final.mimic_electrolyte_trajectory_long
    UNION ALL SELECT stay_id,5,'bicarbonate',time_bin,time_label,bicarbonate::double precision FROM study_make30_final.mimic_electrolyte_trajectory_long
),
eligible_series AS
(
    SELECT stay_id,1 variable_order,'MSI'::text variable_name,time_bin,time_label,msi::double precision value_num FROM study_make30_final.mimic_core4_msi_trajectory_long
    UNION ALL SELECT stay_id,2,'sodium',time_bin,time_label,sodium::double precision FROM study_make30_final.mimic_core4_electrolyte_trajectory_long
    UNION ALL SELECT stay_id,3,'potassium',time_bin,time_label,potassium::double precision FROM study_make30_final.mimic_core4_electrolyte_trajectory_long
    UNION ALL SELECT stay_id,4,'chloride',time_bin,time_label,chloride::double precision FROM study_make30_final.mimic_core4_electrolyte_trajectory_long
    UNION ALL SELECT stay_id,5,'bicarbonate',time_bin,time_label,bicarbonate::double precision FROM study_make30_final.mimic_core4_electrolyte_trajectory_long
),
all_missing AS
(
    SELECT '03_BIN_LEVEL_MISSINGNESS'::text section,3 sort1,(variable_order*100+time_bin)::integer sort2,variable_name item,time_label time_bin,COUNT(*)::bigint n_total,COUNT(value_num)::bigint n_observed,COUNT(*) FILTER(WHERE value_num IS NULL)::bigint n_missing,ROUND(100.0*COUNT(value_num)/NULLIF(COUNT(*),0),2) observed_pct,ROUND(100.0*COUNT(*) FILTER(WHERE value_num IS NULL)/NULLIF(COUNT(*),0),2) missing_pct,NULL::bigint n_meet_requirement,NULL::bigint n_all_bins,NULL::numeric mean_value,NULL::numeric sd_value,NULL::numeric min_value,NULL::numeric p01,NULL::numeric p05,NULL::numeric p25,NULL::numeric median_value,NULL::numeric p75,NULL::numeric p95,NULL::numeric p99,NULL::numeric max_value,NULL::text detail FROM all_series GROUP BY variable_order,variable_name,time_bin,time_label
),
elig_missing AS
(
    SELECT '04_CORE4_ELIGIBLE_BIN_MISSINGNESS'::text section,4 sort1,(variable_order*100+time_bin)::integer sort2,variable_name item,time_label time_bin,COUNT(*)::bigint n_total,COUNT(value_num)::bigint n_observed,COUNT(*) FILTER(WHERE value_num IS NULL)::bigint n_missing,ROUND(100.0*COUNT(value_num)/NULLIF(COUNT(*),0),2) observed_pct,ROUND(100.0*COUNT(*) FILTER(WHERE value_num IS NULL)/NULLIF(COUNT(*),0),2) missing_pct,NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::text FROM eligible_series GROUP BY variable_order,variable_name,time_bin,time_label
),
value_qc AS
(
    SELECT '05_CORE4_VALUE_DISTRIBUTION_QC'::text section,5 sort1,(variable_order*100+time_bin)::integer sort2,variable_name item,time_label time_bin,COUNT(*)::bigint n_total,COUNT(value_num)::bigint n_observed,COUNT(*) FILTER(WHERE value_num IS NULL)::bigint n_missing,ROUND(100.0*COUNT(value_num)/NULLIF(COUNT(*),0),2) observed_pct,ROUND(100.0*COUNT(*) FILTER(WHERE value_num IS NULL)/NULLIF(COUNT(*),0),2) missing_pct,NULL::bigint,NULL::bigint,AVG(value_num)::numeric,STDDEV_SAMP(value_num)::numeric,MIN(value_num)::numeric,PERCENTILE_CONT(.01) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.05) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.25) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.5) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.75) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.95) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.99) WITHIN GROUP(ORDER BY value_num)::numeric,MAX(value_num)::numeric,NULL::text FROM eligible_series GROUP BY variable_order,variable_name,time_bin,time_label
),
threshold_qc AS
(
    SELECT '06_MSI_COVERAGE_THRESHOLD_SENSITIVITY'::text section,6 sort1,x.ord sort2,x.label item,NULL::text time_bin,COUNT(*)::bigint n_total,COUNT(*) FILTER(WHERE x.flag=1)::bigint n_observed,COUNT(*) FILTER(WHERE x.flag=0)::bigint n_missing,ROUND(100.0*COUNT(*) FILTER(WHERE x.flag=1)/NULLIF(COUNT(*),0),2) observed_pct,ROUND(100.0*COUNT(*) FILTER(WHERE x.flag=0)/NULLIF(COUNT(*),0),2) missing_pct,NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::text
    FROM study_make30_final.mimic_trajectory_coverage c CROSS JOIN LATERAL (VALUES (1,'Core-4 if MSI>=6/12',c.core4_eligible_msi6),(2,'Core-4 if MSI>=9/12 [PRIMARY]',c.core4_eligible_msi9),(3,'Core-4 if MSI=12/12',c.core4_eligible_msi12)) x(ord,label,flag) GROUP BY x.ord,x.label
),
anomaly_trace AS
(
    SELECT '07_MSI_GT4_TRACE'::text section,7 sort1,ROW_NUMBER() OVER(ORDER BY hourly_msi DESC,stay_id,hour_bin)::integer sort2,'hourly MSI > 4'::text item,msi_2h_label::text time_bin,1::bigint n_total,1::bigint n_observed,0::bigint n_missing,100::numeric observed_pct,0::numeric missing_pct,NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,
           jsonb_build_object('subject_id',subject_id,'hadm_id',hadm_id,'stay_id',stay_id,'hour_bin',hour_bin,'msi_2h_bin',msi_2h_bin,'hr_median',hr_median,'map_invasive_median',map_invasive_median,'map_noninvasive_median',map_noninvasive_median,'map_abs_diff',map_abs_diff,'map_used',map_used,'map_source',map_source,'hourly_msi',hourly_msi,'msi_2h_max',msi_2h_max)::text detail
    FROM study_make30_final.mimic_msi_gt4_trace
),
structure_qc AS
(
    SELECT '08_STRUCTURE_QC'::text section,8 sort1,1 sort2,'eligible patients'::text item,NULL::text time_bin,COUNT(*)::bigint n_total,NULL::bigint n_observed,NULL::bigint n_missing,NULL::numeric observed_pct,NULL::numeric missing_pct,NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::text FROM study_make30_final.mimic_core4_trajectory_cohort
    UNION ALL SELECT '08_STRUCTURE_QC',8,2,'MSI long rows (expected eligible N x 12)',NULL,COUNT(*),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.mimic_core4_msi_trajectory_long
    UNION ALL SELECT '08_STRUCTURE_QC',8,3,'electrolyte long rows (expected eligible N x 4)',NULL,COUNT(*),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.mimic_core4_electrolyte_trajectory_long
)
SELECT section,item,time_bin,n_total,n_observed,n_missing,observed_pct,missing_pct,n_meet_requirement,n_all_bins,
       ROUND(mean_value,4) mean_value,ROUND(sd_value,4) sd_value,ROUND(min_value,4) min_value,ROUND(p01,4) p01,ROUND(p05,4) p05,ROUND(p25,4) p25,ROUND(median_value,4) median_value,ROUND(p75,4) p75,ROUND(p95,4) p95,ROUND(p99,4) p99,ROUND(max_value,4) max_value,detail
FROM (
    SELECT * FROM implementation_info UNION ALL SELECT * FROM cohort_flow UNION ALL SELECT * FROM patient_coverage UNION ALL SELECT * FROM all_missing UNION ALL SELECT * FROM elig_missing UNION ALL SELECT * FROM value_qc UNION ALL SELECT * FROM threshold_qc UNION ALL SELECT * FROM anomaly_trace UNION ALL SELECT * FROM structure_qc
) z
ORDER BY sort1,sort2;


INSERT INTO study_make30_final.mimic_qc_summary
(section,item,time_bin,n_total,n_observed,n_missing,observed_pct,missing_pct,n_meet_requirement,n_all_bins,mean_value,sd_value,min_value,p01,p05,p25,median_value,p75,p95,p99,max_value,detail)
SELECT '09_MAKE30_QC','MAKE30 evaluable / events',NULL,COUNT(*)::bigint,COUNT(make30)::bigint,COUNT(*) FILTER(WHERE make30 IS NULL)::bigint,
       ROUND(100.0*COUNT(make30)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE make30 IS NULL)/NULLIF(COUNT(*),0),2),
       COUNT(*) FILTER(WHERE make30=1)::bigint,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
       jsonb_build_object('death',COUNT(*) FILTER(WHERE make30_death=1),'new_krt',COUNT(*) FILTER(WHERE make30_new_krt=1),'persistent_kidney_dysfunction',COUNT(*) FILTER(WHERE make30_persistent_kidney_dysfunction=1),'missing_final_creatinine',COUNT(*) FILTER(WHERE final_inpatient_creatinine IS NULL))::text
FROM study_make30_final.mimic_make30;

INSERT INTO study_make30_final.mimic_qc_summary
(section,item,time_bin,n_total,n_observed,n_missing,observed_pct,missing_pct,n_meet_requirement,n_all_bins,mean_value,sd_value,min_value,p01,p05,p25,median_value,p75,p95,p99,max_value,detail)
SELECT '10_COVARIATE_QC',v.var,NULL,COUNT(*)::bigint,COUNT(v.val)::bigint,(COUNT(*)-COUNT(v.val))::bigint,
       ROUND(100.0*COUNT(v.val)/NULLIF(COUNT(*),0),2),ROUND(100.0*(COUNT(*)-COUNT(v.val))/NULLIF(COUNT(*),0),2),NULL,NULL,
       AVG(v.val),STDDEV_SAMP(v.val),MIN(v.val),PERCENTILE_CONT(.01) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.05) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.25) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.5) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.75) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.95) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.99) WITHIN GROUP(ORDER BY v.val),MAX(v.val),NULL
FROM study_make30_final.mimic_covariates c
CROSS JOIN LATERAL (VALUES ('age',c.age),('bmi',c.bmi),('sepsis_sofa_score',c.sepsis_sofa_score),('baseline_creatinine',c.baseline_creatinine),('baseline_egfr_mdrd',c.baseline_egfr_mdrd),('charlson_comorbidity_index',c.charlson_comorbidity_index)) v(var,val)
GROUP BY v.var;

/* FINAL USER-VISIBLE RESULT: one row per eligible patient, including all
   frozen trajectory values + MAKE30 + covariates. QC remains persisted in
   study_make30_final.mimic_qc_summary. */
SELECT * FROM study_make30_final.mimic_analysis_master ORDER BY subject_id,stay_id;
