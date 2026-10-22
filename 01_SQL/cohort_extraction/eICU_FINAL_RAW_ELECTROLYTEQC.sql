/* FINAL FULL-RAW ANALYSIS PIPELINE: reconstructs the complete eICU cohort,
   ORAKLE-style Sepsis/SOFA, frozen MSI/electrolyte trajectories, ORAKLE MAKE30,
   and adjustment covariates from raw eICU tables in one run.
   Outcome events are counted only after Sepsis t0+24h landmark.
   Performance optimization changes execution strategy only, not definitions.
   Generated 2026-08-27. */

/* ========================================================================
   MSI STUDY — eICU-CRD v2.0
   FULL-COHORT FORMAL RESEARCH PIPELINE — MIXED-FREQUENCY CORE-4 VERSION: MSI 2-H + ELECTROLYTES 6-H
   PostgreSQL / directly runnable

   PRIMARY TRAJECTORY:
       MSI + Na + K + Cl + HCO3

   PURPOSE
   -------
   Run the formal study workflow on the full eICU patient cohort:

     1) full cohort at patient (uniquepid) level
     2) retain all ICU stays for all patients
     3) identify ORAKLE-style suspected infection in eICU
     4) pair suspected infection with hourly SOFA >=2
     5) define Sepsis onset (t0)
     6) retain first qualifying Sepsis ICU stay per patient
     7) require observation in ICU/hospital through t0 + 24 h
     8) exclude KRT during Sepsis 0–24 h
     9) derive hourly MSI and twelve 2-h MSI windows
    10) derive Na/K/Cl/HCO3 in four 6-h windows
    11) create patient-level coverage
    12) create the PRIMARY Core-4 trajectory cohort
    13) output cohort flow + coverage + bin-level missingness

   ------------------------------------------------------------------------
   TIME SCALE
   ------------------------------------------------------------------------
   eICU offsets are in minutes relative to ICU admission.

   Exposure window:
       [sepsis_offset, sepsis_offset + 1440)

   Four 6-h bins:
       T1 0–6 h
       T2 6–12 h
       T3 12–18 h
       T4 18–24 h

   ------------------------------------------------------------------------
   ORAKLE-STYLE SEPSIS
   ------------------------------------------------------------------------
   Published ORAKLE rule:
     - suspected infection + SOFA increase >=2
     - baseline SOFA assumed 0 before ICU
     - qualifying SOFA time no more than 24 h before
       or 12 h after suspected infection
     - sepsis onset = earlier of suspected-infection time and SOFA time
     - in eICU, antibiotic administrations are used as the proxy for
       suspected infection because fluid-sampling data are limited

   In eICU, this script follows the published ORAKLE wording directly:
   the OCCURRENCE OF AN ANTIBIOTIC ADMINISTRATION is used to define the
   suspected-infection anchor. No requirement for >=2 administrations is
   imposed.

   Each unique recognized antibiotic start time is therefore a candidate
   suspected-infection time. It must still satisfy the ORAKLE SOFA pairing:
      SOFA >=2 from 24 h before through 12 h after the antibiotic anchor.

   IMPORTANT:
   ORAKLE source code is not public. Therefore the exact author-side raw
   drug-string list cannot be guaranteed. This script keeps the transparent
   explicit 25-pattern antibiotic dictionary so the mapping is auditable.

   ------------------------------------------------------------------------
   SOFA SOURCE
   ------------------------------------------------------------------------
   This script DOES NOT use the abnormal local eicu_derived.sepsis3 table
   and DOES NOT require a pre-existing hourly SOFA table.

   Instead, the standard original SOFA score is reconstructed directly from
   raw eICU tables for the all patients:
      nursecharting       -> GCS
      lab                 -> PaO2/FiO2, platelets, bilirubin, creatinine
      respiratorycare     -> invasive ventilatory support
      vitalperiodic       -> invasive MAP
      vitalaperiodic      -> non-invasive MAP
      treatment           -> dose-coded vasopressor / inotrope categories
      medication          -> dobutamine exposure if treatment coding absent
      intakeoutput        -> urine output

   Each component is scored using the original 1996 SOFA thresholds.
   The worst component score over the current + preceding 23 hourly rows is
   used to form a rolling 24-h hourly SOFA, paralleling the hourly rolling
   approach used in MIMIC-derived SOFA.

   Missing component values contribute 0 points, as in common database SOFA
   implementations after the component-specific score is unavailable.

   ------------------------------------------------------------------------
   MSI
   ------------------------------------------------------------------------
   Hourly:
      median HR from public.vitalperiodic.heartrate
      median invasive MAP from public.vitalperiodic.systemicmean
      median non-invasive MAP from public.vitalaperiodic.noninvasivemean
      MAP selection after hourly summaries:
        both invasive MAP and NIBP <=20 mmHg -> exclude before hourly summaries
        if both valid hourly sources exist and |invasive-NIBP| >=30 mmHg -> use NIBP
        otherwise valid invasive MAP has priority

      MSI_hour = median HR / prioritized median MAP

   2-h:
      MSI_2h = maximum valid hourly MSI in each 2-h window

   MSI-input plausibility:
      HR  >0 and <300
      invasive MAP and NIBP >20 and <300
      (SOFA reconstruction retains its own original MAP handling)

   ------------------------------------------------------------------------
   ELECTROLYTES
   ------------------------------------------------------------------------
      Na    = labname 'sodium'
      K     = labname 'potassium'
      Cl    = labname 'chloride'
      HCO3  = labname 'bicarbonate'

   Electrolyte hard-QC before 6-h aggregation:
      sodium >300 mmol/L is discarded (300 retained).
      chloride <10 mmol/L is discarded (10 retained).
      Filtering is applied to raw lab rows before the 6-h median.

   Value in each 6-h bin:
      median of all available measurements

   Missing bins remain NULL.
   No LOCF / NOCB / interpolation.

   ------------------------------------------------------------------------
   PRIMARY CORE-4 ELIGIBILITY
   ------------------------------------------------------------------------
      MSI >=9 distinct 2-h bins
      Na >=2
      K >=2
      Cl >=2
      HCO3 >=2

   ======================================================================== */


/* ========================================================================
   0. SCHEMA + CLEANUP
   ======================================================================== */

CREATE SCHEMA IF NOT EXISTS study_make30_final;

DROP TABLE IF EXISTS study_make30_final.eicu_analysis_master;
DROP TABLE IF EXISTS study_make30_final.eicu_make30;
DROP TABLE IF EXISTS study_make30_final.eicu_covariates;
DROP TABLE IF EXISTS study_make30_final.eicu_baseline_kidney;
DROP TABLE IF EXISTS study_make30_final.eicu_qc_summary;


DROP TABLE IF EXISTS study_make30_final.eicu_msi_gt4_trace;
DROP TABLE IF EXISTS study_make30_final.eicu_core4_multitrajectory_wide;
DROP TABLE IF EXISTS study_make30_final.eicu_core4_electrolyte_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.eicu_core4_msi_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.eicu_core4_trajectory_cohort;
DROP TABLE IF EXISTS study_make30_final.eicu_trajectory_coverage;
DROP TABLE IF EXISTS study_make30_final.eicu_electrolyte_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.eicu_msi_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.eicu_electrolyte_6h;
DROP TABLE IF EXISTS study_make30_final.eicu_msi_2h;
DROP TABLE IF EXISTS study_make30_final.eicu_msi_hourly;
DROP TABLE IF EXISTS study_make30_final.eicu_trajectory_pool;
DROP TABLE IF EXISTS study_make30_final.eicu_landmark_cohort;
DROP TABLE IF EXISTS study_make30_final.eicu_sepsis_index;
DROP TABLE IF EXISTS study_make30_final.eicu_sepsis_adult;
DROP TABLE IF EXISTS study_make30_final.eicu_sepsis_first_stay;
DROP TABLE IF EXISTS study_make30_final.eicu_sepsis_pairs;
DROP TABLE IF EXISTS study_make30_final.eicu_sofa_standardized;
DROP TABLE IF EXISTS study_make30_final.eicu_sofa_hourly_components;
DROP TABLE IF EXISTS study_make30_final.eicu_suspected_infection;
DROP TABLE IF EXISTS study_make30_final.eicu_antibiotic_admin;
DROP TABLE IF EXISTS study_make30_final.eicu_antibiotic_patterns;
DROP TABLE IF EXISTS study_make30_final.eicu_icustays;
DROP TABLE IF EXISTS study_make30_final.eicu_patients;


/* ========================================================================
   1. FULL UNIQUE-PATIENT COHORT

   Unit = uniquepid
   All ICU stays belonging to all patients are retained.
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_patients AS

SELECT DISTINCT
    p.uniquepid

FROM public.patient p

WHERE p.uniquepid IS NOT NULL;


CREATE UNIQUE INDEX idx_eicu_patients_uniquepid
ON study_make30_final.eicu_patients(uniquepid);


/* ========================================================================
   2. ALL ICU STAYS FOR THE FULL PATIENT COHORT
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_icustays AS

SELECT
    p.uniquepid,
    p.patienthealthsystemstayid,
    p.patientunitstayid,

    p.gender,

    CASE
        WHEN BTRIM(p.age) ~ '^[0-9]+$'
        THEN BTRIM(p.age)::integer

        WHEN BTRIM(p.age) ~ '^>'
        THEN 90

        ELSE NULL
    END AS age_num,

    p.ethnicity,

    p.hospitalid,
    p.unitvisitnumber,
    p.unittype,
    p.unitstaytype,

    p.hospitaladmitoffset,
    p.hospitaldischargeoffset,
    p.hospitaldischargestatus,

    0::integer AS unitadmitoffset,
    p.unitdischargeoffset,
    p.unitdischargestatus

FROM public.patient p

INNER JOIN study_make30_final.eicu_patients s
    ON p.uniquepid = s.uniquepid;


CREATE UNIQUE INDEX idx_eicu_icustays_stay
ON study_make30_final.eicu_icustays(patientunitstayid);

CREATE INDEX idx_eicu_icustays_patient
ON study_make30_final.eicu_icustays(uniquepid);


/* ========================================================================
   3. EXPLICIT 25-PATTERN ANTIBIOTIC DICTIONARY

   Transparent database-specific mapping for the antibiotic-administration
   suspected-infection proxy.

   The pattern is matched case-insensitively against medication.drugname.
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_antibiotic_patterns
(
    pattern_id      integer PRIMARY KEY,
    pattern_regex   text NOT NULL,
    pattern_label   text NOT NULL
);

INSERT INTO study_make30_final.eicu_antibiotic_patterns
(pattern_id, pattern_regex, pattern_label)
VALUES
( 1, 'vancomycin|vancocin',                         'vancomycin'),
( 2, 'piperacillin|tazobactam|zosyn',              'piperacillin-tazobactam'),
( 3, 'cefepime|maxipime',                          'cefepime'),
( 4, 'ceftriaxone|rocephin',                       'ceftriaxone'),
( 5, 'ceftazidime|fortaz',                         'ceftazidime'),
( 6, 'cefotaxime|claforan',                        'cefotaxime'),
( 7, 'cefazolin|ancef',                            'cefazolin'),
( 8, 'cefoxitin|mefoxin',                          'cefoxitin'),
( 9, 'cefuroxime|zinacef|ceftin',                  'cefuroxime'),
(10, 'meropenem|merrem',                           'meropenem'),
(11, 'imipenem|cilastatin|primaxin',               'imipenem-cilastatin'),
(12, 'ertapenem|invanz',                           'ertapenem'),
(13, 'aztreonam|azactam',                          'aztreonam'),
(14, 'ampicillin',                                 'ampicillin'),
(15, 'amoxicillin',                                'amoxicillin'),
(16, '(^|[^a-z])penicillin([^a-z]|$)',             'penicillin'),
(17, 'oxacillin|nafcillin',                        'oxacillin-nafcillin'),
(18, 'ciprofloxacin|cipro',                        'ciprofloxacin'),
(19, 'levofloxacin|levaquin',                      'levofloxacin'),
(20, 'moxifloxacin|avelox',                        'moxifloxacin'),
(21, 'azithromycin|zithromax|clarithromycin|erythromycin',
                                                    'macrolide'),
(22, 'clindamycin|cleocin',                        'clindamycin'),
(23, 'metronidazole|flagyl',                       'metronidazole'),
(24, 'linezolid|zyvox|daptomycin|cubicin',         'linezolid-daptomycin'),
(25, 'gentamicin|tobramycin|amikacin',             'aminoglycoside');


/* ========================================================================
   4. RECOGNIZED ANTIBIOTIC ADMINISTRATION / START EVENTS

   ORAKLE public description for eICU:
      because fluid-sampling information is limited, suspicion of infection
      is defined using the OCCURRENCE OF ANTIBIOTIC ADMINISTRATIONS.

   Therefore:
      - one recognized antibiotic event is sufficient to create a candidate
        suspected-infection anchor;
      - there is NO >=2-dose / >=2-administration requirement;
      - route is NOT used as an additional inclusion restriction here,
        because the published eICU-specific ORAKLE wording does not state
        such a restriction.

   Medication fields used:
      patientunitstayid
      drugstartoffset
      drugordercancelled
      drugname
      routeadmin
      drugivadmixture

   eICU ETL note:
      drugstartoffset = 0 can represent unavailable/imprecise timing in some
      medication rows, so zero is excluded from the time-anchored analysis.

   Duplicate rows at the same stay + offset + normalized drug name are
   removed.

   The route fields are retained for auditing only; they do NOT determine
   eligibility for the suspected-infection anchor in this version.
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_antibiotic_admin AS

SELECT DISTINCT
    m.patientunitstayid,

    m.drugstartoffset::integer
        AS antibiotic_offset,

    LOWER(BTRIM(m.drugname))
        AS drugname_normalized,

    NULLIF(BTRIM(m.routeadmin), '')
        AS routeadmin,

    NULLIF(BTRIM(m.drugivadmixture), '')
        AS drugivadmixture

FROM public.medication m

INNER JOIN study_make30_final.eicu_icustays s
    ON m.patientunitstayid = s.patientunitstayid

WHERE m.drugstartoffset IS NOT NULL
  AND m.drugstartoffset <> 0

  AND COALESCE(m.drugordercancelled, 'No') = 'No'

  AND m.drugname IS NOT NULL

  AND EXISTS
  (
      SELECT 1

      FROM study_make30_final.eicu_antibiotic_patterns ptn

      WHERE LOWER(m.drugname) ~ ptn.pattern_regex
  );


CREATE INDEX idx_eicu_abx_stay_time
ON study_make30_final.eicu_antibiotic_admin
(patientunitstayid, antibiotic_offset);


/* ========================================================================
   5. SUSPECTED INFECTION — ORAKLE ANTIBIOTIC-ADMINISTRATION PROXY

   Published ORAKLE eICU/SICdb wording:
      occurrence of antibiotic administrations is used to define suspicion
      of infection when fluid-sampling information is limited.

   Operationalization here:
      every unique recognized antibiotic administration/start time is a
      candidate suspected-infection anchor.

   antibiotic_admin_count_24h is retained only as a DESCRIPTIVE field:
      it counts recognized antibiotic events from the anchor through the
      following 24 h.

   IMPORTANT:
      antibiotic_admin_count_24h is NOT an inclusion criterion.

   Multiple suspected-infection anchors may exist in one ICU stay.
   The ORAKLE SOFA window is applied to every candidate anchor, and the
   earliest qualifying Sepsis event is selected downstream.
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_suspected_infection AS

SELECT
    a.patientunitstayid,

    a.antibiotic_offset::integer
        AS suspected_infection_offset,

    (
        SELECT COUNT(*)::integer

        FROM study_make30_final.eicu_antibiotic_admin b

        WHERE b.patientunitstayid = a.patientunitstayid

          AND b.antibiotic_offset >= a.antibiotic_offset

          AND b.antibiotic_offset <= a.antibiotic_offset + 1440
    ) AS antibiotic_admin_count_24h

FROM
(
    SELECT DISTINCT
        patientunitstayid,
        antibiotic_offset

    FROM study_make30_final.eicu_antibiotic_admin
) a;


CREATE UNIQUE INDEX idx_eicu_si_stay_time
ON study_make30_final.eicu_suspected_infection
(patientunitstayid, suspected_infection_offset);


/* ========================================================================
   6. BUILD STANDARD HOURLY SOFA DIRECTLY FROM RAW eICU TABLES

   RATIONALE
   ---------
   eICU-CRD does not natively contain a ready-made hourly SOFA table.
   Therefore the score is reconstructed from the raw source tables.

   STANDARD ORIGINAL SOFA COMPONENTS
   ---------------------------------
   Respiratory:
      PaO2/FiO2, with invasive respiratory support for scores 3–4

   Coagulation:
      platelets

   Liver:
      total bilirubin

   Cardiovascular:
      MAP, dopamine, dobutamine, epinephrine, norepinephrine

   CNS:
      Glasgow Coma Scale

   Renal:
      creatinine and 24-h urine output

   The component thresholds are the original Vincent et al. 1996 SOFA
   thresholds. Component scores are calculated at hourly resolution; the
   maximum component score in the current + preceding 23 hourly rows is
   used for the rolling 24-h SOFA score.

   eICU offsets are minutes relative to ICU admission.
   ======================================================================== */


/* ------------------------------------------------------------------------
   6A. COMPLETE HOURLY GRID FOR ALL SAMPLED ICU STAYS

   hr = 0 corresponds to [0, 60) min.
   sofa_offset later uses the END of each hour.
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_hour_grid;

CREATE TEMP TABLE tmp_eicu_hour_grid AS

SELECT
    s.patientunitstayid,
    gs.hr::integer AS hr,

    (gs.hr * 60)::integer
        AS startoffset,

    ((gs.hr + 1) * 60)::integer
        AS endoffset

FROM study_make30_final.eicu_icustays s

CROSS JOIN LATERAL
    generate_series
    (
        0,

        GREATEST
        (
            0,

            CEIL(
                GREATEST(
                    COALESCE(s.unitdischargeoffset, 0),
                    1
                ) / 60.0
            )::integer - 1
        )
    ) AS gs(hr);


CREATE UNIQUE INDEX idx_tmp_eicu_hour_grid
ON tmp_eicu_hour_grid(patientunitstayid, hr);


/* ------------------------------------------------------------------------
   6B. GCS EVENTS FROM NURSECHARTING

   eICU concept labels:
      Glasgow coma score / GCS Total
      Score (Glasgow Coma Scale) / Value
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_gcs_event;

CREATE TEMP TABLE tmp_eicu_gcs_event AS

SELECT
    n.patientunitstayid,

    n.nursingchartoffset::integer
        AS chartoffset,

    CASE
        WHEN BTRIM(n.nursingchartvalue) ~
             '^[0-9]+([.][0-9]+)?$'

        THEN BTRIM(n.nursingchartvalue)::numeric

        ELSE NULL
    END AS gcs

FROM public.nursecharting n

INNER JOIN study_make30_final.eicu_icustays s
    ON n.patientunitstayid = s.patientunitstayid

WHERE n.nursingchartoffset >= 0

  AND
  (
      (
          n.nursingchartcelltypevallabel = 'Glasgow coma score'
          AND n.nursingchartcelltypevalname = 'GCS Total'
      )

      OR

      (
          n.nursingchartcelltypevallabel = 'Score (Glasgow Coma Scale)'
          AND n.nursingchartcelltypevalname = 'Value'
      )
  )

  AND BTRIM(n.nursingchartvalue) ~
      '^[0-9]+([.][0-9]+)?$'

  AND BTRIM(n.nursingchartvalue)::numeric
      BETWEEN 3 AND 15;


CREATE INDEX idx_tmp_eicu_gcs_event
ON tmp_eicu_gcs_event(patientunitstayid, chartoffset);


/* ------------------------------------------------------------------------
   6C. GCS -> HOURLY CNS SCORE
   Worst score within each hour.
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_cns_hour;

CREATE TEMP TABLE tmp_eicu_cns_hour AS

SELECT
    patientunitstayid,

    FLOOR(chartoffset / 60.0)::integer
        AS hr,

    MAX(
        CASE
            WHEN gcs < 6
            THEN 4

            WHEN gcs BETWEEN 6 AND 9
            THEN 3

            WHEN gcs BETWEEN 10 AND 12
            THEN 2

            WHEN gcs BETWEEN 13 AND 14
            THEN 1

            WHEN gcs = 15
            THEN 0

            ELSE NULL
        END
    )::integer AS cns_score_hour

FROM tmp_eicu_gcs_event

GROUP BY
    patientunitstayid,
    FLOOR(chartoffset / 60.0)::integer;


CREATE INDEX idx_tmp_eicu_cns_hour
ON tmp_eicu_cns_hour(patientunitstayid, hr);


/* ------------------------------------------------------------------------
   6D. PLATELET / BILIRUBIN / CREATININE EVENTS
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_sofa_lab_event;

CREATE TEMP TABLE tmp_eicu_sofa_lab_event AS

SELECT
    l.patientunitstayid,

    l.labresultoffset::integer
        AS chartoffset,

    LOWER(BTRIM(l.labname))
        AS labname,

    l.labresult::numeric
        AS labresult

FROM public.lab l

INNER JOIN study_make30_final.eicu_icustays s
    ON l.patientunitstayid = s.patientunitstayid

WHERE l.labresultoffset >= 0

  AND l.labresult IS NOT NULL

  AND LOWER(BTRIM(l.labname)) IN
  (
      'platelets x 1000',
      'platelets',
      'total bilirubin',
      'bilirubin',
      'creatinine'
  );


CREATE INDEX idx_tmp_eicu_sofa_lab_event
ON tmp_eicu_sofa_lab_event
(patientunitstayid, chartoffset, labname);


/* ------------------------------------------------------------------------
   6E. LAB EVENTS -> HOURLY COAGULATION / LIVER / CREATININE SCORES
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_lab_score_hour;

CREATE TEMP TABLE tmp_eicu_lab_score_hour AS

SELECT
    patientunitstayid,

    FLOOR(chartoffset / 60.0)::integer
        AS hr,

    MAX(
        CASE
            WHEN labname IN
                 ('platelets x 1000', 'platelets')

            THEN
                CASE
                    WHEN labresult < 20
                    THEN 4

                    WHEN labresult < 50
                    THEN 3

                    WHEN labresult < 100
                    THEN 2

                    WHEN labresult < 150
                    THEN 1

                    ELSE 0
                END

            ELSE NULL
        END
    )::integer AS coag_score_hour,

    MAX(
        CASE
            WHEN labname IN
                 ('total bilirubin', 'bilirubin')

            THEN
                CASE
                    WHEN labresult >= 12.0
                    THEN 4

                    WHEN labresult >= 6.0
                    THEN 3

                    WHEN labresult >= 2.0
                    THEN 2

                    WHEN labresult >= 1.2
                    THEN 1

                    ELSE 0
                END

            ELSE NULL
        END
    )::integer AS liver_score_hour,

    MAX(
        CASE
            WHEN labname = 'creatinine'

            THEN
                CASE
                    WHEN labresult >= 5.0
                    THEN 4

                    WHEN labresult >= 3.5
                    THEN 3

                    WHEN labresult >= 2.0
                    THEN 2

                    WHEN labresult >= 1.2
                    THEN 1

                    ELSE 0
                END

            ELSE NULL
        END
    )::integer AS creatinine_score_hour

FROM tmp_eicu_sofa_lab_event

GROUP BY
    patientunitstayid,
    FLOOR(chartoffset / 60.0)::integer;


CREATE INDEX idx_tmp_eicu_lab_score_hour
ON tmp_eicu_lab_score_hour(patientunitstayid, hr);


/* ------------------------------------------------------------------------
   6F. PaO2 / FiO2 EVENTS

   Duplicate revised labs at the same time are resolved by keeping the
   latest revised result.

   FiO2 normalization:
      0.20–1.00 -> fraction already
      20–100    -> percent, divide by 100
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_bg;

CREATE TEMP TABLE tmp_eicu_bg AS

WITH ranked_lab AS
(
    SELECT
        l.patientunitstayid,
        l.labresultoffset,
        l.labname,
        l.labresult,

        ROW_NUMBER() OVER
        (
            PARTITION BY
                l.patientunitstayid,
                l.labname,
                l.labresultoffset

            ORDER BY
                l.labresultrevisedoffset DESC NULLS LAST,
                l.labresult DESC NULLS LAST
        ) AS rn

    FROM public.lab l

    INNER JOIN study_make30_final.eicu_icustays s
        ON l.patientunitstayid = s.patientunitstayid

    WHERE l.labresultoffset >= 0
      AND l.labresult IS NOT NULL
      AND l.labname IN ('paO2', 'FiO2')
),

last_lab AS
(
    SELECT
        patientunitstayid,
        labresultoffset,
        labname,
        labresult

    FROM ranked_lab

    WHERE rn = 1
),

pivoted AS
(
    SELECT
        patientunitstayid,
        labresultoffset,

        MAX(labresult)
            FILTER (WHERE labname = 'paO2')
            AS pao2,

        MAX(labresult)
            FILTER (WHERE labname = 'FiO2')
            AS fio2_raw

    FROM last_lab

    GROUP BY
        patientunitstayid,
        labresultoffset
)

SELECT
    patientunitstayid,

    labresultoffset::integer
        AS chartoffset,

    pao2::numeric,

    CASE
        WHEN fio2_raw BETWEEN 0.20 AND 1.00
        THEN fio2_raw

        WHEN fio2_raw BETWEEN 20 AND 100
        THEN fio2_raw / 100.0

        ELSE NULL
    END::numeric AS fio2

FROM pivoted

WHERE pao2 BETWEEN 15 AND 720;


CREATE INDEX idx_tmp_eicu_bg
ON tmp_eicu_bg(patientunitstayid, chartoffset);


/* ------------------------------------------------------------------------
   6G. INVASIVE VENTILATION INTERVALS

   Uses eICU RESPIRATORYCARE vent start/end offsets.
   Current and prior ventilation intervals are retained when valid.
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_vent_interval;

CREATE TEMP TABLE tmp_eicu_vent_interval AS

SELECT DISTINCT
    r.patientunitstayid,

    r.ventstartoffset::integer
        AS startoffset,

    r.ventendoffset::integer
        AS endoffset

FROM public.respiratorycare r

INNER JOIN study_make30_final.eicu_icustays s
    ON r.patientunitstayid = s.patientunitstayid

WHERE r.ventstartoffset IS NOT NULL
  AND r.ventendoffset IS NOT NULL
  AND r.ventendoffset > r.ventstartoffset


UNION


SELECT DISTINCT
    r.patientunitstayid,

    r.priorventstartoffset::integer,

    r.priorventendoffset::integer

FROM public.respiratorycare r

INNER JOIN study_make30_final.eicu_icustays s
    ON r.patientunitstayid = s.patientunitstayid

WHERE r.priorventstartoffset IS NOT NULL
  AND r.priorventendoffset IS NOT NULL
  AND r.priorventendoffset > r.priorventstartoffset;


CREATE INDEX idx_tmp_eicu_vent_interval
ON tmp_eicu_vent_interval
(patientunitstayid, startoffset, endoffset);


/* ------------------------------------------------------------------------
   6H. PaO2/FiO2 -> HOURLY RESPIRATORY SCORE
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_resp_hour;

CREATE TEMP TABLE tmp_eicu_resp_hour AS

SELECT
    bg.patientunitstayid,

    FLOOR(bg.chartoffset / 60.0)::integer
        AS hr,

    MAX(
        CASE
            WHEN bg.fio2 IS NULL
              OR bg.fio2 <= 0
            THEN NULL

            WHEN (bg.pao2 / bg.fio2) < 100

             AND EXISTS
             (
                 SELECT 1
                 FROM tmp_eicu_vent_interval v
                 WHERE v.patientunitstayid = bg.patientunitstayid
                   AND bg.chartoffset >= v.startoffset
                   AND bg.chartoffset <= v.endoffset
             )

            THEN 4

            WHEN (bg.pao2 / bg.fio2) < 200

             AND EXISTS
             (
                 SELECT 1
                 FROM tmp_eicu_vent_interval v
                 WHERE v.patientunitstayid = bg.patientunitstayid
                   AND bg.chartoffset >= v.startoffset
                   AND bg.chartoffset <= v.endoffset
             )

            THEN 3

            WHEN (bg.pao2 / bg.fio2) < 300
            THEN 2

            WHEN (bg.pao2 / bg.fio2) < 400
            THEN 1

            ELSE 0
        END
    )::integer AS resp_score_hour

FROM tmp_eicu_bg bg

WHERE bg.fio2 IS NOT NULL

GROUP BY
    bg.patientunitstayid,
    FLOOR(bg.chartoffset / 60.0)::integer;


CREATE INDEX idx_tmp_eicu_resp_hour
ON tmp_eicu_resp_hour(patientunitstayid, hr);


/* ------------------------------------------------------------------------
   6I. MAP EVENTS -> HOURLY MAP SCORE

   invasive:
      vitalperiodic.systemicmean

   non-invasive:
      vitalaperiodic.noninvasivemean
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_map_hour;

CREATE TEMP TABLE tmp_eicu_map_hour AS

WITH map_events AS
(
    SELECT
        v.patientunitstayid,

        v.observationoffset::integer
            AS chartoffset,

        v.systemicmean::numeric
            AS map_value

    FROM public.vitalperiodic v

    INNER JOIN study_make30_final.eicu_icustays s
        ON v.patientunitstayid = s.patientunitstayid

    WHERE v.observationoffset >= 0
      AND v.systemicmean > 0
      AND v.systemicmean < 300


    UNION ALL


    SELECT
        v.patientunitstayid,

        v.observationoffset::integer,

        v.noninvasivemean::numeric

    FROM public.vitalaperiodic v

    INNER JOIN study_make30_final.eicu_icustays s
        ON v.patientunitstayid = s.patientunitstayid

    WHERE v.observationoffset >= 0
      AND v.noninvasivemean > 0
      AND v.noninvasivemean < 300
)

SELECT
    patientunitstayid,

    FLOOR(chartoffset / 60.0)::integer
        AS hr,

    MAX(
        CASE
            WHEN map_value < 70
            THEN 1
            ELSE 0
        END
    )::integer AS map_score_hour

FROM map_events

GROUP BY
    patientunitstayid,
    FLOOR(chartoffset / 60.0)::integer;


CREATE INDEX idx_tmp_eicu_map_hour
ON tmp_eicu_map_hour(patientunitstayid, hr);


/* ------------------------------------------------------------------------
   6J. DOSE-CODED PRESSOR / INOTROPE TREATMENT SCORE

   eICU treatment strings explicitly encode the original SOFA dose
   categories for dopamine, norepinephrine and epinephrine.

   Original SOFA:
      score 4:
         dopamine >15
         OR epinephrine >0.1
         OR norepinephrine >0.1

      score 3:
         dopamine >5 to <=15
         OR epinephrine <=0.1
         OR norepinephrine <=0.1

      score 2:
         dopamine >0 to <=5
         OR dobutamine any dose
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_cvdrug_hour;

CREATE TEMP TABLE tmp_eicu_cvdrug_hour AS

WITH treatment_score AS
(
    SELECT
        t.patientunitstayid,
        t.treatmentoffset,

        CASE
            WHEN LOWER(t.treatmentstring) LIKE
                 '%norepinephrine > 0.1 micrograms/kg/min%'

              OR LOWER(t.treatmentstring) LIKE
                 '%epinephrine > 0.1 micrograms/kg/min%'

              OR LOWER(t.treatmentstring) LIKE
                 '%dopamine >15 micrograms/kg/min%'

              OR LOWER(t.treatmentstring) LIKE
                 '%dopamine > 15 micrograms/kg/min%'

            THEN 4


            WHEN LOWER(t.treatmentstring) LIKE
                 '%norepinephrine <= 0.1 micrograms/kg/min%'

              OR LOWER(t.treatmentstring) LIKE
                 '%epinephrine <= 0.1 micrograms/kg/min%'

              OR LOWER(t.treatmentstring) LIKE
                 '%dopamine 5-15 micrograms/kg/min%'

            THEN 3


            WHEN LOWER(t.treatmentstring) LIKE
                 '%dopamine <= 5 micrograms/kg/min%'

              OR LOWER(t.treatmentstring) LIKE
                 '%dobutamine%'

            THEN 2

            ELSE NULL
        END AS cv_score

    FROM public.treatment t

    INNER JOIN study_make30_final.eicu_icustays s
        ON t.patientunitstayid = s.patientunitstayid

    WHERE t.treatmentoffset >= 0
      AND t.treatmentstring IS NOT NULL

      AND LOWER(t.treatmentstring) ~
          '(dopamine|dobutamine|norepinephrine|epinephrine)'
),

med_dobutamine AS
(
    SELECT
        m.patientunitstayid,

        FLOOR(m.drugstartoffset / 60.0)::integer
            AS hr,

        2::integer AS cv_score

    FROM public.medication m

    INNER JOIN study_make30_final.eicu_icustays s
        ON m.patientunitstayid = s.patientunitstayid

    WHERE m.drugstartoffset IS NOT NULL
      AND m.drugstartoffset > 0

      AND COALESCE(m.drugordercancelled, 'No') = 'No'

      AND
      (
          m.drughiclseqno IN (8777, 40)
          OR LOWER(COALESCE(m.drugname,'')) LIKE '%dobutamine%'
      )
),

all_scores AS
(
    SELECT
        patientunitstayid,

        FLOOR(treatmentoffset / 60.0)::integer
            AS hr,

        cv_score

    FROM treatment_score

    WHERE cv_score IS NOT NULL


    UNION ALL


    SELECT
        patientunitstayid,
        hr,
        cv_score

    FROM med_dobutamine
)

SELECT
    patientunitstayid,
    hr,

    MAX(cv_score)::integer
        AS cv_drug_score_hour

FROM all_scores

GROUP BY
    patientunitstayid,
    hr;


CREATE INDEX idx_tmp_eicu_cvdrug_hour
ON tmp_eicu_cvdrug_hour(patientunitstayid, hr);


/* ------------------------------------------------------------------------
   6K. URINE OUTPUT EVENTS

   This uses the eICU I/O output hierarchy and urinary labels.
   Dialysis / dialysate output is explicitly excluded.

   IMPORTANT:
      no-record is kept as NULL, not as zero.
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_urine_event;

CREATE TEMP TABLE tmp_eicu_urine_event AS

SELECT
    io.patientunitstayid,

    io.intakeoutputoffset::integer
        AS chartoffset,

    SUM(io.cellvaluenumeric)::numeric
        AS urine_ml

FROM public.intakeoutput io

INNER JOIN study_make30_final.eicu_icustays s
    ON io.patientunitstayid = s.patientunitstayid

WHERE io.intakeoutputoffset >= 0

  AND io.cellvaluenumeric IS NOT NULL
  AND io.cellvaluenumeric >= 0

  AND LOWER(COALESCE(io.cellpath,'')) LIKE
      'flowsheet|flowsheet cell labels|i&o|output (ml)|%'

  AND
  (
      LOWER(COALESCE(io.cellpath,'')) LIKE '%urine%'

      OR LOWER(COALESCE(io.cellpath,'')) LIKE '%foley%'

      OR LOWER(COALESCE(io.cellpath,'')) LIKE
         '%urethral catheter%'

      OR LOWER(COALESCE(io.cellpath,'')) LIKE
         '%urinary catheter%'

      OR LOWER(COALESCE(io.cellpath,'')) LIKE
         '%straight cath%'

      OR LOWER(COALESCE(io.cellpath,'')) LIKE
         '%condom cath%'

      OR LOWER(COALESCE(io.cellpath,'')) LIKE
         '%nephrostomy%'
  )

  AND LOWER(COALESCE(io.cellpath,'')) NOT LIKE
      '%dialysis%'

  AND LOWER(COALESCE(io.cellpath,'')) NOT LIKE
      '%dialysate%'

GROUP BY
    io.patientunitstayid,
    io.intakeoutputoffset;


CREATE INDEX idx_tmp_eicu_urine_event
ON tmp_eicu_urine_event(patientunitstayid, chartoffset);


/* ------------------------------------------------------------------------
   6L. URINE -> HOURLY VOLUME
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_urine_hour;

CREATE TEMP TABLE tmp_eicu_urine_hour AS

SELECT
    patientunitstayid,

    FLOOR(chartoffset / 60.0)::integer
        AS hr,

    SUM(urine_ml)::numeric
        AS urine_ml_hour

FROM tmp_eicu_urine_event

GROUP BY
    patientunitstayid,
    FLOOR(chartoffset / 60.0)::integer;


CREATE INDEX idx_tmp_eicu_urine_hour
ON tmp_eicu_urine_hour(patientunitstayid, hr);


/* ------------------------------------------------------------------------
   6M. ASSEMBLE RAW HOURLY COMPONENTS

   Missing component-specific score remains NULL at this stage.
   ------------------------------------------------------------------------ */

DROP TABLE IF EXISTS tmp_eicu_sofa_raw_hour;

CREATE TEMP TABLE tmp_eicu_sofa_raw_hour AS

SELECT
    g.patientunitstayid,
    g.hr,
    g.startoffset,
    g.endoffset,

    r.resp_score_hour
        AS respiration_hour,

    l.coag_score_hour
        AS coagulation_hour,

    l.liver_score_hour
        AS liver_hour,

    CASE
        WHEN m.map_score_hour IS NULL
         AND v.cv_drug_score_hour IS NULL
        THEN NULL

        ELSE GREATEST(
            COALESCE(m.map_score_hour, 0),
            COALESCE(v.cv_drug_score_hour, 0)
        )
    END AS cardiovascular_hour,

    c.cns_score_hour
        AS cns_hour,

    l.creatinine_score_hour
        AS creatinine_hour,

    u.urine_ml_hour
        AS urine_ml_hour

FROM tmp_eicu_hour_grid g

LEFT JOIN tmp_eicu_resp_hour r
    ON g.patientunitstayid = r.patientunitstayid
   AND g.hr = r.hr

LEFT JOIN tmp_eicu_lab_score_hour l
    ON g.patientunitstayid = l.patientunitstayid
   AND g.hr = l.hr

LEFT JOIN tmp_eicu_map_hour m
    ON g.patientunitstayid = m.patientunitstayid
   AND g.hr = m.hr

LEFT JOIN tmp_eicu_cvdrug_hour v
    ON g.patientunitstayid = v.patientunitstayid
   AND g.hr = v.hr

LEFT JOIN tmp_eicu_cns_hour c
    ON g.patientunitstayid = c.patientunitstayid
   AND g.hr = c.hr

LEFT JOIN tmp_eicu_urine_hour u
    ON g.patientunitstayid = u.patientunitstayid
   AND g.hr = u.hr;


CREATE UNIQUE INDEX idx_tmp_eicu_sofa_raw_hour
ON tmp_eicu_sofa_raw_hour(patientunitstayid, hr);


/* ------------------------------------------------------------------------
   6N. ROLLING 24-H COMPONENTS

   For score components:
      maximum over current + preceding 23 hourly rows.

   Urine:
      sum only actual observed urine-output records.
      Missing urine records are NOT converted to zero.
   ------------------------------------------------------------------------ */

CREATE TABLE study_make30_final.eicu_sofa_hourly_components AS

SELECT
    patientunitstayid,
    hr,
    startoffset,
    endoffset,

    MAX(respiration_hour)
        OVER w24
        AS respiration_24hours_raw,

    MAX(coagulation_hour)
        OVER w24
        AS coagulation_24hours_raw,

    MAX(liver_hour)
        OVER w24
        AS liver_24hours_raw,

    MAX(cardiovascular_hour)
        OVER w24
        AS cardiovascular_24hours_raw,

    MAX(cns_hour)
        OVER w24
        AS cns_24hours_raw,

    MAX(creatinine_hour)
        OVER w24
        AS creatinine_24hours_raw,

    SUM(urine_ml_hour)
        OVER w24
        AS urine_24hours,

    COUNT(urine_ml_hour)
        OVER w24
        AS urine_observed_hours

FROM tmp_eicu_sofa_raw_hour

WINDOW w24 AS
(
    PARTITION BY patientunitstayid

    ORDER BY hr

    ROWS BETWEEN 23 PRECEDING AND CURRENT ROW
);


CREATE UNIQUE INDEX idx_eicu_sofa_hourly_components
ON study_make30_final.eicu_sofa_hourly_components
(patientunitstayid, hr);


/* ------------------------------------------------------------------------
   6O. FINAL STANDARD HOURLY 24-H SOFA

   Missing component scores are imputed to 0 points at the final rolling
   component stage, matching common ICU-database SOFA implementations.

   Renal score:
      creatinine criterion can contribute at all hours.
      urine criterion is used only when:
         a full 24-hour ICU window exists (hr >=23)
         AND at least one urine measurement exists in that window.

   sofa_offset:
      end of the hourly interval, in minutes from ICU admission.
   ------------------------------------------------------------------------ */

CREATE TABLE study_make30_final.eicu_sofa_standardized AS

WITH component_scores AS
(
    SELECT
        patientunitstayid,
        hr,
        startoffset,
        endoffset,

        COALESCE(
            respiration_24hours_raw,
            0
        )::integer AS respiration,

        COALESCE(
            coagulation_24hours_raw,
            0
        )::integer AS coagulation,

        COALESCE(
            liver_24hours_raw,
            0
        )::integer AS liver,

        COALESCE(
            cardiovascular_24hours_raw,
            0
        )::integer AS cardiovascular,

        COALESCE(
            cns_24hours_raw,
            0
        )::integer AS cns,

        GREATEST
        (
            COALESCE(
                creatinine_24hours_raw,
                0
            ),

            CASE
                WHEN hr < 23
                THEN 0

                WHEN urine_observed_hours <= 0
                THEN 0

                WHEN urine_24hours < 200
                THEN 4

                WHEN urine_24hours < 500
                THEN 3

                ELSE 0
            END
        )::integer AS renal

    FROM study_make30_final.eicu_sofa_hourly_components
)

SELECT
    patientunitstayid,
    hr,

    endoffset::integer
        AS sofa_offset,

    respiration,
    coagulation,
    liver,
    cardiovascular,
    cns,
    renal,

    (
        respiration
        + coagulation
        + liver
        + cardiovascular
        + cns
        + renal
    )::integer AS sofa_score

FROM component_scores;


CREATE UNIQUE INDEX idx_eicu_sofa_standardized
ON study_make30_final.eicu_sofa_standardized
(patientunitstayid, sofa_offset);

CREATE INDEX idx_eicu_sofa_ge2
ON study_make30_final.eicu_sofa_standardized
(patientunitstayid, sofa_score, sofa_offset);


/* ========================================================================
   7. ORAKLE-STYLE SUSPECTED INFECTION + SOFA MATCHING

   qualifying sofa_offset:
       suspected infection - 24 h
       through
       suspected infection + 12 h

   In minutes:
       -1440
       +720

   sepsis_offset:
       earlier of suspected infection and qualifying SOFA
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_sepsis_pairs AS

SELECT
    s.uniquepid,
    s.patienthealthsystemstayid,
    s.patientunitstayid,

    si.suspected_infection_offset,
    si.antibiotic_admin_count_24h,

    sofa.sofa_offset,
    sofa.sofa_score,

    LEAST(
        si.suspected_infection_offset,
        sofa.sofa_offset
    )::integer AS sepsis_offset

FROM study_make30_final.eicu_icustays s

INNER JOIN study_make30_final.eicu_suspected_infection si
    ON s.patientunitstayid = si.patientunitstayid

INNER JOIN study_make30_final.eicu_sofa_standardized sofa
    ON s.patientunitstayid = sofa.patientunitstayid

   AND sofa.sofa_score >= 2

   AND sofa.sofa_offset >=
       si.suspected_infection_offset - 1440

   AND sofa.sofa_offset <=
       si.suspected_infection_offset + 720;


CREATE INDEX idx_eicu_sepsis_pairs_stay
ON study_make30_final.eicu_sepsis_pairs
(patientunitstayid, sepsis_offset);


/* ========================================================================
   8. FIRST QUALIFYING SEPSIS EVENT PER ICU STAY
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_sepsis_first_stay AS

SELECT
    uniquepid,
    patienthealthsystemstayid,
    patientunitstayid,

    suspected_infection_offset,
    antibiotic_admin_count_24h,

    sofa_offset,
    sofa_score,

    sepsis_offset

FROM
(
    SELECT
        p.*,

        ROW_NUMBER() OVER
        (
            PARTITION BY p.patientunitstayid

            ORDER BY
                p.sepsis_offset,
                p.suspected_infection_offset,
                p.sofa_offset
        ) AS rn_stay

    FROM study_make30_final.eicu_sepsis_pairs p
) q

WHERE rn_stay = 1;


/* ========================================================================
   9. ADULT SEPSIS WITH ONSET INSIDE ICU

   eICU ICU admission offset = 0.

   Require:
       age >=18
       sepsis_offset >=0
       sepsis_offset < unitdischargeoffset
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_sepsis_adult AS

SELECT
    s.uniquepid,
    s.patienthealthsystemstayid,
    s.patientunitstayid,

    p.gender,
    p.age_num AS age,
    p.ethnicity,

    p.hospitalid,
    p.unitvisitnumber,
    p.unittype,
    p.unitstaytype,

    p.hospitaladmitoffset,
    p.hospitaldischargeoffset,
    p.hospitaldischargestatus,

    p.unitadmitoffset,
    p.unitdischargeoffset,
    p.unitdischargestatus,

    s.suspected_infection_offset,
    s.antibiotic_admin_count_24h,

    s.sofa_offset,
    s.sofa_score,

    s.sepsis_offset,

    (s.sepsis_offset + 1440)::integer
        AS landmark_offset,

    (s.sepsis_offset + 43200)::integer
        AS sepsis_day30_offset

FROM study_make30_final.eicu_sepsis_first_stay s

INNER JOIN study_make30_final.eicu_icustays p
    ON s.patientunitstayid = p.patientunitstayid

WHERE p.age_num >= 18

  AND s.sepsis_offset >= 0

  AND p.unitdischargeoffset IS NOT NULL

  AND s.sepsis_offset < p.unitdischargeoffset;


/* ========================================================================
   10. FIRST QUALIFYING SEPSIS ICU STAY PER UNIQUE PATIENT
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_sepsis_index AS

SELECT
    uniquepid,
    patienthealthsystemstayid,
    patientunitstayid,

    gender,
    age,
    ethnicity,

    hospitalid,
    unitvisitnumber,
    unittype,
    unitstaytype,

    hospitaladmitoffset,
    hospitaldischargeoffset,
    hospitaldischargestatus,

    unitadmitoffset,
    unitdischargeoffset,
    unitdischargestatus,

    suspected_infection_offset,
    antibiotic_admin_count_24h,

    sofa_offset,
    sofa_score,

    sepsis_offset,
    landmark_offset,
    sepsis_day30_offset

FROM
(
    SELECT
        a.*,

        ROW_NUMBER() OVER
        (
            PARTITION BY a.uniquepid

            ORDER BY
                a.sepsis_offset,
                a.patienthealthsystemstayid,
                a.unitvisitnumber,
                a.patientunitstayid
        ) AS rn_patient

    FROM study_make30_final.eicu_sepsis_adult a
) q

WHERE rn_patient = 1;


CREATE UNIQUE INDEX idx_eicu_sepsis_index_stay
ON study_make30_final.eicu_sepsis_index(patientunitstayid);

CREATE UNIQUE INDEX idx_eicu_sepsis_index_patient
ON study_make30_final.eicu_sepsis_index(uniquepid);


/* ========================================================================
   11. +24 H LANDMARK COHORT

   Require:
      still in ICU at landmark
      and hospital discharge not before landmark
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_landmark_cohort AS

SELECT
    s.*,

    CASE
        WHEN s.unitdischargeoffset >= s.landmark_offset

         AND
         (
             s.hospitaldischargeoffset IS NULL
             OR
             s.hospitaldischargeoffset >= s.landmark_offset
         )

        THEN 1
        ELSE 0
    END AS reached_24h_landmark

FROM study_make30_final.eicu_sepsis_index s;


CREATE UNIQUE INDEX idx_eicu_landmark_stay
ON study_make30_final.eicu_landmark_cohort(patientunitstayid);


/* ========================================================================
   12. KRT DURING PHENOTYPE WINDOW [t0, t0+24h)

   Sources:
      public.treatment
      public.intakeoutput

   Treatment string signals:
      dialysis
      hemodialysis
      CRRT / CVVH / CVVHD / CVVHDF
      SLED
      peritoneal dialysis

   Because ANY KRT during the phenotype window materially changes K/Mg/P/
   bicarbonate and kidney physiology, the phenotype cohort excludes such
   treatment regardless of whether chronic or acute status can be proven
   from this single window.

   Formal post-landmark "new KRT" MAKE component will be constructed later.
   ======================================================================== */

DROP TABLE IF EXISTS tmp_eicu_early_krt;

CREATE TEMP TABLE tmp_eicu_early_krt AS

WITH treatment_krt AS
(
    SELECT DISTINCT
        c.patientunitstayid

    FROM study_make30_final.eicu_landmark_cohort c

    INNER JOIN public.treatment t
        ON c.patientunitstayid = t.patientunitstayid

       AND t.treatmentoffset >= c.sepsis_offset
       AND t.treatmentoffset <  c.landmark_offset

    WHERE c.reached_24h_landmark = 1

      AND LOWER(COALESCE(t.treatmentstring, '')) ~
          '(dialysis|hemodialysis|haemodialysis|crrt|cvvh|cvvhd|cvvhdf|sled|peritoneal)'
),

io_krt AS
(
    SELECT DISTINCT
        c.patientunitstayid

    FROM study_make30_final.eicu_landmark_cohort c

    INNER JOIN public.intakeoutput io
        ON c.patientunitstayid = io.patientunitstayid

       AND io.intakeoutputoffset >= c.sepsis_offset
       AND io.intakeoutputoffset <  c.landmark_offset

    WHERE c.reached_24h_landmark = 1

      AND
      (
          COALESCE(io.dialysistotal, 0) <> 0

          OR LOWER(COALESCE(io.celllabel, '')) LIKE '%dialysis%'

          OR LOWER(COALESCE(io.cellpath, '')) LIKE '%dialysis%'
      )
),

all_krt AS
(
    SELECT patientunitstayid
    FROM treatment_krt

    UNION

    SELECT patientunitstayid
    FROM io_krt
)

SELECT
    c.patientunitstayid,

    CASE
        WHEN k.patientunitstayid IS NOT NULL
        THEN 1
        ELSE 0
    END AS krt_0_24h

FROM study_make30_final.eicu_landmark_cohort c

LEFT JOIN all_krt k
    ON c.patientunitstayid = k.patientunitstayid

WHERE c.reached_24h_landmark = 1;


CREATE UNIQUE INDEX idx_tmp_eicu_early_krt
ON tmp_eicu_early_krt(patientunitstayid);


/* ========================================================================
   13. TRAJECTORY POOL AFTER LANDMARK + EARLY-KRT EXCLUSION
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_trajectory_pool AS

SELECT
    c.*,

    COALESCE(k.krt_0_24h, 0)
        AS krt_0_24h

FROM study_make30_final.eicu_landmark_cohort c

LEFT JOIN tmp_eicu_early_krt k
    ON c.patientunitstayid = k.patientunitstayid

WHERE c.reached_24h_landmark = 1

  AND COALESCE(k.krt_0_24h, 0) = 0;


CREATE UNIQUE INDEX idx_eicu_trajcand_stay
ON study_make30_final.eicu_trajectory_pool(patientunitstayid);


/* ========================================================================
   14. HR HOURLY MEDIAN
   ======================================================================== */

DROP TABLE IF EXISTS tmp_eicu_hr_hour;

CREATE TEMP TABLE tmp_eicu_hr_hour AS

SELECT
    b.patientunitstayid,

    FLOOR(
        (v.observationoffset - b.sepsis_offset) / 60.0
    )::integer AS hour_bin,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY v.heartrate)
        AS hr_median

FROM study_make30_final.eicu_trajectory_pool b

INNER JOIN public.vitalperiodic v
    ON b.patientunitstayid = v.patientunitstayid

   AND v.observationoffset >= b.sepsis_offset
   AND v.observationoffset <  b.landmark_offset

WHERE v.heartrate > 0
  AND v.heartrate < 300

GROUP BY
    b.patientunitstayid,

    FLOOR(
        (v.observationoffset - b.sepsis_offset) / 60.0
    )::integer;


CREATE INDEX idx_tmp_eicu_hr_hour
ON tmp_eicu_hr_hour(patientunitstayid, hour_bin);


/* ========================================================================
   15. INVASIVE MAP HOURLY MEDIAN
   ======================================================================== */

DROP TABLE IF EXISTS tmp_eicu_imap_hour;

CREATE TEMP TABLE tmp_eicu_imap_hour AS

SELECT
    b.patientunitstayid,

    FLOOR(
        (v.observationoffset - b.sepsis_offset) / 60.0
    )::integer AS hour_bin,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY v.systemicmean)
        AS map_invasive_median

FROM study_make30_final.eicu_trajectory_pool b

INNER JOIN public.vitalperiodic v
    ON b.patientunitstayid = v.patientunitstayid

   AND v.observationoffset >= b.sepsis_offset
   AND v.observationoffset <  b.landmark_offset

WHERE v.systemicmean > 20
  AND v.systemicmean < 300

GROUP BY
    b.patientunitstayid,

    FLOOR(
        (v.observationoffset - b.sepsis_offset) / 60.0
    )::integer;


CREATE INDEX idx_tmp_eicu_imap_hour
ON tmp_eicu_imap_hour(patientunitstayid, hour_bin);


/* ========================================================================
   16. NON-INVASIVE MAP HOURLY MEDIAN
   ======================================================================== */

DROP TABLE IF EXISTS tmp_eicu_nimap_hour;

CREATE TEMP TABLE tmp_eicu_nimap_hour AS

SELECT
    b.patientunitstayid,

    FLOOR(
        (v.observationoffset - b.sepsis_offset) / 60.0
    )::integer AS hour_bin,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY v.noninvasivemean)
        AS map_noninvasive_median

FROM study_make30_final.eicu_trajectory_pool b

INNER JOIN public.vitalaperiodic v
    ON b.patientunitstayid = v.patientunitstayid

   AND v.observationoffset >= b.sepsis_offset
   AND v.observationoffset <  b.landmark_offset

WHERE v.noninvasivemean > 20
  AND v.noninvasivemean < 300

GROUP BY
    b.patientunitstayid,

    FLOOR(
        (v.observationoffset - b.sepsis_offset) / 60.0
    )::integer;


CREATE INDEX idx_tmp_eicu_nimap_hour
ON tmp_eicu_nimap_hour(patientunitstayid, hour_bin);


/* ========================================================================
   17. VALID HOURLY MSI — UNIFIED MAP QC (FINAL <=20 RULE)

   HR:
      hourly median of the database primary monitored HR.

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

CREATE TABLE study_make30_final.eicu_msi_hourly AS

WITH selected AS
(
    SELECT
        h.patientunitstayid,
        h.hour_bin,
        h.hr_median,
        i.map_invasive_median,
        n.map_noninvasive_median,

        CASE
            WHEN i.map_invasive_median IS NOT NULL
             AND n.map_noninvasive_median IS NOT NULL
            THEN ABS(i.map_invasive_median - n.map_noninvasive_median)
            ELSE NULL
        END::double precision AS map_abs_diff,

        CASE
            WHEN i.map_invasive_median IS NOT NULL
             AND n.map_noninvasive_median IS NOT NULL
             AND ABS(i.map_invasive_median - n.map_noninvasive_median) >= 30
            THEN n.map_noninvasive_median

            WHEN i.map_invasive_median IS NOT NULL
            THEN i.map_invasive_median

            WHEN n.map_noninvasive_median IS NOT NULL
            THEN n.map_noninvasive_median

            ELSE NULL
        END::double precision AS map_used,

        CASE
            WHEN i.map_invasive_median IS NOT NULL
             AND n.map_noninvasive_median IS NOT NULL
             AND ABS(i.map_invasive_median - n.map_noninvasive_median) >= 30
            THEN 'noninvasive_fallback_discordant_ge30'

            WHEN i.map_invasive_median IS NOT NULL
            THEN 'invasive'

            WHEN n.map_noninvasive_median IS NOT NULL
            THEN 'noninvasive'

            ELSE NULL
        END::text AS map_source

    FROM tmp_eicu_hr_hour h

    LEFT JOIN tmp_eicu_imap_hour i
        ON h.patientunitstayid = i.patientunitstayid
       AND h.hour_bin = i.hour_bin

    LEFT JOIN tmp_eicu_nimap_hour n
        ON h.patientunitstayid = n.patientunitstayid
       AND h.hour_bin = n.hour_bin

    WHERE h.hour_bin BETWEEN 0 AND 23
)

SELECT
    patientunitstayid,
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

CREATE UNIQUE INDEX idx_eicu_msi_hourly
ON study_make30_final.eicu_msi_hourly
(patientunitstayid, hour_bin);


/* ========================================================================
   18. HOURLY MSI -> TWELVE 2-H MSI WINDOWS
   ======================================================================== */
CREATE TABLE study_make30_final.eicu_msi_2h AS
SELECT patientunitstayid,FLOOR(hour_bin/2.0)::integer AS time_bin,MAX(msi)::double precision AS msi,COUNT(msi)::integer AS valid_hourly_msi_count
FROM study_make30_final.eicu_msi_hourly WHERE msi IS NOT NULL
GROUP BY patientunitstayid,FLOOR(hour_bin/2.0)::integer;
CREATE UNIQUE INDEX idx_eicu_msi_2h ON study_make30_final.eicu_msi_2h(patientunitstayid,time_bin);

/* ========================================================================
   19. CORE-4 ELECTROLYTES -> FOUR 6-H WINDOWS

   eICU lab names:
      sodium
      potassium
      chloride
      bicarbonate

   Each 6-h value:
      median of available measurements.
   ======================================================================== */

CREATE TABLE study_make30_final.eicu_electrolyte_6h AS

SELECT
    b.patientunitstayid,

    FLOOR(
        (l.labresultoffset - b.sepsis_offset) / 360.0
    )::integer AS time_bin,

    LOWER(l.labname)
        AS variable_name,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY l.labresult)
        AS median_value,

    COUNT(*)::integer
        AS raw_measurement_count

FROM study_make30_final.eicu_trajectory_pool b

INNER JOIN public.lab l
    ON b.patientunitstayid = l.patientunitstayid

   AND l.labresultoffset >= b.sepsis_offset
   AND l.labresultoffset <  b.landmark_offset

WHERE LOWER(l.labname) IN
(
    'sodium',
    'potassium',
    'chloride',
    'bicarbonate'
)

  AND l.labresult IS NOT NULL
  AND NOT (LOWER(l.labname) = 'sodium'   AND l.labresult > 300)
  AND NOT (LOWER(l.labname) = 'chloride' AND l.labresult < 10)

GROUP BY
    b.patientunitstayid,

    FLOOR(
        (l.labresultoffset - b.sepsis_offset) / 360.0
    )::integer,

    LOWER(l.labname);


CREATE UNIQUE INDEX idx_eicu_electrolyte_6h
ON study_make30_final.eicu_electrolyte_6h
(patientunitstayid, time_bin, variable_name);



/* ========================================================================
   20. MSI 2-H TRAJECTORY LONG TABLE — 12 ROWS PER TRAJECTORY-POOL PATIENT
   ======================================================================== */
CREATE TABLE study_make30_final.eicu_msi_trajectory_long AS
WITH patient_bins AS
(
    SELECT b.uniquepid,b.patienthealthsystemstayid,b.patientunitstayid,b.sepsis_offset,b.landmark_offset,gs.time_bin
    FROM study_make30_final.eicu_trajectory_pool b CROSS JOIN generate_series(0,11) gs(time_bin)
)
SELECT p.uniquepid,p.patienthealthsystemstayid,p.patientunitstayid,p.sepsis_offset,p.landmark_offset,p.time_bin,
       ('M'||(p.time_bin+1)::text)::text AS time_label,(p.time_bin*2)::integer AS window_start_hour,((p.time_bin+1)*2)::integer AS window_end_hour,(p.time_bin*2+1)::integer AS time_midpoint_hour,
       m.msi,m.valid_hourly_msi_count
FROM patient_bins p LEFT JOIN study_make30_final.eicu_msi_2h m ON p.patientunitstayid=m.patientunitstayid AND p.time_bin=m.time_bin
ORDER BY p.patientunitstayid,p.time_bin;
CREATE UNIQUE INDEX idx_eicu_msi_trajectory_long ON study_make30_final.eicu_msi_trajectory_long(patientunitstayid,time_bin);

/* ========================================================================
   21. ELECTROLYTE 6-H TRAJECTORY LONG TABLE — 4 ROWS PER TRAJECTORY-POOL PATIENT
   ======================================================================== */
CREATE TABLE study_make30_final.eicu_electrolyte_trajectory_long AS
WITH patient_bins AS
(
    SELECT b.uniquepid,b.patienthealthsystemstayid,b.patientunitstayid,b.sepsis_offset,b.landmark_offset,gs.time_bin
    FROM study_make30_final.eicu_trajectory_pool b CROSS JOIN generate_series(0,3) gs(time_bin)
),ew AS
(
    SELECT patientunitstayid,time_bin,
           MAX(median_value) FILTER(WHERE variable_name='sodium') AS sodium,
           MAX(median_value) FILTER(WHERE variable_name='potassium') AS potassium,
           MAX(median_value) FILTER(WHERE variable_name='chloride') AS chloride,
           MAX(median_value) FILTER(WHERE variable_name='bicarbonate') AS bicarbonate
    FROM study_make30_final.eicu_electrolyte_6h GROUP BY patientunitstayid,time_bin
)
SELECT p.uniquepid,p.patienthealthsystemstayid,p.patientunitstayid,p.sepsis_offset,p.landmark_offset,p.time_bin,
       ('T'||(p.time_bin+1)::text)::text AS time_label,(p.time_bin*6)::integer AS window_start_hour,((p.time_bin+1)*6)::integer AS window_end_hour,(p.time_bin*6+3)::integer AS time_midpoint_hour,
       e.sodium,e.potassium,e.chloride,e.bicarbonate
FROM patient_bins p LEFT JOIN ew e ON p.patientunitstayid=e.patientunitstayid AND p.time_bin=e.time_bin
ORDER BY p.patientunitstayid,p.time_bin;
CREATE UNIQUE INDEX idx_eicu_electrolyte_trajectory_long ON study_make30_final.eicu_electrolyte_trajectory_long(patientunitstayid,time_bin);

/* ========================================================================
   22. PATIENT COVERAGE + FINAL PRIMARY ELIGIBILITY
   ======================================================================== */
CREATE TABLE study_make30_final.eicu_trajectory_coverage AS
WITH m AS
(
    SELECT uniquepid,patienthealthsystemstayid,patientunitstayid,COUNT(msi)::integer msi_bins FROM study_make30_final.eicu_msi_trajectory_long GROUP BY uniquepid,patienthealthsystemstayid,patientunitstayid
),e AS
(
    SELECT uniquepid,patienthealthsystemstayid,patientunitstayid,COUNT(sodium)::integer sodium_bins,COUNT(potassium)::integer potassium_bins,COUNT(chloride)::integer chloride_bins,COUNT(bicarbonate)::integer bicarbonate_bins FROM study_make30_final.eicu_electrolyte_trajectory_long GROUP BY uniquepid,patienthealthsystemstayid,patientunitstayid
)
SELECT m.*,e.sodium_bins,e.potassium_bins,e.chloride_bins,e.bicarbonate_bins,
       CASE WHEN m.msi_bins>=6 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END core4_eligible_msi6,
       CASE WHEN m.msi_bins>=9 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END core4_eligible_msi9,
       CASE WHEN m.msi_bins=12 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END core4_eligible_msi12,
       CASE WHEN m.msi_bins>=9 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END core4_trajectory_eligible
FROM m INNER JOIN e USING(uniquepid,patienthealthsystemstayid,patientunitstayid);
CREATE UNIQUE INDEX idx_eicu_trajectory_coverage ON study_make30_final.eicu_trajectory_coverage(patientunitstayid);

CREATE TABLE study_make30_final.eicu_core4_trajectory_cohort AS
SELECT c.*,v.msi_bins,v.sodium_bins,v.potassium_bins,v.chloride_bins,v.bicarbonate_bins,v.core4_eligible_msi6,v.core4_eligible_msi9,v.core4_eligible_msi12,v.core4_trajectory_eligible
FROM study_make30_final.eicu_trajectory_pool c INNER JOIN study_make30_final.eicu_trajectory_coverage v USING(patientunitstayid)
WHERE v.core4_trajectory_eligible=1;
CREATE UNIQUE INDEX idx_eicu_core4_cohort_stay ON study_make30_final.eicu_core4_trajectory_cohort(patientunitstayid);

CREATE TABLE study_make30_final.eicu_core4_msi_trajectory_long AS SELECT t.* FROM study_make30_final.eicu_msi_trajectory_long t INNER JOIN study_make30_final.eicu_trajectory_coverage c USING(patientunitstayid) WHERE c.core4_trajectory_eligible=1 ORDER BY t.patientunitstayid,t.time_bin;
CREATE UNIQUE INDEX idx_eicu_core4_msi_long ON study_make30_final.eicu_core4_msi_trajectory_long(patientunitstayid,time_bin);
CREATE TABLE study_make30_final.eicu_core4_electrolyte_trajectory_long AS SELECT t.* FROM study_make30_final.eicu_electrolyte_trajectory_long t INNER JOIN study_make30_final.eicu_trajectory_coverage c USING(patientunitstayid) WHERE c.core4_trajectory_eligible=1 ORDER BY t.patientunitstayid,t.time_bin;
CREATE UNIQUE INDEX idx_eicu_core4_electrolyte_long ON study_make30_final.eicu_core4_electrolyte_trajectory_long(patientunitstayid,time_bin);

CREATE TABLE study_make30_final.eicu_core4_multitrajectory_wide AS
WITH mw AS
(
 SELECT patientunitstayid,
 MAX(msi) FILTER(WHERE time_bin=0) msi_m1,MAX(msi) FILTER(WHERE time_bin=1) msi_m2,MAX(msi) FILTER(WHERE time_bin=2) msi_m3,MAX(msi) FILTER(WHERE time_bin=3) msi_m4,
 MAX(msi) FILTER(WHERE time_bin=4) msi_m5,MAX(msi) FILTER(WHERE time_bin=5) msi_m6,MAX(msi) FILTER(WHERE time_bin=6) msi_m7,MAX(msi) FILTER(WHERE time_bin=7) msi_m8,
 MAX(msi) FILTER(WHERE time_bin=8) msi_m9,MAX(msi) FILTER(WHERE time_bin=9) msi_m10,MAX(msi) FILTER(WHERE time_bin=10) msi_m11,MAX(msi) FILTER(WHERE time_bin=11) msi_m12
 FROM study_make30_final.eicu_core4_msi_trajectory_long GROUP BY patientunitstayid
),ew AS
(
 SELECT patientunitstayid,
 MAX(sodium) FILTER(WHERE time_bin=0) sodium_t1,MAX(sodium) FILTER(WHERE time_bin=1) sodium_t2,MAX(sodium) FILTER(WHERE time_bin=2) sodium_t3,MAX(sodium) FILTER(WHERE time_bin=3) sodium_t4,
 MAX(potassium) FILTER(WHERE time_bin=0) potassium_t1,MAX(potassium) FILTER(WHERE time_bin=1) potassium_t2,MAX(potassium) FILTER(WHERE time_bin=2) potassium_t3,MAX(potassium) FILTER(WHERE time_bin=3) potassium_t4,
 MAX(chloride) FILTER(WHERE time_bin=0) chloride_t1,MAX(chloride) FILTER(WHERE time_bin=1) chloride_t2,MAX(chloride) FILTER(WHERE time_bin=2) chloride_t3,MAX(chloride) FILTER(WHERE time_bin=3) chloride_t4,
 MAX(bicarbonate) FILTER(WHERE time_bin=0) bicarbonate_t1,MAX(bicarbonate) FILTER(WHERE time_bin=1) bicarbonate_t2,MAX(bicarbonate) FILTER(WHERE time_bin=2) bicarbonate_t3,MAX(bicarbonate) FILTER(WHERE time_bin=3) bicarbonate_t4
 FROM study_make30_final.eicu_core4_electrolyte_trajectory_long GROUP BY patientunitstayid
)
SELECT c.uniquepid,c.patienthealthsystemstayid,c.patientunitstayid,c.sepsis_offset,c.landmark_offset,
       mw.msi_m1,mw.msi_m2,mw.msi_m3,mw.msi_m4,mw.msi_m5,mw.msi_m6,mw.msi_m7,mw.msi_m8,mw.msi_m9,mw.msi_m10,mw.msi_m11,mw.msi_m12,
       ew.sodium_t1,ew.sodium_t2,ew.sodium_t3,ew.sodium_t4,
       ew.potassium_t1,ew.potassium_t2,ew.potassium_t3,ew.potassium_t4,
       ew.chloride_t1,ew.chloride_t2,ew.chloride_t3,ew.chloride_t4,
       ew.bicarbonate_t1,ew.bicarbonate_t2,ew.bicarbonate_t3,ew.bicarbonate_t4
FROM study_make30_final.eicu_core4_trajectory_cohort c LEFT JOIN mw USING(patientunitstayid) LEFT JOIN ew USING(patientunitstayid);
CREATE UNIQUE INDEX idx_eicu_core4_multitrajectory_wide ON study_make30_final.eicu_core4_multitrajectory_wide(patientunitstayid);

CREATE TABLE study_make30_final.eicu_msi_gt4_trace AS
SELECT c.uniquepid,c.patienthealthsystemstayid,h.patientunitstayid,h.hour_bin,FLOOR(h.hour_bin/2.0)::integer msi_2h_bin,('M'||(FLOOR(h.hour_bin/2.0)::integer+1)::text)::text msi_2h_label,
       h.hr_median,h.map_invasive_median,h.map_noninvasive_median,h.map_abs_diff,h.map_used,h.map_source,h.msi hourly_msi,m.msi msi_2h_max
FROM study_make30_final.eicu_msi_hourly h INNER JOIN study_make30_final.eicu_trajectory_pool c USING(patientunitstayid)
LEFT JOIN study_make30_final.eicu_msi_2h m ON h.patientunitstayid=m.patientunitstayid AND FLOOR(h.hour_bin/2.0)::integer=m.time_bin
WHERE h.msi>4 ORDER BY h.msi DESC,h.patientunitstayid,h.hour_bin;
CREATE INDEX idx_eicu_msi_gt4_trace ON study_make30_final.eicu_msi_gt4_trace(patientunitstayid,hour_bin);



/* ========================================================================
   23. ORAKLE MAKE30 + BASELINE / ADJUSTMENT COVARIATES
       FULL RAW-REBUILD IMPLEMENTATION — PERFORMANCE-OPTIMIZED ONLY

   IMPORTANT IMPLEMENTATION PRINCIPLE
   ----------------------------------
   This section does NOT reuse any previously frozen study tables from an
   older run.  The complete Sepsis/SOFA/trajectory cohort above has been
   rebuilt from raw eICU tables in this same script, exactly as in the
   validated full-cohort extraction pipeline.

   The optimization below changes only execution strategy:
     - materialize the final cohort's hospital-unit mapping once;
     - scan creatinine events once for baseline + MAKE30;
     - scan KRT event sources once;
     - scan vasopressor sources once;
     - index/analyze those small intermediate tables.

   Clinical definitions, time windows, eligibility and final variables are
   unchanged.

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

   eICU baseline creatinine implementation:
     eICU does not provide a reproducible absolute-time 12-month pre-hospital
     laboratory timeline across separate hospitalizations.  Therefore the
     historical-median component cannot be reconstructed here; MDRD
     back-calculation assuming eGFR=75 is used and the reference baseline is
     the lower of that value and the first creatinine in the index hospital
     admission before landmark.

   The frozen trajectory cohort itself is NOT changed here.
   ======================================================================== */

/* ------------------------------------------------------------------------
   23A. Small helper index for hospital-level joins
   ------------------------------------------------------------------------ */
CREATE INDEX IF NOT EXISTS idx_eicu_core4_cohort_hosp
ON study_make30_final.eicu_core4_trajectory_cohort(patienthealthsystemstayid);

ANALYZE study_make30_final.eicu_core4_trajectory_cohort;


/* ------------------------------------------------------------------------
   23B. Final cohort hospital-unit mapping

   eICU offsets are relative to each ICU unit admission.  For hospital-level
   MAKE30 we map every final index stay to all ICU units belonging to the
   same patienthealthsystemstayid, while retaining each unit's own
   hospitaladmitoffset so that event offsets can be transformed to a common
   hospital-relative time axis.
   ------------------------------------------------------------------------ */
DROP TABLE IF EXISTS tmp_eicu_final_units;

CREATE TEMP TABLE tmp_eicu_final_units AS
SELECT
    c.patientunitstayid AS index_stay_id,
    c.patienthealthsystemstayid,
    p.patientunitstayid AS unit_stay_id,
    p.hospitaladmitoffset AS unit_hospitaladmitoffset
FROM study_make30_final.eicu_core4_trajectory_cohort c
INNER JOIN public.patient p
    ON p.patienthealthsystemstayid = c.patienthealthsystemstayid;

CREATE INDEX idx_tmp_eicu_final_units_index
ON tmp_eicu_final_units(index_stay_id);

CREATE INDEX idx_tmp_eicu_final_units_unit
ON tmp_eicu_final_units(unit_stay_id);

ANALYZE tmp_eicu_final_units;


/* ------------------------------------------------------------------------
   23C. One-pass creatinine event extraction for the final cohort

   event_hosp_offset = minutes since hospital admission.
   This one materialized event table is reused for both baseline creatinine
   and the final inpatient creatinine component of MAKE30.
   ------------------------------------------------------------------------ */
DROP TABLE IF EXISTS tmp_eicu_final_creatinine;

CREATE TEMP TABLE tmp_eicu_final_creatinine AS
SELECT
    u.index_stay_id,
    u.patienthealthsystemstayid,
    u.unit_stay_id,
    l.labid,
    (l.labresultoffset - u.unit_hospitaladmitoffset)::integer AS event_hosp_offset,
    l.labresult::double precision AS creatinine
FROM tmp_eicu_final_units u
INNER JOIN public.lab l
    ON l.patientunitstayid = u.unit_stay_id
WHERE LOWER(BTRIM(l.labname)) = 'creatinine'
  AND l.labresult IS NOT NULL
  AND l.labresult BETWEEN 0.1 AND 20;

CREATE INDEX idx_tmp_eicu_final_cr_index_time
ON tmp_eicu_final_creatinine(index_stay_id,event_hosp_offset);

ANALYZE tmp_eicu_final_creatinine;


/* ------------------------------------------------------------------------
   23D. ORAKLE baseline kidney function
   ------------------------------------------------------------------------ */
CREATE TABLE study_make30_final.eicu_baseline_kidney AS
WITH first_admit AS
(
    SELECT DISTINCT ON (c.patientunitstayid)
        c.patientunitstayid,
        cr.event_hosp_offset AS first_admission_creatinine_hosp_offset,
        cr.creatinine AS first_admission_creatinine
    FROM study_make30_final.eicu_core4_trajectory_cohort c
    INNER JOIN tmp_eicu_final_creatinine cr
        ON cr.index_stay_id = c.patientunitstayid
    WHERE cr.event_hosp_offset >= 0
      AND cr.event_hosp_offset < (c.landmark_offset - c.hospitaladmitoffset)
    ORDER BY c.patientunitstayid,cr.event_hosp_offset,cr.labid
),
calc AS
(
    SELECT
        c.uniquepid,
        c.patienthealthsystemstayid,
        c.patientunitstayid,
        c.age,
        c.gender,
        c.ethnicity,
        f.first_admission_creatinine_hosp_offset,
        f.first_admission_creatinine,
        POWER(
            175.0
            * POWER(GREATEST(c.age::double precision,18.0),-0.203)
            * CASE
                WHEN LOWER(COALESCE(c.gender,'')) LIKE 'female%' THEN 0.742
                ELSE 1.0
              END
            * CASE
                WHEN LOWER(COALESCE(c.ethnicity,'')) LIKE '%african american%'
                  OR LOWER(COALESCE(c.ethnicity,'')) = 'black'
                THEN 1.212
                ELSE 1.0
              END
            / 75.0,
            1.0/1.154
        )::double precision AS mdrd_backcalc_creatinine
    FROM study_make30_final.eicu_core4_trajectory_cohort c
    LEFT JOIN first_admit f
        ON f.patientunitstayid = c.patientunitstayid
),
base AS
(
    SELECT
        x.*,
        LEAST(
            x.mdrd_backcalc_creatinine,
            COALESCE(x.first_admission_creatinine,x.mdrd_backcalc_creatinine)
        )::double precision AS baseline_creatinine
    FROM calc x
)
SELECT
    b.*,
    CASE
        WHEN b.first_admission_creatinine IS NOT NULL
         AND b.baseline_creatinine = b.first_admission_creatinine
            THEN 'first admission creatinine'
        ELSE 'MDRD back-calculation eGFR=75'
    END::text AS baseline_creatinine_source,
    CASE
        WHEN b.baseline_creatinine > 0 THEN
            175.0 * POWER(b.baseline_creatinine,-1.154)
            * POWER(GREATEST(b.age::double precision,18.0),-0.203)
            * CASE
                WHEN LOWER(COALESCE(b.gender,'')) LIKE 'female%' THEN 0.742
                ELSE 1.0
              END
            * CASE
                WHEN LOWER(COALESCE(b.ethnicity,'')) LIKE '%african american%'
                  OR LOWER(COALESCE(b.ethnicity,'')) = 'black'
                THEN 1.212
                ELSE 1.0
              END
    END::double precision AS baseline_egfr_mdrd
FROM base b;

CREATE UNIQUE INDEX idx_eicu_baseline_kidney_stay
ON study_make30_final.eicu_baseline_kidney(patientunitstayid);

ANALYZE study_make30_final.eicu_baseline_kidney;


/* ------------------------------------------------------------------------
   23E. Baseline/adjustment covariates
   ------------------------------------------------------------------------ */

/* Vasopressor exposure is materialized once rather than evaluated through
   three correlated EXISTS subqueries for every cohort row. */
DROP TABLE IF EXISTS tmp_eicu_vasopressor_0_24h;

CREATE TEMP TABLE tmp_eicu_vasopressor_0_24h AS
SELECT DISTINCT patientunitstayid
FROM
(
    SELECT c.patientunitstayid
    FROM study_make30_final.eicu_core4_trajectory_cohort c
    INNER JOIN public.infusiondrug i
        ON i.patientunitstayid = c.patientunitstayid
       AND i.infusionoffset >= c.sepsis_offset
       AND i.infusionoffset <  c.landmark_offset
    WHERE LOWER(COALESCE(i.drugname,'')) ~
          '(norepineph|noradren|epineph|adrenalin|phenylephrine|vasopressin|dopamine)'

    UNION ALL

    SELECT c.patientunitstayid
    FROM study_make30_final.eicu_core4_trajectory_cohort c
    INNER JOIN public.treatment t
        ON t.patientunitstayid = c.patientunitstayid
       AND t.treatmentoffset >= c.sepsis_offset
       AND t.treatmentoffset <  c.landmark_offset
    WHERE LOWER(COALESCE(t.treatmentstring,'')) ~
          '(norepineph|noradren|epineph|adrenalin|phenylephrine|vasopressin|dopamine)'

    UNION ALL

    SELECT c.patientunitstayid
    FROM study_make30_final.eicu_core4_trajectory_cohort c
    INNER JOIN public.medication m
        ON m.patientunitstayid = c.patientunitstayid
       AND m.drugstartoffset < c.landmark_offset
       AND COALESCE(m.drugstopoffset,m.drugstartoffset+1) > c.sepsis_offset
    WHERE COALESCE(m.drugordercancelled,'No') = 'No'
      AND LOWER(COALESCE(m.drugname,'')) ~
          '(norepineph|noradren|epineph|adrenalin|phenylephrine|vasopressin|dopamine)'
) s;

CREATE UNIQUE INDEX idx_tmp_eicu_vasopressor_stay
ON tmp_eicu_vasopressor_0_24h(patientunitstayid);

ANALYZE tmp_eicu_vasopressor_0_24h;


CREATE TABLE study_make30_final.eicu_covariates AS
WITH ph AS
(
    SELECT
        c.patientunitstayid,
        MAX(
            CASE
                WHEN LOWER(
                    COALESCE(p.pasthistoryvaluetext,'') || ' ' ||
                    COALESCE(p.pasthistoryvalue,'') || ' ' ||
                    COALESCE(p.pasthistorypath,'')
                ) ~ '(chf|congestive heart failure|heart failure)'
                THEN 1 ELSE 0
            END
        )::integer AS congestive_heart_failure,
        MAX(
            CASE
                WHEN LOWER(
                    COALESCE(p.pasthistoryvaluetext,'') || ' ' ||
                    COALESCE(p.pasthistoryvalue,'') || ' ' ||
                    COALESCE(p.pasthistorypath,'')
                ) ~ '(chronic renal failure|chronic kidney|renal insufficien|renal failure)'
                THEN 1 ELSE 0
            END
        )::integer AS ckd_history,
        MAX(
            CASE
                WHEN LOWER(
                    COALESCE(p.pasthistoryvaluetext,'') || ' ' ||
                    COALESCE(p.pasthistoryvalue,'') || ' ' ||
                    COALESCE(p.pasthistorypath,'')
                ) ~ '(diabet)'
                THEN 1 ELSE 0
            END
        )::integer AS diabetes_history
    FROM study_make30_final.eicu_core4_trajectory_cohort c
    LEFT JOIN public.pasthistory p
        ON p.patientunitstayid = c.patientunitstayid
    GROUP BY c.patientunitstayid
),
vent AS
(
    SELECT DISTINCT
        c.patientunitstayid,
        1::integer AS mechanical_ventilation_0_24h
    FROM study_make30_final.eicu_core4_trajectory_cohort c
    INNER JOIN tmp_eicu_vent_interval v
        ON v.patientunitstayid = c.patientunitstayid
       AND v.startoffset < c.landmark_offset
       AND v.endoffset   > c.sepsis_offset
),
ap AS
(
    SELECT DISTINCT ON (a.patientunitstayid)
        a.patientunitstayid,
        a.diabetes,
        a.hepaticfailure,
        a.cirrhosis,
        a.lymphoma,
        a.metastaticcancer,
        a.leukemia,
        a.immunosuppression,
        a.aids,
        a.ventday1,
        a.oobventday1,
        a.oobintubday1
    FROM public.apachepredvar a
    INNER JOIN study_make30_final.eicu_core4_trajectory_cohort c
        ON c.patientunitstayid = a.patientunitstayid
    ORDER BY a.patientunitstayid
)
SELECT
    c.uniquepid,
    c.patienthealthsystemstayid,
    c.patientunitstayid,
    c.age::double precision AS age,
    c.gender AS sex,
    c.ethnicity AS race_ethnicity,
    CASE
        WHEN p.admissionweight BETWEEN 25 AND 300
            THEN p.admissionweight::double precision
    END AS weight_kg,
    CASE
        WHEN p.admissionheight BETWEEN 100 AND 250
            THEN p.admissionheight::double precision
    END AS height_cm,
    CASE
        WHEN p.admissionweight BETWEEN 25 AND 300
         AND p.admissionheight BETWEEN 100 AND 250
            THEN p.admissionweight / POWER(p.admissionheight/100.0,2)
    END::double precision AS bmi,
    c.sofa_score::double precision AS sepsis_sofa_score,
    bk.first_admission_creatinine,
    bk.mdrd_backcalc_creatinine,
    bk.baseline_creatinine_source,
    bk.baseline_creatinine,
    bk.baseline_egfr_mdrd,
    COALESCE(ph.ckd_history,0)::integer AS ckd_or_renal_disease,
    GREATEST(COALESCE(ph.diabetes_history,0),COALESCE(ap.diabetes,0))::integer AS diabetes,
    COALESCE(ph.congestive_heart_failure,0)::integer AS congestive_heart_failure,
    GREATEST(COALESCE(ap.hepaticfailure,0),COALESCE(ap.cirrhosis,0))::integer AS liver_disease,
    GREATEST(
        COALESCE(ap.lymphoma,0),
        COALESCE(ap.metastaticcancer,0),
        COALESCE(ap.leukemia,0)
    )::integer AS malignancy,
    COALESCE(ap.immunosuppression,0)::integer AS immunosuppression,
    COALESCE(v.mechanical_ventilation_0_24h,0)::integer AS mechanical_ventilation_0_24h,
    CASE WHEN x.patientunitstayid IS NOT NULL THEN 1 ELSE 0 END::integer AS vasopressor_0_24h
FROM study_make30_final.eicu_core4_trajectory_cohort c
LEFT JOIN public.patient p
    ON p.patientunitstayid = c.patientunitstayid
LEFT JOIN study_make30_final.eicu_baseline_kidney bk
    ON bk.patientunitstayid = c.patientunitstayid
LEFT JOIN ph
    ON ph.patientunitstayid = c.patientunitstayid
LEFT JOIN ap
    ON ap.patientunitstayid = c.patientunitstayid
LEFT JOIN vent v
    ON v.patientunitstayid = c.patientunitstayid
LEFT JOIN tmp_eicu_vasopressor_0_24h x
    ON x.patientunitstayid = c.patientunitstayid;

CREATE UNIQUE INDEX idx_eicu_covariates_stay
ON study_make30_final.eicu_covariates(patientunitstayid);

ANALYZE study_make30_final.eicu_covariates;


/* ------------------------------------------------------------------------
   23F. MAKE30 follow-up frame
   ------------------------------------------------------------------------ */
DROP TABLE IF EXISTS tmp_eicu_make30_fup;

CREATE TEMP TABLE tmp_eicu_make30_fup AS
SELECT
    c.uniquepid,
    c.patienthealthsystemstayid,
    c.patientunitstayid,
    c.hospitaladmitoffset,
    c.hospitaldischargeoffset,
    c.hospitaldischargestatus,
    (c.landmark_offset - c.hospitaladmitoffset)::integer AS landmark_hosp_offset,
    LEAST(
        COALESCE(
            c.hospitaldischargeoffset - c.hospitaladmitoffset,
            -c.hospitaladmitoffset + 43200
        ),
        -c.hospitaladmitoffset + 43200
    )::integer AS make30_censor_hosp_offset
FROM study_make30_final.eicu_core4_trajectory_cohort c;

CREATE UNIQUE INDEX idx_tmp_eicu_make30_fup_stay
ON tmp_eicu_make30_fup(patientunitstayid);

ANALYZE tmp_eicu_make30_fup;


/* ------------------------------------------------------------------------
   23G. One-pass KRT event extraction for all hospital units in final cohort
   ------------------------------------------------------------------------ */
DROP TABLE IF EXISTS tmp_eicu_final_krt_events;

CREATE TEMP TABLE tmp_eicu_final_krt_events AS
SELECT
    u.index_stay_id,
    (t.treatmentoffset - u.unit_hospitaladmitoffset)::integer AS event_hosp_offset,
    'treatment'::text AS krt_source
FROM tmp_eicu_final_units u
INNER JOIN public.treatment t
    ON t.patientunitstayid = u.unit_stay_id
WHERE LOWER(COALESCE(t.treatmentstring,'')) ~
      '(dialysis|hemodialysis|haemodialysis|crrt|cvvh|cvvhd|cvvhdf|sled|peritoneal)'

UNION ALL

SELECT
    u.index_stay_id,
    (io.intakeoutputoffset - u.unit_hospitaladmitoffset)::integer AS event_hosp_offset,
    'intakeoutput'::text AS krt_source
FROM tmp_eicu_final_units u
INNER JOIN public.intakeoutput io
    ON io.patientunitstayid = u.unit_stay_id
WHERE COALESCE(io.dialysistotal,0) <> 0
   OR LOWER(COALESCE(io.celllabel,'')) LIKE '%dialysis%'
   OR LOWER(COALESCE(io.cellpath,''))  LIKE '%dialysis%';

CREATE INDEX idx_tmp_eicu_final_krt_index_time
ON tmp_eicu_final_krt_events(index_stay_id,event_hosp_offset);

ANALYZE tmp_eicu_final_krt_events;


/* ------------------------------------------------------------------------
   23H. ORAKLE MAKE30
   ------------------------------------------------------------------------ */
CREATE TABLE study_make30_final.eicu_make30 AS
WITH prior_krt AS
(
    SELECT
        f.patientunitstayid,
        1::integer AS prior_krt_before_landmark
    FROM tmp_eicu_make30_fup f
    INNER JOIN tmp_eicu_final_krt_events e
        ON e.index_stay_id = f.patientunitstayid
       AND e.event_hosp_offset <= f.landmark_hosp_offset
    GROUP BY f.patientunitstayid
),
krt_events AS
(
    SELECT
        f.patientunitstayid,
        MIN(e.event_hosp_offset)::integer AS first_new_krt_hosp_offset
    FROM tmp_eicu_make30_fup f
    LEFT JOIN prior_krt pk
        ON pk.patientunitstayid = f.patientunitstayid
    INNER JOIN tmp_eicu_final_krt_events e
        ON e.index_stay_id = f.patientunitstayid
       AND e.event_hosp_offset >  f.landmark_hosp_offset
       AND e.event_hosp_offset <= f.make30_censor_hosp_offset
    WHERE pk.patientunitstayid IS NULL
    GROUP BY f.patientunitstayid
),
final_cr AS
(
    SELECT DISTINCT ON (f.patientunitstayid)
        f.patientunitstayid,
        cr.event_hosp_offset AS final_inpatient_creatinine_hosp_offset,
        cr.creatinine AS final_inpatient_creatinine
    FROM tmp_eicu_make30_fup f
    INNER JOIN tmp_eicu_final_creatinine cr
        ON cr.index_stay_id = f.patientunitstayid
       AND cr.event_hosp_offset >  f.landmark_hosp_offset
       AND cr.event_hosp_offset <= f.make30_censor_hosp_offset
    ORDER BY f.patientunitstayid,cr.event_hosp_offset DESC,cr.labid DESC
)
SELECT
    f.uniquepid,
    f.patienthealthsystemstayid,
    f.patientunitstayid,
    f.landmark_hosp_offset AS make30_risk_start_hosp_offset,
    f.make30_censor_hosp_offset,
    (f.make30_censor_hosp_offset - f.landmark_hosp_offset)/1440.0
        AS make30_followup_days_after_landmark,
    CASE
        WHEN f.make30_censor_hosp_offset > f.landmark_hosp_offset THEN 1
        ELSE 0
    END::integer AS make30_evaluable_window,
    CASE
        WHEN LOWER(COALESCE(f.hospitaldischargestatus,'')) = 'expired'
         AND (f.hospitaldischargeoffset - f.hospitaladmitoffset) > f.landmark_hosp_offset
         AND (f.hospitaldischargeoffset - f.hospitaladmitoffset) <= f.make30_censor_hosp_offset
            THEN 1
        ELSE 0
    END::integer AS make30_death,
    COALESCE(pk.prior_krt_before_landmark,0)::integer AS prior_krt_before_landmark,
    CASE
        WHEN k.first_new_krt_hosp_offset IS NOT NULL THEN 1
        ELSE 0
    END::integer AS make30_new_krt,
    k.first_new_krt_hosp_offset,
    fc.final_inpatient_creatinine,
    fc.final_inpatient_creatinine_hosp_offset,
    bk.baseline_creatinine,
    CASE
        WHEN fc.final_inpatient_creatinine IS NOT NULL
         AND bk.baseline_creatinine IS NOT NULL
         AND fc.final_inpatient_creatinine >= 2.0*bk.baseline_creatinine
            THEN 1
        WHEN fc.final_inpatient_creatinine IS NOT NULL
         AND bk.baseline_creatinine IS NOT NULL
            THEN 0
        ELSE NULL
    END::integer AS make30_persistent_kidney_dysfunction,
    CASE
        WHEN f.make30_censor_hosp_offset <= f.landmark_hosp_offset
            THEN NULL
        WHEN LOWER(COALESCE(f.hospitaldischargestatus,'')) = 'expired'
         AND (f.hospitaldischargeoffset - f.hospitaladmitoffset) > f.landmark_hosp_offset
         AND (f.hospitaldischargeoffset - f.hospitaladmitoffset) <= f.make30_censor_hosp_offset
            THEN 1
        WHEN k.first_new_krt_hosp_offset IS NOT NULL
            THEN 1
        WHEN fc.final_inpatient_creatinine IS NOT NULL
         AND bk.baseline_creatinine IS NOT NULL
         AND fc.final_inpatient_creatinine >= 2.0*bk.baseline_creatinine
            THEN 1
        WHEN fc.final_inpatient_creatinine IS NOT NULL
         AND bk.baseline_creatinine IS NOT NULL
            THEN 0
        ELSE NULL
    END::integer AS make30
FROM tmp_eicu_make30_fup f
LEFT JOIN prior_krt pk
    ON pk.patientunitstayid = f.patientunitstayid
LEFT JOIN krt_events k
    ON k.patientunitstayid = f.patientunitstayid
LEFT JOIN final_cr fc
    ON fc.patientunitstayid = f.patientunitstayid
LEFT JOIN study_make30_final.eicu_baseline_kidney bk
    ON bk.patientunitstayid = f.patientunitstayid;

CREATE UNIQUE INDEX idx_eicu_make30_stay
ON study_make30_final.eicu_make30(patientunitstayid);

ANALYZE study_make30_final.eicu_make30;


/* ------------------------------------------------------------------------
   23I. Final one-row-per-patient analysis table
   ------------------------------------------------------------------------ */
CREATE TABLE study_make30_final.eicu_analysis_master AS
SELECT
    w.*,
    cv.age,
    cv.sex,
    cv.race_ethnicity,
    cv.weight_kg,
    cv.height_cm,
    cv.bmi,
    cv.sepsis_sofa_score,
    cv.first_admission_creatinine,
    cv.mdrd_backcalc_creatinine,
    cv.baseline_creatinine_source,
    cv.baseline_creatinine,
    cv.baseline_egfr_mdrd,
    cv.ckd_or_renal_disease,
    cv.diabetes,
    cv.congestive_heart_failure,
    cv.liver_disease,
    cv.malignancy,
    cv.immunosuppression,
    cv.mechanical_ventilation_0_24h,
    cv.vasopressor_0_24h,
    o.make30_risk_start_hosp_offset,
    o.make30_censor_hosp_offset,
    o.make30_followup_days_after_landmark,
    o.make30_evaluable_window,
    o.make30_death,
    o.prior_krt_before_landmark,
    o.make30_new_krt,
    o.first_new_krt_hosp_offset,
    o.final_inpatient_creatinine,
    o.final_inpatient_creatinine_hosp_offset,
    o.make30_persistent_kidney_dysfunction,
    o.make30
FROM study_make30_final.eicu_core4_multitrajectory_wide w
LEFT JOIN study_make30_final.eicu_covariates cv
    ON cv.patientunitstayid = w.patientunitstayid
LEFT JOIN study_make30_final.eicu_make30 o
    ON o.patientunitstayid = w.patientunitstayid;

CREATE UNIQUE INDEX idx_eicu_analysis_master_stay
ON study_make30_final.eicu_analysis_master(patientunitstayid);

ANALYZE study_make30_final.eicu_analysis_master;


/* ========================================================================
   24. FINAL QC + ANALYSIS RESULT
   ======================================================================== */
CREATE TABLE study_make30_final.eicu_qc_summary AS
WITH implementation_info AS
(
 SELECT '00_IMPLEMENTATION_INFO'::text section,0 sort1,1 sort2,'FINAL analysis pipeline: trajectories + landmark-start MAKE30 + covariates; SOFA = standard original SOFA reconstructed directly from raw eICU tables; mixed-frequency trajectory MSI 2-h x12 + electrolytes 6-h x4; primary MSI completeness >=9/12'::text item,NULL::text time_bin,NULL::bigint n_total,NULL::bigint n_observed,NULL::bigint n_missing,NULL::numeric observed_pct,NULL::numeric missing_pct,NULL::bigint n_meet_requirement,NULL::bigint n_all_bins,NULL::numeric mean_value,NULL::numeric sd_value,NULL::numeric min_value,NULL::numeric p01,NULL::numeric p05,NULL::numeric p25,NULL::numeric median_value,NULL::numeric p75,NULL::numeric p95,NULL::numeric p99,NULL::numeric max_value,NULL::text detail
),
cohort_flow AS
(
 SELECT '01_COHORT_FLOW'::text section,1 sort1,1 sort2,'all unique ICU patients'::text item,NULL::text time_bin,COUNT(*)::bigint n_total,COUNT(*)::bigint n_observed,NULL::bigint n_missing,100::numeric observed_pct,0::numeric missing_pct,NULL::bigint n_meet_requirement,NULL::bigint n_all_bins,NULL::numeric mean_value,NULL::numeric sd_value,NULL::numeric min_value,NULL::numeric p01,NULL::numeric p05,NULL::numeric p25,NULL::numeric median_value,NULL::numeric p75,NULL::numeric p95,NULL::numeric p99,NULL::numeric max_value,NULL::text detail FROM study_make30_final.eicu_patients
 UNION ALL SELECT '01_COHORT_FLOW',1,2,'adult patients with ORAKLE-style Sepsis onset inside ICU',NULL,(SELECT COUNT(*) FROM study_make30_final.eicu_patients),COUNT(DISTINCT uniquepid),NULL,ROUND(100.0*COUNT(DISTINCT uniquepid)/NULLIF((SELECT COUNT(*) FROM study_make30_final.eicu_patients),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.eicu_sepsis_adult
 UNION ALL SELECT '01_COHORT_FLOW',1,3,'first qualifying Sepsis ICU stay per patient',NULL,(SELECT COUNT(*) FROM study_make30_final.eicu_patients),COUNT(*),NULL,ROUND(100.0*COUNT(*)/NULLIF((SELECT COUNT(*) FROM study_make30_final.eicu_patients),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.eicu_sepsis_index
 UNION ALL SELECT '01_COHORT_FLOW',1,4,'reached Sepsis +24 h landmark in ICU/hospital',NULL,COUNT(*),COUNT(*) FILTER(WHERE reached_24h_landmark=1),COUNT(*) FILTER(WHERE reached_24h_landmark=0),ROUND(100.0*COUNT(*) FILTER(WHERE reached_24h_landmark=1)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE reached_24h_landmark=0)/NULLIF(COUNT(*),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.eicu_landmark_cohort
 UNION ALL SELECT '01_COHORT_FLOW',1,5,'trajectory pool after landmark + early-KRT exclusion',NULL,(SELECT COUNT(*) FROM study_make30_final.eicu_landmark_cohort WHERE reached_24h_landmark=1),COUNT(*),NULL,ROUND(100.0*COUNT(*)/NULLIF((SELECT COUNT(*) FROM study_make30_final.eicu_landmark_cohort WHERE reached_24h_landmark=1),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.eicu_trajectory_pool
 UNION ALL SELECT '01_COHORT_FLOW',1,6,'PRIMARY mixed-frequency Core-4 eligible: MSI>=9/12 + Na/K/Cl/HCO3 each>=2/4',NULL,COUNT(*),COUNT(*) FILTER(WHERE core4_trajectory_eligible=1),COUNT(*) FILTER(WHERE core4_trajectory_eligible=0),ROUND(100.0*COUNT(*) FILTER(WHERE core4_trajectory_eligible=1)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE core4_trajectory_eligible=0)/NULLIF(COUNT(*),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.eicu_trajectory_coverage
),
variable_list AS (SELECT * FROM (VALUES (1,'MSI'::text,9,12),(2,'sodium'::text,2,4),(3,'potassium'::text,2,4),(4,'chloride'::text,2,4),(5,'bicarbonate'::text,2,4)) v(variable_order,variable_name,required_bins,max_bins)),
patient_coverage AS
(
 SELECT '02_PATIENT_LEVEL_COVERAGE'::text section,2 sort1,v.variable_order sort2,v.variable_name item,NULL::text time_bin,COUNT(*)::bigint n_total,
 COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END>=1)::bigint n_observed,
 COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END=0)::bigint n_missing,
 ROUND(100.0*COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END>=1)/NULLIF(COUNT(*),0),2),
 ROUND(100.0*COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END=0)/NULLIF(COUNT(*),0),2),
 COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END>=v.required_bins)::bigint,
 COUNT(*) FILTER(WHERE CASE v.variable_name WHEN 'MSI' THEN c.msi_bins WHEN 'sodium' THEN c.sodium_bins WHEN 'potassium' THEN c.potassium_bins WHEN 'chloride' THEN c.chloride_bins WHEN 'bicarbonate' THEN c.bicarbonate_bins END=v.max_bins)::bigint,
 NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,('required_bins='||v.required_bins||'; max_bins='||v.max_bins)::text
 FROM study_make30_final.eicu_trajectory_coverage c CROSS JOIN variable_list v GROUP BY v.variable_order,v.variable_name,v.required_bins,v.max_bins
),
all_series AS
(
 SELECT patientunitstayid,1 variable_order,'MSI'::text variable_name,time_bin,time_label,msi::double precision value_num FROM study_make30_final.eicu_msi_trajectory_long
 UNION ALL SELECT patientunitstayid,2,'sodium',time_bin,time_label,sodium::double precision FROM study_make30_final.eicu_electrolyte_trajectory_long
 UNION ALL SELECT patientunitstayid,3,'potassium',time_bin,time_label,potassium::double precision FROM study_make30_final.eicu_electrolyte_trajectory_long
 UNION ALL SELECT patientunitstayid,4,'chloride',time_bin,time_label,chloride::double precision FROM study_make30_final.eicu_electrolyte_trajectory_long
 UNION ALL SELECT patientunitstayid,5,'bicarbonate',time_bin,time_label,bicarbonate::double precision FROM study_make30_final.eicu_electrolyte_trajectory_long
),
eligible_series AS
(
 SELECT patientunitstayid,1 variable_order,'MSI'::text variable_name,time_bin,time_label,msi::double precision value_num FROM study_make30_final.eicu_core4_msi_trajectory_long
 UNION ALL SELECT patientunitstayid,2,'sodium',time_bin,time_label,sodium::double precision FROM study_make30_final.eicu_core4_electrolyte_trajectory_long
 UNION ALL SELECT patientunitstayid,3,'potassium',time_bin,time_label,potassium::double precision FROM study_make30_final.eicu_core4_electrolyte_trajectory_long
 UNION ALL SELECT patientunitstayid,4,'chloride',time_bin,time_label,chloride::double precision FROM study_make30_final.eicu_core4_electrolyte_trajectory_long
 UNION ALL SELECT patientunitstayid,5,'bicarbonate',time_bin,time_label,bicarbonate::double precision FROM study_make30_final.eicu_core4_electrolyte_trajectory_long
),
all_missing AS (SELECT '03_BIN_LEVEL_MISSINGNESS'::text section,3 sort1,(variable_order*100+time_bin)::integer sort2,variable_name item,time_label time_bin,COUNT(*)::bigint n_total,COUNT(value_num)::bigint n_observed,COUNT(*) FILTER(WHERE value_num IS NULL)::bigint n_missing,ROUND(100.0*COUNT(value_num)/NULLIF(COUNT(*),0),2) observed_pct,ROUND(100.0*COUNT(*) FILTER(WHERE value_num IS NULL)/NULLIF(COUNT(*),0),2) missing_pct,NULL::bigint n_meet_requirement,NULL::bigint n_all_bins,NULL::numeric mean_value,NULL::numeric sd_value,NULL::numeric min_value,NULL::numeric p01,NULL::numeric p05,NULL::numeric p25,NULL::numeric median_value,NULL::numeric p75,NULL::numeric p95,NULL::numeric p99,NULL::numeric max_value,NULL::text detail FROM all_series GROUP BY variable_order,variable_name,time_bin,time_label),
elig_missing AS (SELECT '04_CORE4_ELIGIBLE_BIN_MISSINGNESS'::text section,4 sort1,(variable_order*100+time_bin)::integer sort2,variable_name item,time_label time_bin,COUNT(*)::bigint,COUNT(value_num)::bigint,COUNT(*) FILTER(WHERE value_num IS NULL)::bigint,ROUND(100.0*COUNT(value_num)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE value_num IS NULL)/NULLIF(COUNT(*),0),2),NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::text FROM eligible_series GROUP BY variable_order,variable_name,time_bin,time_label),
value_qc AS (SELECT '05_CORE4_VALUE_DISTRIBUTION_QC'::text section,5 sort1,(variable_order*100+time_bin)::integer sort2,variable_name item,time_label time_bin,COUNT(*)::bigint,COUNT(value_num)::bigint,COUNT(*) FILTER(WHERE value_num IS NULL)::bigint,ROUND(100.0*COUNT(value_num)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE value_num IS NULL)/NULLIF(COUNT(*),0),2),NULL::bigint,NULL::bigint,AVG(value_num)::numeric,STDDEV_SAMP(value_num)::numeric,MIN(value_num)::numeric,PERCENTILE_CONT(.01) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.05) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.25) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.5) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.75) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.95) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.99) WITHIN GROUP(ORDER BY value_num)::numeric,MAX(value_num)::numeric,NULL::text FROM eligible_series GROUP BY variable_order,variable_name,time_bin,time_label),
threshold_qc AS
(
 SELECT '06_MSI_COVERAGE_THRESHOLD_SENSITIVITY'::text section,6 sort1,x.ord sort2,x.label item,NULL::text time_bin,COUNT(*)::bigint,COUNT(*) FILTER(WHERE x.flag=1)::bigint,COUNT(*) FILTER(WHERE x.flag=0)::bigint,ROUND(100.0*COUNT(*) FILTER(WHERE x.flag=1)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE x.flag=0)/NULLIF(COUNT(*),0),2),NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::text
 FROM study_make30_final.eicu_trajectory_coverage c CROSS JOIN LATERAL (VALUES (1,'Core-4 if MSI>=6/12',c.core4_eligible_msi6),(2,'Core-4 if MSI>=9/12 [PRIMARY]',c.core4_eligible_msi9),(3,'Core-4 if MSI=12/12',c.core4_eligible_msi12)) x(ord,label,flag) GROUP BY x.ord,x.label
),
anomaly_trace AS
(
 SELECT '07_MSI_GT4_TRACE'::text section,7 sort1,ROW_NUMBER() OVER(ORDER BY hourly_msi DESC,patientunitstayid,hour_bin)::integer sort2,'hourly MSI > 4'::text item,msi_2h_label::text time_bin,1::bigint,1::bigint,0::bigint,100::numeric,0::numeric,NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,
 jsonb_build_object('uniquepid',uniquepid,'patienthealthsystemstayid',patienthealthsystemstayid,'patientunitstayid',patientunitstayid,'hour_bin',hour_bin,'msi_2h_bin',msi_2h_bin,'hr_median',hr_median,'map_invasive_median',map_invasive_median,'map_noninvasive_median',map_noninvasive_median,'map_abs_diff',map_abs_diff,'map_used',map_used,'map_source',map_source,'hourly_msi',hourly_msi,'msi_2h_max',msi_2h_max)::text FROM study_make30_final.eicu_msi_gt4_trace
),
structure_qc AS
(
 SELECT '08_STRUCTURE_QC'::text section,8 sort1,1 sort2,'eligible patients'::text item,NULL::text time_bin,COUNT(*)::bigint,NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::text FROM study_make30_final.eicu_core4_trajectory_cohort
 UNION ALL SELECT '08_STRUCTURE_QC',8,2,'MSI long rows (expected eligible N x 12)',NULL,COUNT(*),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.eicu_core4_msi_trajectory_long
 UNION ALL SELECT '08_STRUCTURE_QC',8,3,'electrolyte long rows (expected eligible N x 4)',NULL,COUNT(*),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.eicu_core4_electrolyte_trajectory_long
)
SELECT section,item,time_bin,n_total,n_observed,n_missing,observed_pct,missing_pct,n_meet_requirement,n_all_bins,ROUND(mean_value,4) mean_value,ROUND(sd_value,4) sd_value,ROUND(min_value,4) min_value,ROUND(p01,4) p01,ROUND(p05,4) p05,ROUND(p25,4) p25,ROUND(median_value,4) median_value,ROUND(p75,4) p75,ROUND(p95,4) p95,ROUND(p99,4) p99,ROUND(max_value,4) max_value,detail
FROM (SELECT * FROM implementation_info UNION ALL SELECT * FROM cohort_flow UNION ALL SELECT * FROM patient_coverage UNION ALL SELECT * FROM all_missing UNION ALL SELECT * FROM elig_missing UNION ALL SELECT * FROM value_qc UNION ALL SELECT * FROM threshold_qc UNION ALL SELECT * FROM anomaly_trace UNION ALL SELECT * FROM structure_qc) z ORDER BY sort1,sort2;


INSERT INTO study_make30_final.eicu_qc_summary
(section,item,time_bin,n_total,n_observed,n_missing,observed_pct,missing_pct,n_meet_requirement,n_all_bins,mean_value,sd_value,min_value,p01,p05,p25,median_value,p75,p95,p99,max_value,detail)
SELECT '09_MAKE30_QC','MAKE30 evaluable / events',NULL,COUNT(*)::bigint,COUNT(make30)::bigint,COUNT(*) FILTER(WHERE make30 IS NULL)::bigint,
       ROUND(100.0*COUNT(make30)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE make30 IS NULL)/NULLIF(COUNT(*),0),2),COUNT(*) FILTER(WHERE make30=1)::bigint,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
       jsonb_build_object('death',COUNT(*) FILTER(WHERE make30_death=1),'new_krt',COUNT(*) FILTER(WHERE make30_new_krt=1),'persistent_kidney_dysfunction',COUNT(*) FILTER(WHERE make30_persistent_kidney_dysfunction=1),'missing_final_creatinine',COUNT(*) FILTER(WHERE final_inpatient_creatinine IS NULL))::text
FROM study_make30_final.eicu_make30;

INSERT INTO study_make30_final.eicu_qc_summary
(section,item,time_bin,n_total,n_observed,n_missing,observed_pct,missing_pct,n_meet_requirement,n_all_bins,mean_value,sd_value,min_value,p01,p05,p25,median_value,p75,p95,p99,max_value,detail)
SELECT '10_COVARIATE_QC',v.var,NULL,COUNT(*)::bigint,COUNT(v.val)::bigint,(COUNT(*)-COUNT(v.val))::bigint,
       ROUND(100.0*COUNT(v.val)/NULLIF(COUNT(*),0),2),ROUND(100.0*(COUNT(*)-COUNT(v.val))/NULLIF(COUNT(*),0),2),NULL,NULL,
       AVG(v.val),STDDEV_SAMP(v.val),MIN(v.val),PERCENTILE_CONT(.01) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.05) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.25) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.5) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.75) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.95) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.99) WITHIN GROUP(ORDER BY v.val),MAX(v.val),NULL
FROM study_make30_final.eicu_covariates c CROSS JOIN LATERAL (VALUES ('age',c.age),('bmi',c.bmi),('sepsis_sofa_score',c.sepsis_sofa_score),('baseline_creatinine',c.baseline_creatinine),('baseline_egfr_mdrd',c.baseline_egfr_mdrd)) v(var,val)
GROUP BY v.var;

SELECT * FROM study_make30_final.eicu_analysis_master ORDER BY uniquepid,patientunitstayid;
