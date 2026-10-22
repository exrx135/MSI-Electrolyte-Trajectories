/* FINAL ANALYSIS EXTENSION: preserves frozen trajectory extraction and adds ORAKLE MAKE30 + covariates.
   Outcome events are counted only after Sepsis t0+24h landmark.
   Generated 2026-08-27. */

/* =========================================================================
   MSI STUDY — SICdb v1.0.8
   FULL-COHORT FORMAL RESEARCH PIPELINE — MIXED-FREQUENCY CORE-4 VERSION: MSI 2-H + ELECTROLYTES 6-H
   PostgreSQL — directly runnable from the first line

   PRIMARY TRAJECTORY
   ------------------
       MSI + Na + K + Cl + HCO3

   PURPOSE
   -------
   Run the same formal research workflow used for MIMIC-IV/eICU on
   the full SICdb patient cohort:

      1) full cohort at PatientID level
      2) retain all CaseID admissions of all patients
      3) ORAKLE-style suspected infection:
           occurrence of a recognized antibiotic administration
      4) pair with SOFA >=2:
           SOFA time from -24 h through +12 h relative to antibiotic
      5) Sepsis onset t0:
           earlier of antibiotic anchor and qualifying SOFA time
      6) adult + t0 occurring after SICdb ICUOffset
      7) first qualifying Sepsis case per patient
      8) require observation through t0 +24 h
      9) exclude KRT/CRRT activity during t0 to t0+24 h
     10) extract 24 hourly MSI values
     11) summarize MSI into twelve 2-h bins
     12) extract Na/K/Cl/HCO3 into four 6-h bins
     13) require:
           MSI >=9 distinct 2-h bins
           Na  >=2
           K   >=2
           Cl  >=2
           HCO3>=2
     14) create Core-4 GBMTM long table
     15) output:
           implementation/source mapping
           cohort flow
           patient-level coverage
           all-candidate bin missingness
           Core-4-eligible bin missingness

   -------------------------------------------------------------------------
   IMPORTANT DATABASE-SPECIFIC FACTS
   -------------------------------------------------------------------------
   LOCAL POSTGRESQL COLUMN-NAME COMPATIBILITY
   ------------------------------------------
   The source-locking block matches SICdb logical fields after removing
   case differences and separators. Thus all of these are accepted:
       CaseID / caseid / case_id
       PatientID / patientid / patient_id
       TimeOfStay / timeofstay / time_of_stay
       ICUOffset / icuoffset / icu_offset
       AgeOnAdmission / ageonadmission / age_on_admission
   The same normalization is applied to all source-table column mappings.

   Semantic name matching is used for all CRITICAL TIME fields.
   Ordinal fallback is NEVER allowed for the time offset of data_float_h,
   medication, laboratory or data_range. If a critical time field cannot be
   identified, the script stops rather than guessing.

   Non-time identifiers/optional fields may still use documented fallback
   logic where necessary. The selected physical column names and the mapping
   method are emitted in IMPLEMENTATION_INFO.

   Your local SICdb import uses names such as offset_value; these are now
   explicitly recognized as time offsets.

   SICdb time offsets are in SECONDS from the primary MetaVision admission.

   cases.ICUOffset:
       seconds from primary admission to first actual ICU/intermediate-care
       bed assignment. SICdb may contain preceding surgery.

   data_float_h:
       one-hour aggregated signals.
       The public field "Val" is the SICdb-provided hourly aggregate
       (arithmetic mean for quantitative monitor signals).

   Therefore the SICdb MSI implementation is:

       hourly HR  = database-provided hourly aggregate
       hourly MAP = database-provided hourly aggregate
       hourly MSI = hourly HR / hourly MAP
       2-h MSI    = maximum valid hourly MSI in each 2-h interval

   This is a DATABASE-NATIVE hourly implementation. The minute-level source
   is serialized in data_float_h.rawdata and is not required for this
   feasibility/transport pipeline.

   MAP QC / source selection:
       both invasive MAP and NIBP <=20 mmHg are excluded before source selection;
       if both valid hourly MAP sources differ by >=30 mmHg, use NIBP;
       otherwise valid invasive MAP has priority

   HR priority:
       ECG HR > SpO2 HR > arterial-pressure HR

   -------------------------------------------------------------------------
   ORAKLE-STYLE SEPSIS
   -------------------------------------------------------------------------
   Published ORAKLE rule:
       suspected infection + SOFA increase >=2

   Baseline SOFA before ICU:
       assumed 0

   SICdb/eICU suspected infection:
       occurrence of antibiotic administrations because fluid-sampling data
       are limited.

   Qualifying SOFA:
       no more than 24 h before and no more than 12 h after antibiotic time.

   t0:
       earlier of suspected-infection time and qualifying SOFA time.

   IMPORTANT:
       ORAKLE author-side source code/raw antibiotic dictionary is not public.
       This script therefore uses a transparent auditable antibiotic name
       dictionary against SICdb d_references.ReferenceValue.

   -------------------------------------------------------------------------
   SICdb IDs LOCKED FROM THE v1.0.8 DICTIONARY / PRIOR STRUCTURE AUDIT
   -------------------------------------------------------------------------

   SOFA:
       3139  SOFA Score

   HR signals:
       707   HeartRateECG
       708   HeartRateSPO2
       724   HeartRateABP

   MAP signals:
       703   BloodPressureArterialMAP
       706   BloodPressureNIMAP

   Core-4 laboratory references:
       Sodium:
          469 central laboratory
          686 blood-gas candidate/fallback

       Potassium:
          463 central laboratory
          685 blood-gas candidate/fallback

       Chloride:
          450 central laboratory
          683 blood-gas candidate/fallback

       Bicarbonate:
          456 central laboratory
          666 blood-gas candidate/fallback

   LAB SOURCE RULE:
       within each 6-h bin:
         use median CENTRAL-lab value if available;
         otherwise use median blood-gas fallback value.

       Central and blood-gas values are NOT pooled when central laboratory
       data are available in that bin.

   Known CRRT/KRT signal references:
       723
       730
       731
       732
       2022
       3071

   In addition, this script dynamically includes data_float_h/data_range
   references whose dictionary labels contain:
       CRRT, CVVH, dialysis, haemodialysis/hemodialysis, hemofiltration

   ========================================================================= */


/* =========================================================================
   0. OUTPUT SCHEMA + CLEANUP
   ========================================================================= */

CREATE SCHEMA IF NOT EXISTS study_make30_final;

DROP TABLE IF EXISTS study_make30_final.sic_analysis_master;
DROP TABLE IF EXISTS study_make30_final.sic_make30;
DROP TABLE IF EXISTS study_make30_final.sic_covariates;
DROP TABLE IF EXISTS study_make30_final.sic_baseline_kidney;
DROP TABLE IF EXISTS study_make30_final.sic_case_extra;
DROP TABLE IF EXISTS study_make30_final.sic_src_data_ref;
DROP TABLE IF EXISTS study_make30_final.sic_creatinine_ids;
DROP TABLE IF EXISTS study_make30_final.sic_qc_summary;


DROP TABLE IF EXISTS study_make30_final.sic_msi_gt4_trace;
DROP TABLE IF EXISTS study_make30_final.sic_core4_multitrajectory_wide;
DROP TABLE IF EXISTS study_make30_final.sic_core4_electrolyte_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.sic_core4_msi_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.sic_core4_trajectory_cohort;
DROP TABLE IF EXISTS study_make30_final.sic_trajectory_coverage;
DROP TABLE IF EXISTS study_make30_final.sic_electrolyte_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.sic_msi_trajectory_long;
DROP TABLE IF EXISTS study_make30_final.sic_electrolyte_6h;
DROP TABLE IF EXISTS study_make30_final.sic_msi_2h;
DROP TABLE IF EXISTS study_make30_final.sic_msi_hourly;
DROP TABLE IF EXISTS study_make30_final.sic_trajectory_pool;
DROP TABLE IF EXISTS study_make30_final.sic_landmark_cohort;
DROP TABLE IF EXISTS study_make30_final.sic_sepsis_index;
DROP TABLE IF EXISTS study_make30_final.sic_sepsis_adult;
DROP TABLE IF EXISTS study_make30_final.sic_sepsis_first_case;
DROP TABLE IF EXISTS study_make30_final.sic_sepsis_pairs;
DROP TABLE IF EXISTS study_make30_final.sic_sofa_ge2;
DROP TABLE IF EXISTS study_make30_final.sic_suspected_infection;
DROP TABLE IF EXISTS study_make30_final.sic_antibiotic_admin;
DROP TABLE IF EXISTS study_make30_final.sic_antibiotic_patterns;
DROP TABLE IF EXISTS study_make30_final.sic_krt_ids;
DROP TABLE IF EXISTS study_make30_final.sic_lab_mapping;
DROP TABLE IF EXISTS study_make30_final.sic_signal_mapping;
DROP TABLE IF EXISTS study_make30_final.sic_cases;
DROP TABLE IF EXISTS study_make30_final.sic_patients;
DROP TABLE IF EXISTS study_make30_final.sic_src_data_range;
DROP TABLE IF EXISTS study_make30_final.sic_src_laboratory;
DROP TABLE IF EXISTS study_make30_final.sic_src_medication;
DROP TABLE IF EXISTS study_make30_final.sic_src_data_float_h;
DROP TABLE IF EXISTS study_make30_final.sic_src_refs;
DROP TABLE IF EXISTS study_make30_final.sic_src_column_map;
DROP TABLE IF EXISTS study_make30_final.sic_src_cases;


/* =========================================================================
   1. STANDARDIZE LOCAL SICdb TABLE/COLUMN NAMES

   This block searches the local PostgreSQL installation for one schema
   containing:
       cases
       d_references
       data_float_h
       medication
       laboratory

   It then creates standardized study_make30_final.sic_src_* tables so the rest
   of the script does not depend on upper/lower-case CSV header handling.

   laboratory reference column:
       accepts LaboratoryID (official SQL examples)
       and DrugID (appears in some documentation/build variants)

   data_range is optional for the KRT signal check.
   ========================================================================= */

CREATE TABLE study_make30_final.sic_src_column_map
(
    source_table    text NOT NULL,
    logical_field  text NOT NULL,
    physical_column text NOT NULL,
    mapping_method text NOT NULL
);


DO $$
DECLARE
    v_schema              text;

    c_caseid              text;
    c_patientid           text;
    c_timeofstay          text;
    c_icuoffset           text;
    c_age                 text;
    c_offsetfirst         text;
    c_offsetdeath         text;
    c_hospday             text;
    c_hospstaydays        text;

    r_id                  text;
    r_value               text;
    r_name                text;
    r_unit                text;

    f_caseid              text;
    f_dataid              text;
    f_offset              text;
    f_val                 text;

    f_caseid_method       text;
    f_dataid_method       text;
    f_offset_method       text;
    f_val_method          text;

    f_table_actual        text;

    m_caseid              text;
    m_drugid              text;
    m_offset              text;
    m_offsetend           text;
    m_single              text;
    m_amount              text;
    m_amountpm            text;

    m_caseid_method       text;
    m_drugid_method       text;
    m_offset_method       text;
    m_offsetend_method    text;
    m_single_method       text;
    m_amount_method       text;
    m_amountpm_method     text;

    m_table_actual        text;

    l_caseid              text;
    l_patientid           text;
    l_labid               text;
    l_offset              text;
    l_value               text;
    l_type                text;

    l_caseid_method       text;
    l_patientid_method    text;
    l_labid_method        text;
    l_offset_method       text;
    l_value_method        text;
    l_type_method         text;

    l_table_actual        text;

    dr_caseid             text;
    dr_dataid             text;
    dr_offset             text;
    dr_offsetend          text;

    dr_caseid_method      text;
    dr_dataid_method      text;
    dr_offset_method      text;
    dr_offsetend_method   text;

    dr_table_actual       text;
BEGIN

    /* ---------------------------------------------------------------------
       Identify SICdb schema.
       --------------------------------------------------------------------- */

    SELECT n.nspname
    INTO v_schema

    FROM pg_namespace n

    WHERE EXISTS
          (
              SELECT 1
              FROM pg_class c
              WHERE c.relnamespace = n.oid
                AND LOWER(c.relname) = 'cases'
                AND c.relkind IN ('r','p','v','m')
          )

      AND EXISTS
          (
              SELECT 1
              FROM pg_class c
              WHERE c.relnamespace = n.oid
                AND LOWER(c.relname) = 'd_references'
                AND c.relkind IN ('r','p','v','m')
          )

      AND EXISTS
          (
              SELECT 1
              FROM pg_class c
              WHERE c.relnamespace = n.oid
                AND LOWER(c.relname) = 'data_float_h'
                AND c.relkind IN ('r','p','v','m')
          )

      AND EXISTS
          (
              SELECT 1
              FROM pg_class c
              WHERE c.relnamespace = n.oid
                AND LOWER(c.relname) = 'medication'
                AND c.relkind IN ('r','p','v','m')
          )

      AND EXISTS
          (
              SELECT 1
              FROM pg_class c
              WHERE c.relnamespace = n.oid
                AND LOWER(c.relname) = 'laboratory'
                AND c.relkind IN ('r','p','v','m')
          )

      /* The schema must also have the essential SICdb CASES columns,
         allowing CamelCase, lowercase and snake_case imports. */
      AND
      (
          SELECT COUNT(DISTINCT
              REGEXP_REPLACE(
                  LOWER(ic.column_name),
                  '[^a-z0-9]',
                  '',
                  'g'
              )
          )
          FROM information_schema.columns ic
          WHERE ic.table_schema = n.nspname
            AND LOWER(ic.table_name) = 'cases'
            AND REGEXP_REPLACE(
                    LOWER(ic.column_name),
                    '[^a-z0-9]',
                    '',
                    'g'
                )
                IN
                (
                    'caseid',
                    'patientid',
                    'timeofstay',
                    'icuoffset',
                    'ageonadmission'
                )
      ) = 5

    ORDER BY
        CASE WHEN n.nspname = 'public' THEN 0 ELSE 1 END,
        n.nspname

    LIMIT 1;


    IF v_schema IS NULL THEN
        RAISE EXCEPTION
        'Could not find a SICdb schema containing cases, d_references, data_float_h, medication and laboratory WITH the essential CASES fields (CaseID, PatientID, TimeOfStay, ICUOffset, AgeOnAdmission), after normalization for snake_case/CamelCase names.';
    END IF;


    /* ---------------------------------------------------------------------
       Helper: CASES actual column names.
       --------------------------------------------------------------------- */

    SELECT column_name INTO c_caseid
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'cases'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'caseid'
    LIMIT 1;

    SELECT column_name INTO c_patientid
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'cases'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'patientid'
    LIMIT 1;

    SELECT column_name INTO c_timeofstay
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'cases'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'timeofstay'
    LIMIT 1;

    SELECT column_name INTO c_icuoffset
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'cases'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'icuoffset'
    LIMIT 1;

    SELECT column_name INTO c_age
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'cases'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'ageonadmission'
    LIMIT 1;

    SELECT column_name INTO c_offsetfirst
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'cases'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'offsetafterfirstadmission'
    LIMIT 1;

    SELECT column_name INTO c_offsetdeath
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'cases'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'offsetofdeath'
    LIMIT 1;

    SELECT column_name INTO c_hospday
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'cases'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'hospitaldischargeday'
    LIMIT 1;

    SELECT column_name INTO c_hospstaydays
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'cases'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'hospitalstaydays'
    LIMIT 1;


    IF c_caseid IS NULL
       OR c_patientid IS NULL
       OR c_timeofstay IS NULL
       OR c_icuoffset IS NULL
       OR c_age IS NULL
    THEN
        RAISE EXCEPTION
        'Selected SICdb cases table could not map one of the required logical fields after normalized-name matching: CaseID, PatientID, TimeOfStay, ICUOffset, AgeOnAdmission.';
    END IF;


    EXECUTE FORMAT(
        $q$
        CREATE TABLE study_make30_final.sic_src_cases AS
        SELECT
            %I::bigint          AS case_id,
            %I::bigint          AS patient_id,
            %I::bigint          AS time_of_stay_sec,
            COALESCE(%I,0)::bigint AS icu_offset_sec,
            %I::numeric         AS age_years,
            %s                  AS offset_after_first_admission_sec,
            %s                  AS offset_of_death_sec,
            %s                  AS hospital_discharge_day,
            %s                  AS hospital_stay_days
        FROM %I.%I
        WHERE %I IS NOT NULL
          AND %I IS NOT NULL
        $q$,

        c_caseid,
        c_patientid,
        c_timeofstay,
        c_icuoffset,
        c_age,

        CASE
            WHEN c_offsetfirst IS NULL
            THEN '0::bigint'
            ELSE FORMAT('COALESCE(%I,0)::bigint', c_offsetfirst)
        END,

        CASE
            WHEN c_offsetdeath IS NULL
            THEN 'NULL::bigint'
            ELSE FORMAT('%I::bigint', c_offsetdeath)
        END,

        CASE
            WHEN c_hospday IS NULL
            THEN 'NULL::integer'
            ELSE FORMAT('%I::integer', c_hospday)
        END,

        CASE
            WHEN c_hospstaydays IS NULL
            THEN 'NULL::integer'
            ELSE FORMAT('%I::integer', c_hospstaydays)
        END,

        v_schema,
        (
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = v_schema
              AND LOWER(table_name) = 'cases'
            LIMIT 1
        ),

        c_caseid,
        c_patientid
    );


    /* ---------------------------------------------------------------------
       REFERENCES.
       --------------------------------------------------------------------- */

    SELECT column_name INTO r_id
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'd_references'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'referenceglobalid'
    LIMIT 1;

    SELECT column_name INTO r_value
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'd_references'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'referencevalue'
    LIMIT 1;

    SELECT column_name INTO r_name
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'd_references'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'referencename'
    LIMIT 1;

    SELECT column_name INTO r_unit
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'd_references'
      AND REGEXP_REPLACE(LOWER(column_name), '[^a-z0-9]', '', 'g') = 'referenceunit'
    LIMIT 1;


    IF r_id IS NULL OR r_value IS NULL THEN
        RAISE EXCEPTION
        'SICdb d_references table is missing ReferenceGlobalID or ReferenceValue.';
    END IF;


    EXECUTE FORMAT(
        $q$
        CREATE TABLE study_make30_final.sic_src_refs AS
        SELECT
            %I::integer AS ref_id,
            %I::text    AS ref_value,
            %s          AS ref_name,
            %s          AS ref_unit
        FROM %I.%I
        $q$,

        r_id,
        r_value,

        CASE
            WHEN r_name IS NULL
            THEN 'NULL::text'
            ELSE FORMAT('%I::text', r_name)
        END,

        CASE
            WHEN r_unit IS NULL
            THEN 'NULL::text'
            ELSE FORMAT('%I::text', r_unit)
        END,

        v_schema,
        (
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = v_schema
              AND LOWER(table_name) = 'd_references'
            LIMIT 1
        )
    );


    /* ---------------------------------------------------------------------
       DATA_FLOAT_H.

       Official SICdb logical schema:
           CaseID | DataID | Offset | Val | cnt | rawdata

       Local PostgreSQL imports may rename these columns. We therefore use
       TWO-STAGE mapping:

       Stage 1:
           normalized semantic aliases

       Stage 2:
           official SICdb column-order fallback:
               1 = CaseID
               2 = DataID
               3 = Offset
               4 = Val

       The order fallback is used ONLY inside the table named data_float_h,
       whose official schema is fixed by SICdb documentation.

       The selected physical column names and mapping method are written to:
           study_make30_final.sic_src_column_map
       --------------------------------------------------------------------- */

    SELECT table_name
    INTO f_table_actual

    FROM information_schema.tables

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'data_float_h'

    LIMIT 1;


    /* ---------------- CaseID ---------------- */

    SELECT column_name
    INTO f_caseid

    FROM information_schema.columns

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'data_float_h'

      AND REGEXP_REPLACE(
              LOWER(column_name),
              '[^a-z0-9]',
              '',
              'g'
          )
          IN
          (
              'caseid',
              'case',
              'admissionid',
              'stayid'
          )

    ORDER BY
        CASE
            WHEN REGEXP_REPLACE(
                     LOWER(column_name),
                     '[^a-z0-9]',
                     '',
                     'g'
                 ) = 'caseid'
            THEN 0
            ELSE 1
        END,
        ordinal_position

    LIMIT 1;


    IF f_caseid IS NOT NULL THEN
        f_caseid_method := 'normalized semantic name';
    ELSE
        SELECT column_name
        INTO f_caseid

        FROM information_schema.columns

        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'data_float_h'
          AND ordinal_position = 1

        LIMIT 1;

        f_caseid_method := 'official SICdb ordinal fallback #1';
    END IF;


    /* ---------------- DataID ---------------- */

    SELECT column_name
    INTO f_dataid

    FROM information_schema.columns

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'data_float_h'

      AND REGEXP_REPLACE(
              LOWER(column_name),
              '[^a-z0-9]',
              '',
              'g'
          )
          IN
          (
              'dataid',
              'fieldid',
              'referenceid',
              'referenceglobalid',
              'variableid',
              'itemid'
          )

    ORDER BY
        CASE REGEXP_REPLACE(
                 LOWER(column_name),
                 '[^a-z0-9]',
                 '',
                 'g'
             )
            WHEN 'dataid' THEN 0
            WHEN 'fieldid' THEN 1
            WHEN 'referenceglobalid' THEN 2
            WHEN 'referenceid' THEN 3
            ELSE 9
        END,
        ordinal_position

    LIMIT 1;


    IF f_dataid IS NOT NULL THEN
        f_dataid_method := 'normalized semantic name';
    ELSE
        SELECT column_name
        INTO f_dataid

        FROM information_schema.columns

        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'data_float_h'
          AND ordinal_position = 2

        LIMIT 1;

        f_dataid_method := 'official SICdb ordinal fallback #2';
    END IF;


    /* ---------------- Offset ---------------- */

    SELECT column_name
    INTO f_offset

    FROM information_schema.columns

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'data_float_h'

      AND REGEXP_REPLACE(
              LOWER(column_name),
              '[^a-z0-9]',
              '',
              'g'
          )
          IN
          (
              'offset',
              'offsetvalue',
              'offsetsec',
              'timeoffset',
              'chartoffset',
              'secondsfromadmission',
              'time'
          )

    ORDER BY
        CASE REGEXP_REPLACE(
                 LOWER(column_name),
                 '[^a-z0-9]',
                 '',
                 'g'
             )
            WHEN 'offset' THEN 0
            WHEN 'offsetvalue' THEN 1
            WHEN 'offsetsec' THEN 2
            WHEN 'timeoffset' THEN 3
            WHEN 'chartoffset' THEN 4
            ELSE 9
        END,
        ordinal_position

    LIMIT 1;


    IF f_offset IS NOT NULL THEN
        f_offset_method := 'normalized semantic name';
    ELSE
        RAISE EXCEPTION
        'Could not identify data_float_h time offset by semantic name. Refusing ordinal guessing for a critical time field.';
    END IF;


    /* ---------------- Val ---------------- */

    SELECT column_name
    INTO f_val

    FROM information_schema.columns

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'data_float_h'

      AND REGEXP_REPLACE(
              LOWER(column_name),
              '[^a-z0-9]',
              '',
              'g'
          )
          IN
          (
              'val',
              'value',
              'valuefloat',
              'floatvalue',
              'mean',
              'average',
              'avg'
          )

    ORDER BY
        CASE REGEXP_REPLACE(
                 LOWER(column_name),
                 '[^a-z0-9]',
                 '',
                 'g'
             )
            WHEN 'val' THEN 0
            WHEN 'value' THEN 1
            WHEN 'valuefloat' THEN 2
            WHEN 'floatvalue' THEN 3
            WHEN 'mean' THEN 4
            WHEN 'average' THEN 5
            WHEN 'avg' THEN 6
            ELSE 9
        END,
        ordinal_position

    LIMIT 1;


    IF f_val IS NOT NULL THEN
        f_val_method := 'normalized semantic name';
    ELSE
        SELECT column_name
        INTO f_val

        FROM information_schema.columns

        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'data_float_h'
          AND ordinal_position = 4

        LIMIT 1;

        f_val_method := 'official SICdb ordinal fallback #4';
    END IF;


    IF f_table_actual IS NULL
       OR f_caseid IS NULL
       OR f_dataid IS NULL
       OR f_offset IS NULL
       OR f_val IS NULL
    THEN

        RAISE EXCEPTION
        'Could not map the four required SICdb data_float_h logical fields (CaseID, DataID, Offset, Val), even after semantic-name and official-ordinal fallback mapping.';

    END IF;


    /* Record the exact local mapping for later auditing. */

    INSERT INTO study_make30_final.sic_src_column_map
    (source_table, logical_field, physical_column, mapping_method)
    VALUES
    ('data_float_h', 'CaseID', f_caseid, f_caseid_method),
    ('data_float_h', 'DataID', f_dataid, f_dataid_method),
    ('data_float_h', 'Offset', f_offset, f_offset_method),
    ('data_float_h', 'Val',    f_val,    f_val_method);


    EXECUTE FORMAT(
        $q$
        CREATE TABLE study_make30_final.sic_src_data_float_h AS
        SELECT
            %I::bigint             AS case_id,
            %I::integer            AS data_id,
            %I::bigint             AS offset_sec,
            %I::double precision   AS val
        FROM %I.%I
        WHERE %I IS NOT NULL
          AND %I IS NOT NULL
        $q$,

        f_caseid,
        f_dataid,
        f_offset,
        f_val,

        v_schema,
        f_table_actual,

        f_caseid,
        f_dataid
    );


    /* ---------------------------------------------------------------------
       MEDICATION.

       Official SICdb order:
           #1 id
           #2 CaseID
           #3 DrugID
           #4 Offset
           #5 OffsetDrugEnd
           #6 IsSingleDose
           #7 Amount
           #8 AmountPerMinute

       Two-stage mapping:
           1) normalized semantic column name
           2) official SICdb ordinal fallback

       Exact physical mappings are written to sic_src_column_map.
       --------------------------------------------------------------------- */

    SELECT table_name
    INTO m_table_actual

    FROM information_schema.tables

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'medication'

    LIMIT 1;


    /* ---------------- CaseID ---------------- */

    SELECT column_name
    INTO m_caseid

    FROM information_schema.columns

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'medication'

      AND REGEXP_REPLACE(
              LOWER(column_name),
              '[^a-z0-9]',
              '',
              'g'
          )
          IN ('caseid','case','admissionid','stayid')

    ORDER BY
        CASE
            WHEN REGEXP_REPLACE(
                     LOWER(column_name),
                     '[^a-z0-9]',
                     '',
                     'g'
                 ) = 'caseid'
            THEN 0
            ELSE 1
        END,
        ordinal_position

    LIMIT 1;


    IF m_caseid IS NOT NULL THEN
        m_caseid_method := 'normalized semantic name';
    ELSE
        SELECT column_name
        INTO m_caseid

        FROM information_schema.columns

        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'medication'
          AND ordinal_position = 2

        LIMIT 1;

        m_caseid_method := 'official SICdb ordinal fallback #2';
    END IF;


    /* ---------------- DrugID ---------------- */

    SELECT column_name
    INTO m_drugid

    FROM information_schema.columns

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'medication'

      AND REGEXP_REPLACE(
              LOWER(column_name),
              '[^a-z0-9]',
              '',
              'g'
          )
          IN
          (
              'drugid',
              'medicationid',
              'referenceid',
              'referenceglobalid',
              'itemid'
          )

    ORDER BY
        CASE REGEXP_REPLACE(
                 LOWER(column_name),
                 '[^a-z0-9]',
                 '',
                 'g'
             )
            WHEN 'drugid' THEN 0
            WHEN 'medicationid' THEN 1
            WHEN 'referenceglobalid' THEN 2
            WHEN 'referenceid' THEN 3
            ELSE 9
        END,
        ordinal_position

    LIMIT 1;


    IF m_drugid IS NOT NULL THEN
        m_drugid_method := 'normalized semantic name';
    ELSE
        SELECT column_name
        INTO m_drugid

        FROM information_schema.columns

        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'medication'
          AND ordinal_position = 3

        LIMIT 1;

        m_drugid_method := 'official SICdb ordinal fallback #3';
    END IF;


    /* ---------------- Offset ---------------- */

    SELECT column_name
    INTO m_offset

    FROM information_schema.columns

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'medication'

      AND REGEXP_REPLACE(
              LOWER(column_name),
              '[^a-z0-9]',
              '',
              'g'
          )
          IN
          (
              'offset',
              'offsetvalue',
              'offsetsec',
              'timeoffset',
              'startoffset',
              'drugstartoffset'
          )

    ORDER BY
        CASE REGEXP_REPLACE(
                 LOWER(column_name),
                 '[^a-z0-9]',
                 '',
                 'g'
             )
            WHEN 'offset' THEN 0
            WHEN 'offsetvalue' THEN 1
            WHEN 'offsetsec' THEN 2
            WHEN 'timeoffset' THEN 3
            WHEN 'startoffset' THEN 4
            WHEN 'drugstartoffset' THEN 5
            ELSE 9
        END,
        ordinal_position

    LIMIT 1;


    IF m_offset IS NOT NULL THEN
        m_offset_method := 'normalized semantic name';
    ELSE
        RAISE EXCEPTION
        'Could not identify medication administration time offset by semantic name. Refusing ordinal guessing for a critical time field.';
    END IF;


    /* ---------------- OffsetDrugEnd ---------------- */

    SELECT column_name
    INTO m_offsetend

    FROM information_schema.columns

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'medication'

      AND REGEXP_REPLACE(
              LOWER(column_name),
              '[^a-z0-9]',
              '',
              'g'
          )
          IN
          (
              'offsetdrugend',
              'offsetend',
              'endoffset',
              'drugendoffset'
          )

    ORDER BY
        CASE REGEXP_REPLACE(
                 LOWER(column_name),
                 '[^a-z0-9]',
                 '',
                 'g'
             )
            WHEN 'offsetdrugend' THEN 0
            WHEN 'offsetend' THEN 1
            WHEN 'endoffset' THEN 2
            WHEN 'drugendoffset' THEN 3
            ELSE 9
        END,
        ordinal_position

    LIMIT 1;


    IF m_offsetend IS NOT NULL THEN
        m_offsetend_method := 'normalized semantic name';
    ELSE
        SELECT column_name
        INTO m_offsetend

        FROM information_schema.columns

        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'medication'
          AND ordinal_position = 5

        LIMIT 1;

        m_offsetend_method := 'official SICdb ordinal fallback #5';
    END IF;


    /* ---------------- IsSingleDose ---------------- */

    SELECT column_name
    INTO m_single

    FROM information_schema.columns

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'medication'

      AND REGEXP_REPLACE(
              LOWER(column_name),
              '[^a-z0-9]',
              '',
              'g'
          )
          IN
          (
              'issingledose',
              'singledose',
              'isbolus',
              'bolus'
          )

    ORDER BY ordinal_position
    LIMIT 1;


    IF m_single IS NOT NULL THEN
        m_single_method := 'normalized semantic name';
    ELSE
        SELECT column_name
        INTO m_single

        FROM information_schema.columns

        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'medication'
          AND ordinal_position = 6

        LIMIT 1;

        m_single_method := 'official SICdb ordinal fallback #6';
    END IF;


    /* ---------------- Amount ---------------- */

    SELECT column_name
    INTO m_amount

    FROM information_schema.columns

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'medication'

      AND REGEXP_REPLACE(
              LOWER(column_name),
              '[^a-z0-9]',
              '',
              'g'
          ) = 'amount'

    ORDER BY ordinal_position
    LIMIT 1;


    IF m_amount IS NOT NULL THEN
        m_amount_method := 'normalized semantic name';
    ELSE
        SELECT column_name
        INTO m_amount

        FROM information_schema.columns

        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'medication'
          AND ordinal_position = 7

        LIMIT 1;

        m_amount_method := 'official SICdb ordinal fallback #7';
    END IF;


    /* ---------------- AmountPerMinute ---------------- */

    SELECT column_name
    INTO m_amountpm

    FROM information_schema.columns

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'medication'

      AND REGEXP_REPLACE(
              LOWER(column_name),
              '[^a-z0-9]',
              '',
              'g'
          )
          IN ('amountperminute','rate','doseperminute')

    ORDER BY
        CASE REGEXP_REPLACE(
                 LOWER(column_name),
                 '[^a-z0-9]',
                 '',
                 'g'
             )
            WHEN 'amountperminute' THEN 0
            WHEN 'doseperminute' THEN 1
            WHEN 'rate' THEN 2
            ELSE 9
        END,
        ordinal_position

    LIMIT 1;


    IF m_amountpm IS NOT NULL THEN
        m_amountpm_method := 'normalized semantic name';
    ELSE
        SELECT column_name
        INTO m_amountpm

        FROM information_schema.columns

        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'medication'
          AND ordinal_position = 8

        LIMIT 1;

        m_amountpm_method := 'official SICdb ordinal fallback #8';
    END IF;


    IF m_table_actual IS NULL
       OR m_caseid IS NULL
       OR m_drugid IS NULL
       OR m_offset IS NULL
    THEN
        RAISE EXCEPTION
        'Could not map SICdb medication CaseID/DrugID/Offset after semantic-name and official-ordinal fallback.';
    END IF;


    INSERT INTO study_make30_final.sic_src_column_map
    (source_table, logical_field, physical_column, mapping_method)
    VALUES
    ('medication','CaseID',          m_caseid,    m_caseid_method),
    ('medication','DrugID',          m_drugid,    m_drugid_method),
    ('medication','Offset',          m_offset,    m_offset_method),
    ('medication','OffsetDrugEnd',   COALESCE(m_offsetend,'<not found>'),
                                     COALESCE(m_offsetend_method,'not mapped')),
    ('medication','IsSingleDose',    COALESCE(m_single,'<not found>'),
                                     COALESCE(m_single_method,'not mapped')),
    ('medication','Amount',          COALESCE(m_amount,'<not found>'),
                                     COALESCE(m_amount_method,'not mapped')),
    ('medication','AmountPerMinute', COALESCE(m_amountpm,'<not found>'),
                                     COALESCE(m_amountpm_method,'not mapped'));


    EXECUTE FORMAT(
        $q$
        CREATE TABLE study_make30_final.sic_src_medication AS
        SELECT
            %I::bigint   AS case_id,
            %I::integer  AS drug_id,
            %I::bigint   AS offset_sec,
            %s           AS offset_end_sec,
            %s           AS is_single_dose,
            %s           AS amount,
            %s           AS amount_per_minute
        FROM %I.%I
        WHERE %I IS NOT NULL
          AND %I IS NOT NULL
          AND %I IS NOT NULL
        $q$,

        m_caseid,
        m_drugid,
        m_offset,

        CASE
            WHEN m_offsetend IS NULL
            THEN 'NULL::bigint'
            ELSE FORMAT('%I::bigint', m_offsetend)
        END,

        CASE
            WHEN m_single IS NULL
            THEN 'NULL::integer'
            ELSE FORMAT('%I::integer', m_single)
        END,

        CASE
            WHEN m_amount IS NULL
            THEN 'NULL::double precision'
            ELSE FORMAT('%I::double precision', m_amount)
        END,

        CASE
            WHEN m_amountpm IS NULL
            THEN 'NULL::double precision'
            ELSE FORMAT('%I::double precision', m_amountpm)
        END,

        v_schema,
        m_table_actual,

        m_caseid,
        m_drugid,
        m_offset
    );


    /* ---------------------------------------------------------------------
       LABORATORY.

       Official SICdb order:
           #1 id
           #2 CaseID
           #3 PatientID
           #4 DrugID
           #5 Offset
           #6 LaboratoryValue
           #7 LaboratoryType

       Some prior local dictionaries label the reference field
       LaboratoryID. Both DrugID/LaboratoryID are accepted semantically.
       If name matching fails, official ordinal fallback is used.
       --------------------------------------------------------------------- */

    SELECT table_name
    INTO l_table_actual

    FROM information_schema.tables

    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'laboratory'

    LIMIT 1;


    /* CaseID */
    SELECT column_name INTO l_caseid
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'laboratory'
      AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
          IN ('caseid','case','admissionid','stayid')
    ORDER BY
        CASE WHEN REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='caseid'
             THEN 0 ELSE 1 END,
        ordinal_position
    LIMIT 1;

    IF l_caseid IS NOT NULL THEN
        l_caseid_method := 'normalized semantic name';
    ELSE
        SELECT column_name INTO l_caseid
        FROM information_schema.columns
        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'laboratory'
          AND ordinal_position = 2
        LIMIT 1;
        l_caseid_method := 'official SICdb ordinal fallback #2';
    END IF;


    /* PatientID */
    SELECT column_name INTO l_patientid
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'laboratory'
      AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
          IN ('patientid','patient')
    ORDER BY
        CASE WHEN REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='patientid'
             THEN 0 ELSE 1 END,
        ordinal_position
    LIMIT 1;

    IF l_patientid IS NOT NULL THEN
        l_patientid_method := 'normalized semantic name';
    ELSE
        l_patientid := NULL;
        l_patientid_method := 'not present/not required';
    END IF;


    /* DrugID / LaboratoryID reference */
    SELECT column_name INTO l_labid
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'laboratory'
      AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
          IN
          (
              'drugid',
              'laboratoryid',
              'labid',
              'referenceid',
              'referenceglobalid',
              'itemid'
          )
    ORDER BY
        CASE REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
            WHEN 'drugid' THEN 0
            WHEN 'laboratoryid' THEN 1
            WHEN 'labid' THEN 2
            WHEN 'referenceglobalid' THEN 3
            ELSE 9
        END,
        ordinal_position
    LIMIT 1;

    IF l_labid IS NOT NULL THEN
        l_labid_method := 'normalized semantic name';
    ELSE
        SELECT column_name INTO l_labid
        FROM information_schema.columns
        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'laboratory'
          AND ordinal_position = 4
        LIMIT 1;
        l_labid_method := 'official SICdb ordinal fallback #4';
    END IF;


    /* Offset */
    SELECT column_name INTO l_offset
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'laboratory'
      AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
          IN ('offset','offsetvalue','offsetsec','timeoffset','chartoffset','labresultoffset')
    ORDER BY
        CASE REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
            WHEN 'offset' THEN 0
            WHEN 'offsetvalue' THEN 1
            WHEN 'offsetsec' THEN 2
            WHEN 'timeoffset' THEN 3
            WHEN 'chartoffset' THEN 4
            ELSE 9
        END,
        ordinal_position
    LIMIT 1;

    IF l_offset IS NOT NULL THEN
        l_offset_method := 'normalized semantic name';
    ELSE
        RAISE EXCEPTION
        'Could not identify laboratory measurement time offset by semantic name. Refusing ordinal guessing for a critical time field.';
    END IF;


    /* LaboratoryValue */
    SELECT column_name INTO l_value
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'laboratory'
      AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
          IN ('laboratoryvalue','labvalue','value','val','labresult')
    ORDER BY
        CASE REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
            WHEN 'laboratoryvalue' THEN 0
            WHEN 'labvalue' THEN 1
            WHEN 'labresult' THEN 2
            WHEN 'value' THEN 3
            WHEN 'val' THEN 4
            ELSE 9
        END,
        ordinal_position
    LIMIT 1;

    IF l_value IS NOT NULL THEN
        l_value_method := 'normalized semantic name';
    ELSE
        SELECT column_name INTO l_value
        FROM information_schema.columns
        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'laboratory'
          AND ordinal_position = 6
        LIMIT 1;
        l_value_method := 'official SICdb ordinal fallback #6';
    END IF;


    /* LaboratoryType */
    SELECT column_name INTO l_type
    FROM information_schema.columns
    WHERE table_schema = v_schema
      AND LOWER(table_name) = 'laboratory'
      AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
          IN ('laboratorytype','labtype','type')
    ORDER BY ordinal_position
    LIMIT 1;

    IF l_type IS NOT NULL THEN
        l_type_method := 'normalized semantic name';
    ELSE
        SELECT column_name INTO l_type
        FROM information_schema.columns
        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'laboratory'
          AND ordinal_position = 7
        LIMIT 1;
        l_type_method := 'official SICdb ordinal fallback #7';
    END IF;


    IF l_table_actual IS NULL
       OR l_caseid IS NULL
       OR l_labid IS NULL
       OR l_offset IS NULL
       OR l_value IS NULL
    THEN
        RAISE EXCEPTION
        'Could not map SICdb laboratory CaseID/DrugID(or LaboratoryID)/Offset/LaboratoryValue after semantic-name and official-ordinal fallback.';
    END IF;


    INSERT INTO study_make30_final.sic_src_column_map
    (source_table, logical_field, physical_column, mapping_method)
    VALUES
    ('laboratory','CaseID',          l_caseid,    l_caseid_method),
    ('laboratory','PatientID',       COALESCE(l_patientid,'<not found>'),
                                    COALESCE(l_patientid_method,'not mapped')),
    ('laboratory','DrugID',          l_labid,     l_labid_method),
    ('laboratory','Offset',          l_offset,    l_offset_method),
    ('laboratory','LaboratoryValue', l_value,     l_value_method),
    ('laboratory','LaboratoryType',  COALESCE(l_type,'<not found>'),
                                    COALESCE(l_type_method,'not mapped'));


    EXECUTE FORMAT(
        $q$
        CREATE TABLE study_make30_final.sic_src_laboratory AS
        SELECT
            %I::bigint            AS case_id,
            %s                    AS patient_id,
            %I::integer           AS lab_id,
            %I::bigint            AS offset_sec,
            %I::double precision  AS lab_value,
            %s                    AS lab_type
        FROM %I.%I
        WHERE %I IS NOT NULL
          AND %I IS NOT NULL
          AND %I IS NOT NULL
          AND %I IS NOT NULL
        $q$,

        l_caseid,

        CASE
            WHEN l_patientid IS NULL
            THEN 'NULL::bigint'
            ELSE FORMAT('%I::bigint', l_patientid)
        END,

        l_labid,
        l_offset,
        l_value,

        CASE
            WHEN l_type IS NULL
            THEN 'NULL::integer'
            ELSE FORMAT('%I::integer', l_type)
        END,

        v_schema,
        l_table_actual,

        l_caseid,
        l_labid,
        l_offset,
        l_value
    );


    /* ---------------------------------------------------------------------
       DATA_RANGE — optional but used when available for KRT activity.

       Official SICdb order:
           #1 id
           #2 CaseID
           #3 DataID
           #4 Offset
           #5 OffsetEnd
           #6 Data

       Semantic-name mapping is attempted first, then official ordinal
       fallback. If data_range does not exist, an empty standardized table
       is created.
       --------------------------------------------------------------------- */

    IF EXISTS
    (
        SELECT 1
        FROM information_schema.tables
        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'data_range'
    )
    THEN

        SELECT table_name
        INTO dr_table_actual

        FROM information_schema.tables

        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'data_range'

        LIMIT 1;


        /* CaseID */
        SELECT column_name INTO dr_caseid
        FROM information_schema.columns
        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'data_range'
          AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
              IN ('caseid','case','admissionid','stayid')
        ORDER BY
            CASE WHEN REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='caseid'
                 THEN 0 ELSE 1 END,
            ordinal_position
        LIMIT 1;

        IF dr_caseid IS NOT NULL THEN
            dr_caseid_method := 'normalized semantic name';
        ELSE
            SELECT column_name INTO dr_caseid
            FROM information_schema.columns
            WHERE table_schema = v_schema
              AND LOWER(table_name) = 'data_range'
              AND ordinal_position = 2
            LIMIT 1;
            dr_caseid_method := 'official SICdb ordinal fallback #2';
        END IF;


        /* DataID */
        SELECT column_name INTO dr_dataid
        FROM information_schema.columns
        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'data_range'
          AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
              IN ('dataid','fieldid','referenceid','referenceglobalid','itemid')
        ORDER BY
            CASE REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
                WHEN 'dataid' THEN 0
                WHEN 'fieldid' THEN 1
                WHEN 'referenceglobalid' THEN 2
                ELSE 9
            END,
            ordinal_position
        LIMIT 1;

        IF dr_dataid IS NOT NULL THEN
            dr_dataid_method := 'normalized semantic name';
        ELSE
            SELECT column_name INTO dr_dataid
            FROM information_schema.columns
            WHERE table_schema = v_schema
              AND LOWER(table_name) = 'data_range'
              AND ordinal_position = 3
            LIMIT 1;
            dr_dataid_method := 'official SICdb ordinal fallback #3';
        END IF;


        /* Offset */
        SELECT column_name INTO dr_offset
        FROM information_schema.columns
        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'data_range'
          AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
              IN ('offset','offsetvalue','offsetsec','timeoffset','startoffset')
        ORDER BY
            CASE REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
                WHEN 'offset' THEN 0
                WHEN 'offsetvalue' THEN 1
                WHEN 'offsetsec' THEN 2
                WHEN 'timeoffset' THEN 3
                WHEN 'startoffset' THEN 4
                ELSE 9
            END,
            ordinal_position
        LIMIT 1;

        IF dr_offset IS NOT NULL THEN
            dr_offset_method := 'normalized semantic name';
        ELSE
            RAISE EXCEPTION
            'Could not identify data_range start offset by semantic name. Refusing ordinal guessing for a critical time field.';
        END IF;


        /* OffsetEnd */
        SELECT column_name INTO dr_offsetend
        FROM information_schema.columns
        WHERE table_schema = v_schema
          AND LOWER(table_name) = 'data_range'
          AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')
              IN ('offsetend','endoffset','endtimeoffset')
        ORDER BY ordinal_position
        LIMIT 1;

        IF dr_offsetend IS NOT NULL THEN
            dr_offsetend_method := 'normalized semantic name';
        ELSE
            SELECT column_name INTO dr_offsetend
            FROM information_schema.columns
            WHERE table_schema = v_schema
              AND LOWER(table_name) = 'data_range'
              AND ordinal_position = 5
            LIMIT 1;
            dr_offsetend_method := 'official SICdb ordinal fallback #5';
        END IF;


        IF dr_table_actual IS NOT NULL
           AND dr_caseid IS NOT NULL
           AND dr_dataid IS NOT NULL
           AND dr_offset IS NOT NULL
        THEN

            INSERT INTO study_make30_final.sic_src_column_map
            (source_table, logical_field, physical_column, mapping_method)
            VALUES
            ('data_range','CaseID',    dr_caseid, dr_caseid_method),
            ('data_range','DataID',    dr_dataid, dr_dataid_method),
            ('data_range','Offset',    dr_offset, dr_offset_method),
            ('data_range','OffsetEnd', COALESCE(dr_offsetend,'<not found>'),
                                      COALESCE(dr_offsetend_method,'not mapped'));


            EXECUTE FORMAT(
                $q$
                CREATE TABLE study_make30_final.sic_src_data_range AS
                SELECT
                    %I::bigint    AS case_id,
                    %I::integer   AS data_id,
                    %I::bigint    AS offset_sec,
                    %s            AS offset_end_sec
                FROM %I.%I
                WHERE %I IS NOT NULL
                  AND %I IS NOT NULL
                  AND %I IS NOT NULL
                $q$,

                dr_caseid,
                dr_dataid,
                dr_offset,

                CASE
                    WHEN dr_offsetend IS NULL
                    THEN FORMAT('%I::bigint', dr_offset)
                    ELSE FORMAT('%I::bigint', dr_offsetend)
                END,

                v_schema,
                dr_table_actual,

                dr_caseid,
                dr_dataid,
                dr_offset
            );

        ELSE

            CREATE TABLE study_make30_final.sic_src_data_range
            (
                case_id        bigint,
                data_id        integer,
                offset_sec     bigint,
                offset_end_sec bigint
            );

        END IF;

    ELSE

        CREATE TABLE study_make30_final.sic_src_data_range
        (
            case_id        bigint,
            data_id        integer,
            offset_sec     bigint,
            offset_end_sec bigint
        );

    END IF;


    RAISE NOTICE
    'SICdb standardized from local schema: %', v_schema;

END $$;


/* source indexes */

CREATE UNIQUE INDEX idx_sic_src_cases_case
ON study_make30_final.sic_src_cases(case_id);

CREATE INDEX idx_sic_src_cases_patient
ON study_make30_final.sic_src_cases(patient_id);

CREATE UNIQUE INDEX idx_sic_src_refs
ON study_make30_final.sic_src_refs(ref_id);

CREATE INDEX idx_sic_src_float_case_data_time
ON study_make30_final.sic_src_data_float_h(case_id, data_id, offset_sec);

CREATE INDEX idx_sic_src_med_case_drug_time
ON study_make30_final.sic_src_medication(case_id, drug_id, offset_sec);

CREATE INDEX idx_sic_src_lab_case_lab_time
ON study_make30_final.sic_src_laboratory(case_id, lab_id, offset_sec);

CREATE INDEX idx_sic_src_range_case_data_time
ON study_make30_final.sic_src_data_range(case_id, data_id, offset_sec);


/* =========================================================================
   1A. MANDATORY SOURCE-MAPPING VALIDATION

   Purpose:
   Prevent a wrong local-column mapping from silently producing an empty
   trajectory dataset.

   The script stops if:
     - time offset was accidentally mapped to DataID / DrugID / LabValue;
     - key SICdb signals have no temporal variation;
     - Core-4 laboratory IDs have no temporal variation.

   This validation occurs BEFORE the full research cohort is constructed.
   ========================================================================= */

DO $$
DECLARE
    v_df_distinct_offsets   bigint;
    v_df_distinct_dataids   bigint;
    v_med_distinct_offsets  bigint;
    v_lab_distinct_offsets  bigint;
    v_vital_rows            bigint;
    v_lab_rows              bigint;
BEGIN

    SELECT
        COUNT(DISTINCT offset_sec),
        COUNT(DISTINCT data_id),
        COUNT(*) FILTER
        (
            WHERE data_id IN
            (
                3139,
                707,708,724,
                703,706
            )
        )
    INTO
        v_df_distinct_offsets,
        v_df_distinct_dataids,
        v_vital_rows
    FROM study_make30_final.sic_src_data_float_h;


    SELECT
        COUNT(DISTINCT offset_sec)
    INTO
        v_med_distinct_offsets
    FROM study_make30_final.sic_src_medication;


    SELECT
        COUNT(DISTINCT offset_sec),
        COUNT(*) FILTER
        (
            WHERE lab_id IN
            (
                469,686,
                463,685,
                450,683,
                456,666
            )
        )
    INTO
        v_lab_distinct_offsets,
        v_lab_rows
    FROM study_make30_final.sic_src_laboratory;


    IF v_df_distinct_offsets IS NULL
       OR v_df_distinct_offsets < 100
    THEN
        RAISE EXCEPTION
        'SOURCE VALIDATION FAILED: data_float_h offset has too little temporal variation (% distinct values). Check the physical Offset mapping.',
        COALESCE(v_df_distinct_offsets,0);
    END IF;


    IF v_df_distinct_offsets = v_df_distinct_dataids
       AND EXISTS
       (
           SELECT 1
           FROM study_make30_final.sic_src_data_float_h
           WHERE offset_sec = data_id
           LIMIT 1
       )
    THEN
        RAISE EXCEPTION
        'SOURCE VALIDATION FAILED: data_float_h offset appears to be mapped to DataID.';
    END IF;


    IF v_med_distinct_offsets IS NULL
       OR v_med_distinct_offsets < 100
    THEN
        RAISE EXCEPTION
        'SOURCE VALIDATION FAILED: medication offset has too little temporal variation (% distinct values). Check the physical medication Offset mapping.',
        COALESCE(v_med_distinct_offsets,0);
    END IF;


    IF v_lab_distinct_offsets IS NULL
       OR v_lab_distinct_offsets < 100
    THEN
        RAISE EXCEPTION
        'SOURCE VALIDATION FAILED: laboratory offset has too little temporal variation (% distinct values). Check the physical laboratory Offset mapping.',
        COALESCE(v_lab_distinct_offsets,0);
    END IF;


    IF COALESCE(v_vital_rows,0) = 0
    THEN
        RAISE EXCEPTION
        'SOURCE VALIDATION FAILED: no rows found for SICdb SOFA/HR/MAP IDs 3139,707,708,724,703,706.';
    END IF;


    IF COALESCE(v_lab_rows,0) = 0
    THEN
        RAISE EXCEPTION
        'SOURCE VALIDATION FAILED: no rows found for Core-4 lab IDs 469/686,463/685,450/683,456/666.';
    END IF;


    RAISE NOTICE
    'SOURCE VALIDATION PASSED: data_float_h offsets=%, medication offsets=%, laboratory offsets=%, vital/SOFA rows=%, Core-4 lab rows=%',
    v_df_distinct_offsets,
    v_med_distinct_offsets,
    v_lab_distinct_offsets,
    v_vital_rows,
    v_lab_rows;

END $$;


/* =========================================================================
   2. FULL SICdb PATIENT COHORT

   Unit = PatientID
   All CaseID admissions for all patients are retained.
   ========================================================================= */

CREATE TABLE study_make30_final.sic_patients AS

SELECT DISTINCT
    patient_id

FROM study_make30_final.sic_src_cases

WHERE patient_id IS NOT NULL;


CREATE UNIQUE INDEX idx_sic_patients
ON study_make30_final.sic_patients(patient_id);


/* =========================================================================
   3. ALL CASES FOR THE FULL PATIENT COHORT
   ========================================================================= */

CREATE TABLE study_make30_final.sic_cases AS

SELECT
    c.*

FROM study_make30_final.sic_src_cases c

INNER JOIN study_make30_final.sic_patients p
    ON c.patient_id = p.patient_id;


CREATE UNIQUE INDEX idx_sic_cases_case
ON study_make30_final.sic_cases(case_id);

CREATE INDEX idx_sic_cases_patient
ON study_make30_final.sic_cases(patient_id);


/* =========================================================================
   4. LOCK / AUDIT SIGNAL MAPPING
   ========================================================================= */

CREATE TABLE study_make30_final.sic_signal_mapping
(
    variable_name   text,
    source_priority integer,
    data_id         integer,
    source_type     text
);

INSERT INTO study_make30_final.sic_signal_mapping
(variable_name, source_priority, data_id, source_type)
VALUES
('SOFA',        1, 3139, 'data_float_h'),

('HR',          1,  707, 'data_float_h'),
('HR',          2,  708, 'data_float_h'),
('HR',          3,  724, 'data_float_h'),

('MAP',         1,  703, 'data_float_h'),
('MAP',         2,  706, 'data_float_h');


/* =========================================================================
   5. LOCK / AUDIT CORE-4 LAB MAPPING

   source_priority:
       1 = central laboratory
       2 = blood-gas fallback
   ========================================================================= */

CREATE TABLE study_make30_final.sic_lab_mapping
(
    variable_name   text,
    source_priority integer,
    lab_id          integer
);

INSERT INTO study_make30_final.sic_lab_mapping
(variable_name, source_priority, lab_id)
VALUES
('sodium',       1, 469),
('sodium',       2, 686),

('potassium',    1, 463),
('potassium',    2, 685),

('chloride',     1, 450),
('chloride',     2, 683),

('bicarbonate',  1, 456),
('bicarbonate',  2, 666);


/* =========================================================================
   6. KRT / CRRT DICTIONARY

   Starts with previously confirmed SICdb CRRT references, then augments them
   dynamically using reference labels.
   ========================================================================= */

CREATE TABLE study_make30_final.sic_krt_ids AS

SELECT DISTINCT
    r.ref_id AS data_id,
    r.ref_value,
    r.ref_name,
    r.ref_unit

FROM study_make30_final.sic_src_refs r

WHERE r.ref_id IN
(
    723,
    730,
    731,
    732,
    2022,
    3071
)

OR LOWER(COALESCE(r.ref_value,'')) ~
   '(crrt|cvvh|cvvhd|cvvhdf|dialysis|dialyse|hemodial|haemodial|hemofiltr|haemofiltr)';


CREATE UNIQUE INDEX idx_sic_krt_ids
ON study_make30_final.sic_krt_ids(data_id);


/* =========================================================================
   7. TRANSPARENT ANTIBIOTIC DICTIONARY

   ReferenceValue is the SICdb medication name.

   The list includes common international generic names and several European/
   Austrian proprietary names frequently encountered in hospital medication
   systems.

   It intentionally excludes antifungals/antivirals from the bacterial
   antibiotic proxy.
   ========================================================================= */

CREATE TABLE study_make30_final.sic_antibiotic_patterns
(
    pattern_id      integer PRIMARY KEY,
    pattern_regex   text NOT NULL,
    pattern_label   text NOT NULL
);

INSERT INTO study_make30_final.sic_antibiotic_patterns
(pattern_id, pattern_regex, pattern_label)
VALUES
( 1, 'vancomycin|vancocin',                              'vancomycin'),
( 2, 'teicoplanin|targocid',                             'teicoplanin'),
( 3, 'piperacillin|tazobactam|tazobac|zosyn',           'piperacillin-tazobactam'),
( 4, 'ampicillin|sulbactam|unasyn',                      'ampicillin-sulbactam'),
( 5, 'amoxicillin|clavulan|augmentin',                   'amoxicillin-clavulanate'),
( 6, 'benzylpenicillin|penicillin[ -]?g|penicillin',     'penicillin'),
( 7, 'flucloxacillin|oxacillin|nafcillin',               'anti-staphylococcal penicillin'),
( 8, 'cefazolin|ancef',                                  'cefazolin'),
( 9, 'cefuroxime|zinacef|ceftin',                        'cefuroxime'),
(10, 'cefoxitin|mefoxin',                                'cefoxitin'),
(11, 'cefotaxime|claforan',                              'cefotaxime'),
(12, 'ceftriaxone|rocephin',                             'ceftriaxone'),
(13, 'ceftazidime|fortum',                               'ceftazidime'),
(14, 'cefepime|maxipime',                                'cefepime'),
(15, 'ceftaroline|zinforo',                              'ceftaroline'),
(16, 'meropenem|meronem',                                'meropenem'),
(17, 'imipenem|cilastatin|tienam|primaxin',              'imipenem-cilastatin'),
(18, 'ertapenem|invanz',                                 'ertapenem'),
(19, 'aztreonam|azactam',                                'aztreonam'),
(20, 'ciprofloxacin|ciproxin|cipro',                     'ciprofloxacin'),
(21, 'levofloxacin|tavanic|levaquin',                    'levofloxacin'),
(22, 'moxifloxacin|avalox|avelox',                       'moxifloxacin'),
(23, 'azithromycin|zithromax',                           'azithromycin'),
(24, 'clarithromycin|klacid|biaxin',                     'clarithromycin'),
(25, 'erythromycin',                                     'erythromycin'),
(26, 'clindamycin|dalacin|cleocin',                      'clindamycin'),
(27, 'metronidazole|flagyl',                             'metronidazole'),
(28, 'linezolid|zyvox|zyvoxid',                         'linezolid'),
(29, 'daptomycin|cubicin',                               'daptomycin'),
(30, 'gentamicin',                                       'gentamicin'),
(31, 'tobramycin',                                       'tobramycin'),
(32, 'amikacin',                                         'amikacin'),
(33, 'trimethoprim|sulfamethox|cotrimox|co-trimox',      'trimethoprim-sulfamethoxazole'),
(34, 'fosfomycin',                                       'fosfomycin'),
(35, 'tigecycline|tygacil',                              'tigecycline'),
(36, 'colistin|colistimethate|polymyxin',                'polymyxin');


/* =========================================================================
   8. RECOGNIZED ANTIBIOTIC ADMINISTRATION EVENTS

   SICdb medication.Offset:
       documented start time of drug administration in seconds.

   ORAKLE SICdb rule:
       occurrence of antibiotic administration = suspected infection proxy.

   No >=2-dose requirement.
   No additional route requirement is imposed because route is not a SICdb
   medication field and ORAKLE's SICdb-specific public wording is occurrence
   of antibiotic administrations.
   ========================================================================= */

CREATE TABLE study_make30_final.sic_antibiotic_admin AS

SELECT DISTINCT
    m.case_id,

    m.offset_sec
        AS antibiotic_offset_sec,

    m.drug_id,

    r.ref_value
        AS drug_name

FROM study_make30_final.sic_src_medication m

INNER JOIN study_make30_final.sic_cases c
    ON m.case_id = c.case_id

INNER JOIN study_make30_final.sic_src_refs r
    ON m.drug_id = r.ref_id

WHERE m.offset_sec IS NOT NULL

  AND EXISTS
  (
      SELECT 1

      FROM study_make30_final.sic_antibiotic_patterns p

      WHERE LOWER(COALESCE(r.ref_value,'')) ~
            p.pattern_regex
  );


CREATE INDEX idx_sic_abx_case_time
ON study_make30_final.sic_antibiotic_admin
(case_id, antibiotic_offset_sec);


/* =========================================================================
   9. SUSPECTED INFECTION ANCHORS

   Every unique recognized antibiotic administration time is a candidate.

   antibiotic_admin_count_24h:
       descriptive only;
       NOT an inclusion criterion.
   ========================================================================= */

CREATE TABLE study_make30_final.sic_suspected_infection AS

SELECT
    a.case_id,

    a.antibiotic_offset_sec
        AS suspected_infection_offset_sec,

    (
        SELECT COUNT(*)::integer

        FROM study_make30_final.sic_antibiotic_admin b

        WHERE b.case_id = a.case_id

          AND b.antibiotic_offset_sec
              >= a.antibiotic_offset_sec

          AND b.antibiotic_offset_sec
              <= a.antibiotic_offset_sec + 86400
    ) AS antibiotic_admin_count_24h

FROM
(
    SELECT DISTINCT
        case_id,
        antibiotic_offset_sec

    FROM study_make30_final.sic_antibiotic_admin
) a;


CREATE UNIQUE INDEX idx_sic_si_case_time
ON study_make30_final.sic_suspected_infection
(case_id, suspected_infection_offset_sec);


/* =========================================================================
   10. SICdb SOFA >=2

   Uses the v1.0.7+ SICdb SOFA signal:
       DataID 3139

   If duplicate rows occur in the same case/hour, maximum SOFA is retained.
   ========================================================================= */

CREATE TABLE study_make30_final.sic_sofa_ge2 AS

SELECT
    f.case_id,

    f.offset_sec
        AS sofa_offset_sec,

    MAX(f.val)::double precision
        AS sofa_score

FROM study_make30_final.sic_src_data_float_h f

INNER JOIN study_make30_final.sic_cases c
    ON f.case_id = c.case_id

WHERE f.data_id = 3139
  AND f.val IS NOT NULL
  AND f.val >= 2

GROUP BY
    f.case_id,
    f.offset_sec;


CREATE UNIQUE INDEX idx_sic_sofa_case_time
ON study_make30_final.sic_sofa_ge2
(case_id, sofa_offset_sec);


/* =========================================================================
   11. ORAKLE-STYLE ANTIBIOTIC + SOFA MATCHING

   SOFA window:
       antibiotic -24 h  through antibiotic +12 h

   seconds:
       -86400
       +43200

   sepsis_offset:
       earlier of antibiotic anchor and qualifying SOFA.
   ========================================================================= */

CREATE TABLE study_make30_final.sic_sepsis_pairs AS

SELECT
    c.patient_id,
    c.case_id,

    si.suspected_infection_offset_sec,
    si.antibiotic_admin_count_24h,

    s.sofa_offset_sec,
    s.sofa_score,

    LEAST(
        si.suspected_infection_offset_sec,
        s.sofa_offset_sec
    )::bigint AS sepsis_offset_sec

FROM study_make30_final.sic_cases c

INNER JOIN study_make30_final.sic_suspected_infection si
    ON c.case_id = si.case_id

INNER JOIN study_make30_final.sic_sofa_ge2 s
    ON c.case_id = s.case_id

   AND s.sofa_offset_sec
       >= si.suspected_infection_offset_sec - 86400

   AND s.sofa_offset_sec
       <= si.suspected_infection_offset_sec + 43200;


CREATE INDEX idx_sic_sepsis_pairs_case
ON study_make30_final.sic_sepsis_pairs
(case_id, sepsis_offset_sec);


/* =========================================================================
   12. FIRST QUALIFYING SEPSIS EVENT PER CASE
   ========================================================================= */

CREATE TABLE study_make30_final.sic_sepsis_first_case AS

SELECT
    patient_id,
    case_id,

    suspected_infection_offset_sec,
    antibiotic_admin_count_24h,

    sofa_offset_sec,
    sofa_score,

    sepsis_offset_sec

FROM
(
    SELECT
        p.*,

        ROW_NUMBER() OVER
        (
            PARTITION BY p.case_id

            ORDER BY
                p.sepsis_offset_sec,
                p.suspected_infection_offset_sec,
                p.sofa_offset_sec
        ) AS rn_case

    FROM study_make30_final.sic_sepsis_pairs p
) q

WHERE rn_case = 1;


/* =========================================================================
   13. ADULT SEPSIS OCCURRING AFTER ACTUAL ICU ADMISSION

   SICdb can include preceding surgery before ICUOffset.

   Require:
       age >=18
       t0 >= ICUOffset
       t0 < TimeOfStay
   ========================================================================= */

CREATE TABLE study_make30_final.sic_sepsis_adult AS

SELECT
    s.patient_id,
    s.case_id,

    c.age_years,

    c.icu_offset_sec,
    c.time_of_stay_sec,

    c.offset_after_first_admission_sec,

    c.offset_of_death_sec,
    c.hospital_discharge_day,
    c.hospital_stay_days,

    s.suspected_infection_offset_sec,
    s.antibiotic_admin_count_24h,

    s.sofa_offset_sec,
    s.sofa_score,

    s.sepsis_offset_sec,

    (s.sepsis_offset_sec + 86400)::bigint
        AS landmark_offset_sec,

    (s.sepsis_offset_sec + 2592000)::bigint
        AS sepsis_day30_offset_sec,

    (
        c.offset_after_first_admission_sec
        + s.sepsis_offset_sec
    )::bigint AS absolute_patient_sepsis_offset_sec

FROM study_make30_final.sic_sepsis_first_case s

INNER JOIN study_make30_final.sic_cases c
    ON s.case_id = c.case_id

WHERE c.age_years >= 18

  AND s.sepsis_offset_sec >= c.icu_offset_sec

  AND s.sepsis_offset_sec < c.time_of_stay_sec;


/* =========================================================================
   14. FIRST QUALIFYING SEPSIS CASE PER PATIENT

   Readmissions are ordered using:
       OffsetAfterFirstAdmission + within-case sepsis offset
   rather than comparing within-case offsets alone.
   ========================================================================= */

CREATE TABLE study_make30_final.sic_sepsis_index AS

SELECT
    patient_id,
    case_id,

    age_years,

    icu_offset_sec,
    time_of_stay_sec,

    offset_after_first_admission_sec,

    offset_of_death_sec,
    hospital_discharge_day,
    hospital_stay_days,

    suspected_infection_offset_sec,
    antibiotic_admin_count_24h,

    sofa_offset_sec,
    sofa_score,

    sepsis_offset_sec,
    landmark_offset_sec,
    sepsis_day30_offset_sec,

    absolute_patient_sepsis_offset_sec

FROM
(
    SELECT
        a.*,

        ROW_NUMBER() OVER
        (
            PARTITION BY a.patient_id

            ORDER BY
                a.absolute_patient_sepsis_offset_sec,
                a.case_id
        ) AS rn_patient

    FROM study_make30_final.sic_sepsis_adult a
) q

WHERE rn_patient = 1;


CREATE UNIQUE INDEX idx_sic_sepsis_index_patient
ON study_make30_final.sic_sepsis_index(patient_id);

CREATE UNIQUE INDEX idx_sic_sepsis_index_case
ON study_make30_final.sic_sepsis_index(case_id);


/* =========================================================================
   15. +24-H LANDMARK

   Require the SICdb case to remain observable through t0+24 h:
       TimeOfStay >= landmark.

   Hospital-discharge-day is day-level only and is not used to artificially
   create a more precise ICU timestamp.
   ========================================================================= */

CREATE TABLE study_make30_final.sic_landmark_cohort AS

SELECT
    s.*,

    CASE
        WHEN s.time_of_stay_sec
             >= s.landmark_offset_sec
        THEN 1
        ELSE 0
    END AS reached_24h_landmark

FROM study_make30_final.sic_sepsis_index s;


CREATE UNIQUE INDEX idx_sic_landmark_case
ON study_make30_final.sic_landmark_cohort(case_id);


/* =========================================================================
   16. EARLY KRT/CRRT DURING PHENOTYPE WINDOW [t0, t0+24 h)

   Signal evidence:
       data_float_h with KRT/CRRT references and a non-zero value

   Range evidence:
       data_range interval with KRT/CRRT reference overlapping the window

   This is the phenotype-contamination exclusion.
   Post-landmark "new KRT" for MAKE will be constructed separately later.
   ========================================================================= */

DROP TABLE IF EXISTS tmp_sic_early_krt;

CREATE TEMP TABLE tmp_sic_early_krt AS

WITH float_krt AS
(
    SELECT DISTINCT
        c.case_id

    FROM study_make30_final.sic_landmark_cohort c

    INNER JOIN study_make30_final.sic_src_data_float_h f
        ON c.case_id = f.case_id

       AND f.offset_sec >= c.sepsis_offset_sec
       AND f.offset_sec <  c.landmark_offset_sec

    INNER JOIN study_make30_final.sic_krt_ids k
        ON f.data_id = k.data_id

    WHERE c.reached_24h_landmark = 1

      AND f.val IS NOT NULL
      AND ABS(f.val) > 0
),

range_krt AS
(
    SELECT DISTINCT
        c.case_id

    FROM study_make30_final.sic_landmark_cohort c

    INNER JOIN study_make30_final.sic_src_data_range r
        ON c.case_id = r.case_id

    INNER JOIN study_make30_final.sic_krt_ids k
        ON r.data_id = k.data_id

    WHERE c.reached_24h_landmark = 1

      AND COALESCE(r.offset_end_sec, r.offset_sec)
          > c.sepsis_offset_sec

      AND r.offset_sec
          < c.landmark_offset_sec
),

all_krt AS
(
    SELECT case_id
    FROM float_krt

    UNION

    SELECT case_id
    FROM range_krt
)

SELECT
    c.case_id,

    CASE
        WHEN k.case_id IS NOT NULL
        THEN 1
        ELSE 0
    END AS krt_0_24h

FROM study_make30_final.sic_landmark_cohort c

LEFT JOIN all_krt k
    ON c.case_id = k.case_id

WHERE c.reached_24h_landmark = 1;


CREATE UNIQUE INDEX idx_tmp_sic_early_krt
ON tmp_sic_early_krt(case_id);


/* =========================================================================
   17. TRAJECTORY POOL
   ========================================================================= */

CREATE TABLE study_make30_final.sic_trajectory_pool AS

SELECT
    c.*,

    COALESCE(k.krt_0_24h,0)
        AS krt_0_24h

FROM study_make30_final.sic_landmark_cohort c

LEFT JOIN tmp_sic_early_krt k
    ON c.case_id = k.case_id

WHERE c.reached_24h_landmark = 1

  AND COALESCE(k.krt_0_24h,0) = 0;


CREATE UNIQUE INDEX idx_sic_trajcand_case
ON study_make30_final.sic_trajectory_pool(case_id);


/* =========================================================================
   18. SICdb HOURLY HR / MAP

   IMPORTANT:
       data_float_h.Val is SICdb's database-provided one-hour aggregate.

   HR priority:
       707 ECG
       708 SpO2
       724 arterial-pressure-derived HR

   MAP sources:
       703 invasive arterial MAP
       706 non-invasive MAP

   MAP preprocessing / source selection:
       both invasive MAP and NIBP <=20 mmHg -> exclude before source selection
       if both valid hourly sources exist and |invasive-NIBP| >=30 mmHg -> use NIBP
       otherwise valid invasive MAP has priority

   MSI-input plausibility:
       HR  >0 and <300
       MAP >20 and <300 for both invasive MAP and NIBP
   ========================================================================= */

DROP TABLE IF EXISTS tmp_sic_vital_hourly;

CREATE TEMP TABLE tmp_sic_vital_hourly AS

SELECT
    b.case_id,

    FLOOR(
        (f.offset_sec - b.sepsis_offset_sec)
        / 3600.0
    )::integer AS hour_bin,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY f.val)
    FILTER
    (
        WHERE f.data_id = 707
          AND f.val > 0
          AND f.val < 300
    ) AS hr_ecg,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY f.val)
    FILTER
    (
        WHERE f.data_id = 708
          AND f.val > 0
          AND f.val < 300
    ) AS hr_spo2,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY f.val)
    FILTER
    (
        WHERE f.data_id = 724
          AND f.val > 0
          AND f.val < 300
    ) AS hr_abp,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY f.val)
    FILTER
    (
        WHERE f.data_id = 703
          AND f.val > 20
          AND f.val < 300
    ) AS map_invasive,

    PERCENTILE_CONT(0.5)
    WITHIN GROUP (ORDER BY f.val)
    FILTER
    (
        WHERE f.data_id = 706
          AND f.val > 20
          AND f.val < 300
    ) AS map_noninvasive

FROM study_make30_final.sic_trajectory_pool b

INNER JOIN study_make30_final.sic_src_data_float_h f
    ON b.case_id = f.case_id

   AND f.offset_sec >= b.sepsis_offset_sec
   AND f.offset_sec <  b.landmark_offset_sec

WHERE f.data_id IN
(
    707,708,724,
    703,706
)

GROUP BY
    b.case_id,

    FLOOR(
        (f.offset_sec - b.sepsis_offset_sec)
        / 3600.0
    )::integer;


CREATE INDEX idx_tmp_sic_vital_hourly
ON tmp_sic_vital_hourly(case_id, hour_bin);


/* =========================================================================
   19. VALID HOURLY MSI — UNIFIED MAP QC (FINAL <=20 RULE)

   HR source hierarchy:
      ECG HR > SpO2 HR > ABP-derived HR.

   MAP preprocessing and selection order:
      1) SICdb data_float_h provides database-native hourly MAP aggregates;
      2) discard every invasive MAP value <=20 mmHg and every NIBP MAP
         value <=20 mmHg before source selection;
      3) if only one valid hourly MAP source exists, use that source;
      4) if both valid hourly sources exist:
           - if |invasive MAP - NIBP| >=30 mmHg, use NIBP;
           - otherwise use invasive MAP;
      5) if neither valid source exists, MAP is missing.

   No minimum raw-measurement count is imposed within an hour.

   NOTE:
      The <=20-mmHg exclusion is applied symmetrically to invasive MAP and
      NIBP, before source-priority and discordance rules are evaluated.
   ========================================================================= */

CREATE TABLE study_make30_final.sic_msi_hourly AS

WITH selected AS
(
    SELECT
        case_id,
        hour_bin,

        hr_ecg::double precision AS hr_ecg,
        hr_spo2::double precision AS hr_spo2,
        hr_abp::double precision AS hr_abp,

        COALESCE(
            hr_ecg,
            hr_spo2,
            hr_abp
        )::double precision AS hr_used,

        CASE
            WHEN hr_ecg IS NOT NULL
            THEN 'ECG'

            WHEN hr_spo2 IS NOT NULL
            THEN 'SpO2'

            WHEN hr_abp IS NOT NULL
            THEN 'ABP'

            ELSE NULL
        END::text AS hr_source,

        map_invasive::double precision AS map_invasive,
        map_noninvasive::double precision AS map_noninvasive,

        CASE
            WHEN map_invasive IS NOT NULL
             AND map_noninvasive IS NOT NULL
            THEN ABS(map_invasive - map_noninvasive)
            ELSE NULL
        END::double precision AS map_abs_diff,

        CASE
            WHEN map_invasive IS NOT NULL
             AND map_noninvasive IS NOT NULL
             AND ABS(map_invasive - map_noninvasive) >= 30
            THEN map_noninvasive

            WHEN map_invasive IS NOT NULL
            THEN map_invasive

            WHEN map_noninvasive IS NOT NULL
            THEN map_noninvasive

            ELSE NULL
        END::double precision AS map_used,

        CASE
            WHEN map_invasive IS NOT NULL
             AND map_noninvasive IS NOT NULL
             AND ABS(map_invasive - map_noninvasive) >= 30
            THEN 'noninvasive_fallback_discordant_ge30'

            WHEN map_invasive IS NOT NULL
            THEN 'invasive'

            WHEN map_noninvasive IS NOT NULL
            THEN 'noninvasive'

            ELSE NULL
        END::text AS map_source

    FROM tmp_sic_vital_hourly

    WHERE hour_bin BETWEEN 0 AND 23
)

SELECT
    case_id,
    hour_bin,
    hr_ecg,
    hr_spo2,
    hr_abp,
    hr_used,
    hr_source,
    map_invasive,
    map_noninvasive,
    map_abs_diff,
    map_used,
    map_source,

    CASE
        WHEN hr_used IS NOT NULL
         AND map_used IS NOT NULL
         AND map_used > 0
        THEN hr_used / map_used
        ELSE NULL
    END::double precision AS msi

FROM selected;

CREATE UNIQUE INDEX idx_sic_msi_hourly
ON study_make30_final.sic_msi_hourly(case_id, hour_bin);


/* =========================================================================
   20. HOURLY MSI -> TWELVE 2-H MSI WINDOWS
   ========================================================================= */
CREATE TABLE study_make30_final.sic_msi_2h AS
SELECT case_id,FLOOR(hour_bin/2.0)::integer time_bin,MAX(msi)::double precision msi,COUNT(msi)::integer valid_hourly_msi_count
FROM study_make30_final.sic_msi_hourly WHERE msi IS NOT NULL GROUP BY case_id,FLOOR(hour_bin/2.0)::integer;
CREATE UNIQUE INDEX idx_sic_msi_2h ON study_make30_final.sic_msi_2h(case_id,time_bin);

/* =========================================================================
   21. CORE-4 ELECTROLYTES

   SOURCE-PRIORITY RULE
   --------------------
   For every case + variable + 6-h bin:

       1) median of CENTRAL-lab measurements, if present;
       2) otherwise median of blood-gas fallback measurements.

   Central and BGA values are not pooled when central values are available.

   Electrolyte hard-QC is applied to each raw lab measurement before
   source-specific medians are calculated:
       sodium >300 mmol/L -> discarded (300 retained);
       chloride <10 mmol/L -> discarded (10 retained).
   ========================================================================= */

CREATE TABLE study_make30_final.sic_electrolyte_6h AS

WITH lab_raw AS
(
    SELECT
        b.case_id,

        FLOOR(
            (l.offset_sec - b.sepsis_offset_sec)
            / 21600.0
        )::integer AS time_bin,

        m.variable_name,
        m.source_priority,

        l.lab_value

    FROM study_make30_final.sic_trajectory_pool b

    INNER JOIN study_make30_final.sic_src_laboratory l
        ON b.case_id = l.case_id

       AND l.offset_sec >= b.sepsis_offset_sec
       AND l.offset_sec <  b.landmark_offset_sec

    INNER JOIN study_make30_final.sic_lab_mapping m
        ON l.lab_id = m.lab_id

    WHERE l.lab_value IS NOT NULL
      AND NOT (m.variable_name = 'sodium'   AND l.lab_value > 300)
      AND NOT (m.variable_name = 'chloride' AND l.lab_value < 10)
),

by_source AS
(
    SELECT
        case_id,
        time_bin,
        variable_name,
        source_priority,

        PERCENTILE_CONT(0.5)
        WITHIN GROUP (ORDER BY lab_value)
            AS source_median,

        COUNT(*)::integer
            AS source_measurement_count

    FROM lab_raw

    GROUP BY
        case_id,
        time_bin,
        variable_name,
        source_priority
)

SELECT
    case_id,
    time_bin,
    variable_name,

    COALESCE
    (
        MAX(source_median)
            FILTER (WHERE source_priority = 1),

        MAX(source_median)
            FILTER (WHERE source_priority = 2)
    )::double precision AS median_value,

    CASE
        WHEN MAX(source_median)
             FILTER (WHERE source_priority = 1)
             IS NOT NULL

        THEN 1

        WHEN MAX(source_median)
             FILTER (WHERE source_priority = 2)
             IS NOT NULL

        THEN 2

        ELSE NULL
    END::integer AS source_priority_used,

    COALESCE
    (
        MAX(source_measurement_count)
            FILTER (WHERE source_priority = 1),

        MAX(source_measurement_count)
            FILTER (WHERE source_priority = 2),

        0
    )::integer AS raw_measurement_count

FROM by_source

GROUP BY
    case_id,
    time_bin,
    variable_name;


CREATE UNIQUE INDEX idx_sic_electrolyte_6h
ON study_make30_final.sic_electrolyte_6h
(case_id, time_bin, variable_name);



/* =========================================================================
   22. MSI 2-H TRAJECTORY LONG TABLE — 12 ROWS PER TRAJECTORY-POOL PATIENT
   ========================================================================= */
CREATE TABLE study_make30_final.sic_msi_trajectory_long AS
WITH patient_bins AS
(
 SELECT b.patient_id,b.case_id,b.sepsis_offset_sec,b.landmark_offset_sec,gs.time_bin FROM study_make30_final.sic_trajectory_pool b CROSS JOIN generate_series(0,11) gs(time_bin)
)
SELECT p.patient_id,p.case_id,p.sepsis_offset_sec,p.landmark_offset_sec,p.time_bin,('M'||(p.time_bin+1)::text)::text time_label,(p.time_bin*2)::integer window_start_hour,((p.time_bin+1)*2)::integer window_end_hour,(p.time_bin*2+1)::integer time_midpoint_hour,m.msi,m.valid_hourly_msi_count
FROM patient_bins p LEFT JOIN study_make30_final.sic_msi_2h m ON p.case_id=m.case_id AND p.time_bin=m.time_bin ORDER BY p.case_id,p.time_bin;
CREATE UNIQUE INDEX idx_sic_msi_trajectory_long ON study_make30_final.sic_msi_trajectory_long(case_id,time_bin);

/* =========================================================================
   23. ELECTROLYTE 6-H TRAJECTORY LONG TABLE — 4 ROWS PER TRAJECTORY-POOL PATIENT
   ========================================================================= */
CREATE TABLE study_make30_final.sic_electrolyte_trajectory_long AS
WITH patient_bins AS
(
 SELECT b.patient_id,b.case_id,b.sepsis_offset_sec,b.landmark_offset_sec,gs.time_bin FROM study_make30_final.sic_trajectory_pool b CROSS JOIN generate_series(0,3) gs(time_bin)
),ew AS
(
 SELECT case_id,time_bin,MAX(median_value) FILTER(WHERE variable_name='sodium') sodium,MAX(median_value) FILTER(WHERE variable_name='potassium') potassium,MAX(median_value) FILTER(WHERE variable_name='chloride') chloride,MAX(median_value) FILTER(WHERE variable_name='bicarbonate') bicarbonate
 FROM study_make30_final.sic_electrolyte_6h GROUP BY case_id,time_bin
)
SELECT p.patient_id,p.case_id,p.sepsis_offset_sec,p.landmark_offset_sec,p.time_bin,('T'||(p.time_bin+1)::text)::text time_label,(p.time_bin*6)::integer window_start_hour,((p.time_bin+1)*6)::integer window_end_hour,(p.time_bin*6+3)::integer time_midpoint_hour,e.sodium,e.potassium,e.chloride,e.bicarbonate
FROM patient_bins p LEFT JOIN ew e ON p.case_id=e.case_id AND p.time_bin=e.time_bin ORDER BY p.case_id,p.time_bin;
CREATE UNIQUE INDEX idx_sic_electrolyte_trajectory_long ON study_make30_final.sic_electrolyte_trajectory_long(case_id,time_bin);

/* =========================================================================
   24. PATIENT COVERAGE + FINAL PRIMARY ELIGIBILITY
   ========================================================================= */
CREATE TABLE study_make30_final.sic_trajectory_coverage AS
WITH m AS (SELECT patient_id,case_id,COUNT(msi)::integer msi_bins FROM study_make30_final.sic_msi_trajectory_long GROUP BY patient_id,case_id),
e AS (SELECT patient_id,case_id,COUNT(sodium)::integer sodium_bins,COUNT(potassium)::integer potassium_bins,COUNT(chloride)::integer chloride_bins,COUNT(bicarbonate)::integer bicarbonate_bins FROM study_make30_final.sic_electrolyte_trajectory_long GROUP BY patient_id,case_id)
SELECT m.*,e.sodium_bins,e.potassium_bins,e.chloride_bins,e.bicarbonate_bins,
CASE WHEN m.msi_bins>=6 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END core4_eligible_msi6,
CASE WHEN m.msi_bins>=9 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END core4_eligible_msi9,
CASE WHEN m.msi_bins=12 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END core4_eligible_msi12,
CASE WHEN m.msi_bins>=9 AND e.sodium_bins>=2 AND e.potassium_bins>=2 AND e.chloride_bins>=2 AND e.bicarbonate_bins>=2 THEN 1 ELSE 0 END core4_trajectory_eligible
FROM m INNER JOIN e USING(patient_id,case_id);
CREATE UNIQUE INDEX idx_sic_trajectory_coverage ON study_make30_final.sic_trajectory_coverage(case_id);

CREATE TABLE study_make30_final.sic_core4_trajectory_cohort AS
SELECT c.*,v.msi_bins,v.sodium_bins,v.potassium_bins,v.chloride_bins,v.bicarbonate_bins,v.core4_eligible_msi6,v.core4_eligible_msi9,v.core4_eligible_msi12,v.core4_trajectory_eligible
FROM study_make30_final.sic_trajectory_pool c INNER JOIN study_make30_final.sic_trajectory_coverage v USING(case_id) WHERE v.core4_trajectory_eligible=1;
CREATE UNIQUE INDEX idx_sic_core4_cohort_case ON study_make30_final.sic_core4_trajectory_cohort(case_id);

CREATE TABLE study_make30_final.sic_core4_msi_trajectory_long AS SELECT t.* FROM study_make30_final.sic_msi_trajectory_long t INNER JOIN study_make30_final.sic_trajectory_coverage c USING(case_id) WHERE c.core4_trajectory_eligible=1 ORDER BY t.case_id,t.time_bin;
CREATE UNIQUE INDEX idx_sic_core4_msi_long ON study_make30_final.sic_core4_msi_trajectory_long(case_id,time_bin);
CREATE TABLE study_make30_final.sic_core4_electrolyte_trajectory_long AS SELECT t.* FROM study_make30_final.sic_electrolyte_trajectory_long t INNER JOIN study_make30_final.sic_trajectory_coverage c USING(case_id) WHERE c.core4_trajectory_eligible=1 ORDER BY t.case_id,t.time_bin;
CREATE UNIQUE INDEX idx_sic_core4_electrolyte_long ON study_make30_final.sic_core4_electrolyte_trajectory_long(case_id,time_bin);

CREATE TABLE study_make30_final.sic_core4_multitrajectory_wide AS
WITH mw AS
(
 SELECT case_id,MAX(msi) FILTER(WHERE time_bin=0) msi_m1,MAX(msi) FILTER(WHERE time_bin=1) msi_m2,MAX(msi) FILTER(WHERE time_bin=2) msi_m3,MAX(msi) FILTER(WHERE time_bin=3) msi_m4,MAX(msi) FILTER(WHERE time_bin=4) msi_m5,MAX(msi) FILTER(WHERE time_bin=5) msi_m6,MAX(msi) FILTER(WHERE time_bin=6) msi_m7,MAX(msi) FILTER(WHERE time_bin=7) msi_m8,MAX(msi) FILTER(WHERE time_bin=8) msi_m9,MAX(msi) FILTER(WHERE time_bin=9) msi_m10,MAX(msi) FILTER(WHERE time_bin=10) msi_m11,MAX(msi) FILTER(WHERE time_bin=11) msi_m12 FROM study_make30_final.sic_core4_msi_trajectory_long GROUP BY case_id
),ew AS
(
 SELECT case_id,MAX(sodium) FILTER(WHERE time_bin=0) sodium_t1,MAX(sodium) FILTER(WHERE time_bin=1) sodium_t2,MAX(sodium) FILTER(WHERE time_bin=2) sodium_t3,MAX(sodium) FILTER(WHERE time_bin=3) sodium_t4,MAX(potassium) FILTER(WHERE time_bin=0) potassium_t1,MAX(potassium) FILTER(WHERE time_bin=1) potassium_t2,MAX(potassium) FILTER(WHERE time_bin=2) potassium_t3,MAX(potassium) FILTER(WHERE time_bin=3) potassium_t4,MAX(chloride) FILTER(WHERE time_bin=0) chloride_t1,MAX(chloride) FILTER(WHERE time_bin=1) chloride_t2,MAX(chloride) FILTER(WHERE time_bin=2) chloride_t3,MAX(chloride) FILTER(WHERE time_bin=3) chloride_t4,MAX(bicarbonate) FILTER(WHERE time_bin=0) bicarbonate_t1,MAX(bicarbonate) FILTER(WHERE time_bin=1) bicarbonate_t2,MAX(bicarbonate) FILTER(WHERE time_bin=2) bicarbonate_t3,MAX(bicarbonate) FILTER(WHERE time_bin=3) bicarbonate_t4 FROM study_make30_final.sic_core4_electrolyte_trajectory_long GROUP BY case_id
)
SELECT c.patient_id,c.case_id,c.sepsis_offset_sec,c.landmark_offset_sec,
       mw.msi_m1,mw.msi_m2,mw.msi_m3,mw.msi_m4,mw.msi_m5,mw.msi_m6,mw.msi_m7,mw.msi_m8,mw.msi_m9,mw.msi_m10,mw.msi_m11,mw.msi_m12,
       ew.sodium_t1,ew.sodium_t2,ew.sodium_t3,ew.sodium_t4,
       ew.potassium_t1,ew.potassium_t2,ew.potassium_t3,ew.potassium_t4,
       ew.chloride_t1,ew.chloride_t2,ew.chloride_t3,ew.chloride_t4,
       ew.bicarbonate_t1,ew.bicarbonate_t2,ew.bicarbonate_t3,ew.bicarbonate_t4
FROM study_make30_final.sic_core4_trajectory_cohort c LEFT JOIN mw USING(case_id) LEFT JOIN ew USING(case_id);
CREATE UNIQUE INDEX idx_sic_core4_multitrajectory_wide ON study_make30_final.sic_core4_multitrajectory_wide(case_id);

CREATE TABLE study_make30_final.sic_msi_gt4_trace AS
SELECT c.patient_id,h.case_id,h.hour_bin,FLOOR(h.hour_bin/2.0)::integer msi_2h_bin,('M'||(FLOOR(h.hour_bin/2.0)::integer+1)::text)::text msi_2h_label,
       h.hr_ecg,h.hr_spo2,h.hr_abp,h.hr_used,h.hr_source,h.map_invasive,h.map_noninvasive,h.map_abs_diff,h.map_used,h.map_source,h.msi hourly_msi,m.msi msi_2h_max
FROM study_make30_final.sic_msi_hourly h INNER JOIN study_make30_final.sic_trajectory_pool c USING(case_id)
LEFT JOIN study_make30_final.sic_msi_2h m ON h.case_id=m.case_id AND FLOOR(h.hour_bin/2.0)::integer=m.time_bin
WHERE h.msi>4 ORDER BY h.msi DESC,h.case_id,h.hour_bin;
CREATE INDEX idx_sic_msi_gt4_trace ON study_make30_final.sic_msi_gt4_trace(case_id,hour_bin);



/* =========================================================================
   25. ORAKLE MAKE30 + BASELINE / ADJUSTMENT COVARIATES

   Exposure: Sepsis t0 through t0+24 h (frozen trajectory definition).
   Outcome risk begins strictly after t0+24 h.
   MAKE30 components/censoring follow ORAKLE.
   ========================================================================= */

/* Extra case fields are obtained by semantic column matching; critical IDs
   and time columns continue to come from the already validated base layer. */
CREATE TABLE study_make30_final.sic_case_extra
(
    case_id bigint PRIMARY KEY,
    sex_ref_id integer,
    weight_kg double precision,
    height_cm double precision,
    saps3 double precision,
    hospital_discharge_type_ref integer,
    icd10_main text,
    icd10_main_text text
);

DO $$
DECLARE
    v_schema text; t_cases text; c_caseid text; c_sex text; c_wt text; c_ht text; c_saps text; c_hdt text; c_icd text; c_icdt text;
BEGIN
    SELECT n.nspname INTO v_schema
    FROM pg_namespace n
    WHERE EXISTS (SELECT 1 FROM pg_class c WHERE c.relnamespace=n.oid AND LOWER(c.relname)='cases' AND c.relkind IN('r','p','v','m'))
      AND EXISTS (SELECT 1 FROM pg_class c WHERE c.relnamespace=n.oid AND LOWER(c.relname)='d_references' AND c.relkind IN('r','p','v','m'))
      AND EXISTS (SELECT 1 FROM pg_class c WHERE c.relnamespace=n.oid AND LOWER(c.relname)='laboratory' AND c.relkind IN('r','p','v','m'))
    ORDER BY CASE WHEN n.nspname='public' THEN 0 ELSE 1 END,n.nspname LIMIT 1;
    SELECT table_name INTO t_cases FROM information_schema.tables WHERE table_schema=v_schema AND LOWER(table_name)='cases' LIMIT 1;
    SELECT column_name INTO c_caseid FROM information_schema.columns WHERE table_schema=v_schema AND table_name=t_cases AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='caseid' LIMIT 1;
    SELECT column_name INTO c_sex FROM information_schema.columns WHERE table_schema=v_schema AND table_name=t_cases AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='sex' LIMIT 1;
    SELECT column_name INTO c_wt FROM information_schema.columns WHERE table_schema=v_schema AND table_name=t_cases AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='weightonadmission' LIMIT 1;
    SELECT column_name INTO c_ht FROM information_schema.columns WHERE table_schema=v_schema AND table_name=t_cases AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='heightonadmission' LIMIT 1;
    SELECT column_name INTO c_saps FROM information_schema.columns WHERE table_schema=v_schema AND table_name=t_cases AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='saps3' LIMIT 1;
    SELECT column_name INTO c_hdt FROM information_schema.columns WHERE table_schema=v_schema AND table_name=t_cases AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='hospitaldischargetype' LIMIT 1;
    SELECT column_name INTO c_icd FROM information_schema.columns WHERE table_schema=v_schema AND table_name=t_cases AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='icd10main' LIMIT 1;
    SELECT column_name INTO c_icdt FROM information_schema.columns WHERE table_schema=v_schema AND table_name=t_cases AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='icd10maintext' LIMIT 1;
    IF c_caseid IS NULL THEN RAISE EXCEPTION 'SICdb cases.CaseID not found for final covariate extraction.'; END IF;
    EXECUTE FORMAT(
        'INSERT INTO study_make30_final.sic_case_extra(case_id,sex_ref_id,weight_kg,height_cm,saps3,hospital_discharge_type_ref,icd10_main,icd10_main_text) SELECT %I::bigint,%s,%s,%s,%s,%s,%s,%s FROM %I.%I',
        c_caseid,
        CASE WHEN c_sex IS NULL THEN 'NULL::integer' ELSE FORMAT('%I::integer',c_sex) END,
        CASE WHEN c_wt IS NULL THEN 'NULL::double precision' ELSE FORMAT('%I::double precision',c_wt) END,
        CASE WHEN c_ht IS NULL THEN 'NULL::double precision' ELSE FORMAT('%I::double precision',c_ht) END,
        CASE WHEN c_saps IS NULL THEN 'NULL::double precision' ELSE FORMAT('%I::double precision',c_saps) END,
        CASE WHEN c_hdt IS NULL THEN 'NULL::integer' ELSE FORMAT('%I::integer',c_hdt) END,
        CASE WHEN c_icd IS NULL THEN 'NULL::text' ELSE FORMAT('%I::text',c_icd) END,
        CASE WHEN c_icdt IS NULL THEN 'NULL::text' ELSE FORMAT('%I::text',c_icdt) END,
        v_schema,t_cases
    );
END $$;

/* Standardize SICdb data_ref for structured preconditions. */
CREATE TABLE study_make30_final.sic_src_data_ref(case_id bigint,field_id integer,ref_id integer);
DO $$
DECLARE v_schema text;t text;c_case text;c_field text;c_ref text;
BEGIN
    SELECT n.nspname INTO v_schema FROM pg_namespace n
    WHERE EXISTS(SELECT 1 FROM pg_class c WHERE c.relnamespace=n.oid AND LOWER(c.relname)='data_ref' AND c.relkind IN('r','p','v','m'))
      AND EXISTS(SELECT 1 FROM pg_class c WHERE c.relnamespace=n.oid AND LOWER(c.relname)='cases' AND c.relkind IN('r','p','v','m'))
    ORDER BY CASE WHEN n.nspname='public' THEN 0 ELSE 1 END,n.nspname LIMIT 1;
    IF v_schema IS NOT NULL THEN
      SELECT table_name INTO t FROM information_schema.tables WHERE table_schema=v_schema AND LOWER(table_name)='data_ref' LIMIT 1;
      SELECT column_name INTO c_case FROM information_schema.columns WHERE table_schema=v_schema AND table_name=t AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='caseid' LIMIT 1;
      SELECT column_name INTO c_field FROM information_schema.columns WHERE table_schema=v_schema AND table_name=t AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='fieldid' LIMIT 1;
      SELECT column_name INTO c_ref FROM information_schema.columns WHERE table_schema=v_schema AND table_name=t AND REGEXP_REPLACE(LOWER(column_name),'[^a-z0-9]','','g')='refid' LIMIT 1;
      IF c_case IS NOT NULL AND c_field IS NOT NULL AND c_ref IS NOT NULL THEN
        EXECUTE FORMAT('INSERT INTO study_make30_final.sic_src_data_ref SELECT %I::bigint,%I::integer,%I::integer FROM %I.%I WHERE %I IS NOT NULL',c_case,c_field,c_ref,v_schema,t,c_case);
      END IF;
    END IF;
END $$;
CREATE INDEX idx_sic_data_ref_case ON study_make30_final.sic_src_data_ref(case_id);

/* Creatinine dictionary: official SICdb examples use ID 367; semantic matches
   are also accepted so locally equivalent dictionary builds remain auditable. */
CREATE TABLE study_make30_final.sic_creatinine_ids AS
SELECT DISTINCT r.ref_id AS lab_id,r.ref_value,r.ref_unit
FROM study_make30_final.sic_src_refs r
WHERE r.ref_id=367 OR LOWER(COALESCE(r.ref_value,'')) ~ '(creatinine|kreatinin)';
CREATE UNIQUE INDEX idx_sic_creatinine_ids ON study_make30_final.sic_creatinine_ids(lab_id);

CREATE TABLE study_make30_final.sic_baseline_kidney AS
WITH cr AS
(
    SELECT
        l.case_id,sc.patient_id,l.offset_sec,
        CASE
            WHEN LOWER(COALESCE(ci.ref_unit,'')) ~ '(umol|µmol|micromol)' THEN l.lab_value/88.4
            ELSE l.lab_value
        END::double precision AS creatinine_mg_dl,
        sc.offset_after_first_admission_sec
    FROM study_make30_final.sic_src_laboratory l
    INNER JOIN study_make30_final.sic_src_cases sc ON sc.case_id=l.case_id
    INNER JOIN study_make30_final.sic_creatinine_ids ci ON ci.lab_id=l.lab_id
),
hist AS
(
    SELECT
        c.case_id,
        PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY x.creatinine_mg_dl)::double precision AS historical_creatinine_median
    FROM study_make30_final.sic_core4_trajectory_cohort c
    INNER JOIN cr x ON x.patient_id=c.patient_id
    WHERE x.creatinine_mg_dl BETWEEN 0.1 AND 20
      AND x.offset_after_first_admission_sec+x.offset_sec >= c.offset_after_first_admission_sec-31536000
      AND x.offset_after_first_admission_sec+x.offset_sec <  c.offset_after_first_admission_sec
    GROUP BY c.case_id
),
first_admit AS
(
    SELECT DISTINCT ON(c.case_id)
        c.case_id,x.offset_sec AS first_admission_creatinine_offset_sec,x.creatinine_mg_dl AS first_admission_creatinine
    FROM study_make30_final.sic_core4_trajectory_cohort c
    INNER JOIN cr x ON x.case_id=c.case_id
    WHERE x.creatinine_mg_dl BETWEEN 0.1 AND 20
      AND x.offset_sec>=0
      AND x.offset_sec<c.landmark_offset_sec
    ORDER BY c.case_id,x.offset_sec
),
sex AS
(
    SELECT e.case_id,LOWER(COALESCE(r.ref_value,'')) AS sex_text
    FROM study_make30_final.sic_case_extra e LEFT JOIN study_make30_final.sic_src_refs r ON r.ref_id=e.sex_ref_id
),
calc AS
(
    SELECT
        c.patient_id,c.case_id,c.age_years,
        s.sex_text,h.historical_creatinine_median,f.first_admission_creatinine_offset_sec,f.first_admission_creatinine,
        POWER(175.0*POWER(GREATEST(c.age_years::double precision,18.0),-0.203)
              *CASE WHEN s.sex_text ~ '(female|weiblich|^f$)' THEN 0.742 ELSE 1.0 END/75.0,1.0/1.154)::double precision AS mdrd_backcalc_creatinine
    FROM study_make30_final.sic_core4_trajectory_cohort c
    LEFT JOIN hist h USING(case_id)
    LEFT JOIN first_admit f USING(case_id)
    LEFT JOIN sex s USING(case_id)
),base AS
(
    SELECT x.*,
        LEAST(COALESCE(x.historical_creatinine_median,x.mdrd_backcalc_creatinine),
              COALESCE(x.first_admission_creatinine,x.historical_creatinine_median,x.mdrd_backcalc_creatinine))::double precision AS baseline_creatinine
    FROM calc x
)
SELECT b.*,
    CASE WHEN b.historical_creatinine_median IS NOT NULL AND b.baseline_creatinine=b.historical_creatinine_median THEN '12mo historical median'
         WHEN b.first_admission_creatinine IS NOT NULL AND b.baseline_creatinine=b.first_admission_creatinine THEN 'first admission creatinine'
         ELSE 'MDRD back-calculation eGFR=75' END::text AS baseline_creatinine_source,
    CASE WHEN b.baseline_creatinine>0 THEN
        175.0*POWER(b.baseline_creatinine,-1.154)*POWER(GREATEST(b.age_years::double precision,18.0),-0.203)
        *CASE WHEN b.sex_text ~ '(female|weiblich|^f$)' THEN 0.742 ELSE 1.0 END
    END::double precision AS baseline_egfr_mdrd
FROM base b;
CREATE UNIQUE INDEX idx_sic_baseline_kidney_case ON study_make30_final.sic_baseline_kidney(case_id);

CREATE TABLE study_make30_final.sic_covariates AS
WITH pre AS
(
    SELECT
        d.case_id,
        MAX(CASE WHEN REGEXP_REPLACE(LOWER(COALESCE(fr.ref_value,'')),'[^a-z0-9]','','g')='preconditiondiabetes'
                  THEN CASE WHEN LOWER(COALESCE(vr.ref_value,'')) ~ '^(yes|ja|true|1|present|positive)$' THEN 1 WHEN LOWER(COALESCE(vr.ref_value,'')) ~ '^(no|nein|false|0|absent|negative)$' THEN 0 END END)::integer AS diabetes,
        MAX(CASE WHEN REGEXP_REPLACE(LOWER(COALESCE(fr.ref_value,'')),'[^a-z0-9]','','g')='preconditionrenaldysfunction'
                  THEN CASE WHEN LOWER(COALESCE(vr.ref_value,'')) ~ '^(yes|ja|true|1|present|positive)$' THEN 1 WHEN LOWER(COALESCE(vr.ref_value,'')) ~ '^(no|nein|false|0|absent|negative)$' THEN 0 END END)::integer AS ckd_or_renal_disease
    FROM study_make30_final.sic_src_data_ref d
    LEFT JOIN study_make30_final.sic_src_refs fr ON fr.ref_id=d.field_id
    LEFT JOIN study_make30_final.sic_src_refs vr ON vr.ref_id=d.ref_id
    GROUP BY d.case_id
),
vaso_ids AS
(
    SELECT DISTINCT ref_id
    FROM study_make30_final.sic_src_refs
    WHERE LOWER(COALESCE(ref_value,'')||' '||COALESCE(ref_name,'')) ~ '(norepineph|noradren|epineph|adrenalin|vasopressin|dopamine|phenylephrine|dobutamine)'
),
vaso_status AS
(
    SELECT COUNT(*)::integer AS recognized_vaso_ids FROM vaso_ids
),
vaso AS
(
    SELECT c.case_id,1::integer AS vasopressor_0_24h
    FROM study_make30_final.sic_core4_trajectory_cohort c
    WHERE EXISTS
    (
       SELECT 1 FROM study_make30_final.sic_src_medication m INNER JOIN vaso_ids v ON v.ref_id=m.drug_id
       WHERE m.case_id=c.case_id
         AND COALESCE(m.offset_end_sec,m.offset_sec+60)>c.sepsis_offset_sec
         AND m.offset_sec<c.landmark_offset_sec
    )
),
mv_ids AS
(
    SELECT DISTINCT ref_id
    FROM study_make30_final.sic_src_refs
    WHERE LOWER(COALESCE(ref_value,'')||' '||COALESCE(ref_name,'')) ~ '(mechanical ventil|invasive ventil|intubat|endotracheal)'
),
mv_status AS
(
    SELECT COUNT(*)::integer AS recognized_mv_ids FROM mv_ids
),
mv AS
(
    SELECT c.case_id,1::integer AS mechanical_ventilation_0_24h
    FROM study_make30_final.sic_core4_trajectory_cohort c
    WHERE EXISTS
    (
       SELECT 1 FROM study_make30_final.sic_src_data_range r INNER JOIN mv_ids v ON v.ref_id=r.data_id
       WHERE r.case_id=c.case_id
         AND COALESCE(r.offset_end_sec,r.offset_sec)>c.sepsis_offset_sec
         AND r.offset_sec<c.landmark_offset_sec
    )
),
sex AS
(
    SELECT e.case_id,r.ref_value AS sex
    FROM study_make30_final.sic_case_extra e LEFT JOIN study_make30_final.sic_src_refs r ON r.ref_id=e.sex_ref_id
)
SELECT
    c.patient_id,c.case_id,c.age_years::double precision AS age,s.sex,
    NULL::text AS race_ethnicity,
    CASE WHEN e.weight_kg BETWEEN 25 AND 300 THEN e.weight_kg END::double precision AS weight_kg,
    CASE WHEN e.height_cm BETWEEN 100 AND 250 THEN e.height_cm END::double precision AS height_cm,
    CASE WHEN e.weight_kg BETWEEN 25 AND 300 AND e.height_cm BETWEEN 100 AND 250 THEN e.weight_kg/POWER(e.height_cm/100.0,2) END::double precision AS bmi,
    c.sofa_score::double precision AS sepsis_sofa_score,
    e.saps3,
    bk.historical_creatinine_median,bk.first_admission_creatinine,bk.mdrd_backcalc_creatinine,bk.baseline_creatinine_source,
    bk.baseline_creatinine,bk.baseline_egfr_mdrd,
    p.ckd_or_renal_disease,p.diabetes,
    NULL::integer AS congestive_heart_failure,
    NULL::integer AS liver_disease,
    NULL::integer AS malignancy,
    CASE WHEN ms.recognized_mv_ids=0 THEN NULL ELSE COALESCE(m.mechanical_ventilation_0_24h,0)::integer END AS mechanical_ventilation_0_24h,
    CASE WHEN vs.recognized_vaso_ids=0 THEN NULL ELSE COALESCE(v.vasopressor_0_24h,0)::integer END AS vasopressor_0_24h,
    e.icd10_main,e.icd10_main_text
FROM study_make30_final.sic_core4_trajectory_cohort c
LEFT JOIN study_make30_final.sic_case_extra e USING(case_id)
LEFT JOIN sex s USING(case_id)
LEFT JOIN pre p USING(case_id)
LEFT JOIN study_make30_final.sic_baseline_kidney bk USING(case_id)
LEFT JOIN mv m USING(case_id)
CROSS JOIN mv_status ms
CROSS JOIN vaso_status vs
LEFT JOIN vaso v USING(case_id);
CREATE UNIQUE INDEX idx_sic_covariates_case ON study_make30_final.sic_covariates(case_id);

CREATE TABLE study_make30_final.sic_make30 AS
WITH fup AS
(
    SELECT
        c.patient_id,c.case_id,c.landmark_offset_sec,c.icu_offset_sec,c.offset_of_death_sec,c.hospital_discharge_day,c.hospital_stay_days,c.time_of_stay_sec,
        LEAST(
            COALESCE(c.hospital_stay_days::bigint*86400,c.hospital_discharge_day::bigint*86400,c.time_of_stay_sec),
            c.icu_offset_sec+2592000
        )::bigint AS make30_censor_offset_sec
    FROM study_make30_final.sic_core4_trajectory_cohort c
),
all_krt_events AS
(
    SELECT c.case_id,df.offset_sec AS event_time
    FROM study_make30_final.sic_core4_trajectory_cohort c
    INNER JOIN study_make30_final.sic_src_data_float_h df ON df.case_id=c.case_id
    INNER JOIN study_make30_final.sic_krt_ids k ON k.data_id=df.data_id
    WHERE df.val IS NOT NULL AND ABS(df.val)>0
    UNION ALL
    SELECT c.case_id,r.offset_sec AS event_time
    FROM study_make30_final.sic_core4_trajectory_cohort c
    INNER JOIN study_make30_final.sic_src_data_range r ON r.case_id=c.case_id
    INNER JOIN study_make30_final.sic_krt_ids k ON k.data_id=r.data_id
),
prior_krt AS
(
    SELECT f.case_id,1::integer AS prior_krt_before_landmark
    FROM fup f
    WHERE EXISTS
    (
        SELECT 1 FROM all_krt_events e
        WHERE e.case_id=f.case_id AND e.event_time<=f.landmark_offset_sec
    )
),
krt AS
(
    SELECT f.case_id,MIN(e.event_time)::bigint AS first_new_krt_offset_sec
    FROM fup f
    LEFT JOIN prior_krt pk USING(case_id)
    INNER JOIN all_krt_events e ON e.case_id=f.case_id
    WHERE pk.case_id IS NULL
      AND e.event_time>f.landmark_offset_sec
      AND e.event_time<=f.make30_censor_offset_sec
    GROUP BY f.case_id
),
cr AS
(
    SELECT l.case_id,l.offset_sec,
           CASE WHEN LOWER(COALESCE(ci.ref_unit,'')) ~ '(umol|µmol|micromol)' THEN l.lab_value/88.4 ELSE l.lab_value END::double precision AS creatinine_mg_dl
    FROM study_make30_final.sic_src_laboratory l INNER JOIN study_make30_final.sic_creatinine_ids ci ON ci.lab_id=l.lab_id
),
final_cr AS
(
    SELECT DISTINCT ON(f.case_id)
        f.case_id,x.offset_sec AS final_inpatient_creatinine_offset_sec,x.creatinine_mg_dl AS final_inpatient_creatinine
    FROM fup f INNER JOIN cr x ON x.case_id=f.case_id
    WHERE x.offset_sec>f.landmark_offset_sec AND x.offset_sec<=f.make30_censor_offset_sec AND x.creatinine_mg_dl BETWEEN 0.1 AND 20
    ORDER BY f.case_id,x.offset_sec DESC
),
discharge AS
(
    SELECT
        f.case_id,
        COALESCE(f.hospital_stay_days::bigint*86400,
                 f.hospital_discharge_day::bigint*86400,
                 f.time_of_stay_sec)::bigint AS hospital_discharge_offset_sec
    FROM fup f
)
SELECT
    f.patient_id,f.case_id,
    f.landmark_offset_sec AS make30_risk_start_offset_sec,f.make30_censor_offset_sec,d.hospital_discharge_offset_sec,
    (f.make30_censor_offset_sec-f.landmark_offset_sec)/86400.0 AS make30_followup_days_after_landmark,
    CASE WHEN f.make30_censor_offset_sec>f.landmark_offset_sec THEN 1 ELSE 0 END::integer AS make30_evaluable_window,
    CASE WHEN f.offset_of_death_sec IS NOT NULL
                   AND f.offset_of_death_sec>f.landmark_offset_sec
                   AND f.offset_of_death_sec<=f.make30_censor_offset_sec
         THEN 1 ELSE 0 END::integer AS make30_death,
    COALESCE(pk.prior_krt_before_landmark,0)::integer AS prior_krt_before_landmark,
    CASE WHEN k.first_new_krt_offset_sec IS NOT NULL THEN 1 ELSE 0 END::integer AS make30_new_krt,
    k.first_new_krt_offset_sec,
    fc.final_inpatient_creatinine,fc.final_inpatient_creatinine_offset_sec,bk.baseline_creatinine,
    CASE WHEN fc.final_inpatient_creatinine IS NOT NULL AND bk.baseline_creatinine IS NOT NULL AND fc.final_inpatient_creatinine>=2.0*bk.baseline_creatinine THEN 1
         WHEN fc.final_inpatient_creatinine IS NOT NULL AND bk.baseline_creatinine IS NOT NULL THEN 0 ELSE NULL END::integer AS make30_persistent_kidney_dysfunction,
    CASE
      WHEN f.make30_censor_offset_sec<=f.landmark_offset_sec THEN NULL
      WHEN f.offset_of_death_sec IS NOT NULL
           AND f.offset_of_death_sec>f.landmark_offset_sec
           AND f.offset_of_death_sec<=f.make30_censor_offset_sec THEN 1
      WHEN k.first_new_krt_offset_sec IS NOT NULL THEN 1
      WHEN fc.final_inpatient_creatinine IS NOT NULL AND bk.baseline_creatinine IS NOT NULL AND fc.final_inpatient_creatinine>=2.0*bk.baseline_creatinine THEN 1
      WHEN fc.final_inpatient_creatinine IS NOT NULL AND bk.baseline_creatinine IS NOT NULL THEN 0
      ELSE NULL
    END::integer AS make30
FROM fup f
LEFT JOIN prior_krt pk USING(case_id)
LEFT JOIN krt k USING(case_id)
LEFT JOIN final_cr fc USING(case_id)
LEFT JOIN study_make30_final.sic_baseline_kidney bk USING(case_id)
LEFT JOIN discharge d USING(case_id);
CREATE UNIQUE INDEX idx_sic_make30_case ON study_make30_final.sic_make30(case_id);

CREATE TABLE study_make30_final.sic_analysis_master AS
SELECT
    w.*,
    cv.age,cv.sex,cv.race_ethnicity,cv.weight_kg,cv.height_cm,cv.bmi,cv.sepsis_sofa_score,cv.saps3,
    cv.historical_creatinine_median,cv.first_admission_creatinine,cv.mdrd_backcalc_creatinine,cv.baseline_creatinine_source,
    cv.baseline_creatinine,cv.baseline_egfr_mdrd,cv.ckd_or_renal_disease,cv.diabetes,
    cv.congestive_heart_failure,cv.liver_disease,cv.malignancy,
    cv.mechanical_ventilation_0_24h,cv.vasopressor_0_24h,cv.icd10_main,cv.icd10_main_text,
    o.make30_risk_start_offset_sec,o.make30_censor_offset_sec,o.hospital_discharge_offset_sec,o.make30_followup_days_after_landmark,o.make30_evaluable_window,
    o.make30_death,o.prior_krt_before_landmark,o.make30_new_krt,o.first_new_krt_offset_sec,
    o.final_inpatient_creatinine,o.final_inpatient_creatinine_offset_sec,o.make30_persistent_kidney_dysfunction,o.make30
FROM study_make30_final.sic_core4_multitrajectory_wide w
LEFT JOIN study_make30_final.sic_covariates cv USING(case_id)
LEFT JOIN study_make30_final.sic_make30 o USING(case_id);
CREATE UNIQUE INDEX idx_sic_analysis_master_case ON study_make30_final.sic_analysis_master(case_id);


/* =========================================================================
   26. FINAL QC + ANALYSIS RESULT
   ========================================================================= */
CREATE TABLE study_make30_final.sic_qc_summary AS
WITH implementation_info AS
(
 SELECT '00_IMPLEMENTATION_INFO'::text section,0 sort1,1 sort2,'SICdb mixed-frequency trajectory: MSI 2-h x12 from data_float_h; Na/K/Cl/HCO3 6-h x4 with central-lab priority and BGA fallback; primary MSI completeness >=9/12'::text item,NULL::text time_bin,NULL::bigint n_total,NULL::bigint n_observed,NULL::bigint n_missing,NULL::numeric observed_pct,NULL::numeric missing_pct,NULL::bigint n_meet_requirement,NULL::bigint n_all_bins,NULL::numeric mean_value,NULL::numeric sd_value,NULL::numeric min_value,NULL::numeric p01,NULL::numeric p05,NULL::numeric p25,NULL::numeric median_value,NULL::numeric p75,NULL::numeric p95,NULL::numeric p99,NULL::numeric max_value,NULL::text detail
),
cohort_flow AS
(
 SELECT '01_COHORT_FLOW'::text section,1 sort1,1 sort2,'all unique SICdb patients'::text item,NULL::text time_bin,COUNT(*)::bigint n_total,COUNT(*)::bigint n_observed,NULL::bigint n_missing,100::numeric observed_pct,0::numeric missing_pct,NULL::bigint n_meet_requirement,NULL::bigint n_all_bins,NULL::numeric mean_value,NULL::numeric sd_value,NULL::numeric min_value,NULL::numeric p01,NULL::numeric p05,NULL::numeric p25,NULL::numeric median_value,NULL::numeric p75,NULL::numeric p95,NULL::numeric p99,NULL::numeric max_value,NULL::text detail FROM study_make30_final.sic_patients
 UNION ALL SELECT '01_COHORT_FLOW',1,2,'adult patients with ORAKLE-style Sepsis onset after ICUOffset',NULL,(SELECT COUNT(*) FROM study_make30_final.sic_patients),COUNT(DISTINCT patient_id),NULL,ROUND(100.0*COUNT(DISTINCT patient_id)/NULLIF((SELECT COUNT(*) FROM study_make30_final.sic_patients),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.sic_sepsis_adult
 UNION ALL SELECT '01_COHORT_FLOW',1,3,'first qualifying Sepsis case per patient',NULL,(SELECT COUNT(*) FROM study_make30_final.sic_patients),COUNT(*),NULL,ROUND(100.0*COUNT(*)/NULLIF((SELECT COUNT(*) FROM study_make30_final.sic_patients),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.sic_sepsis_index
 UNION ALL SELECT '01_COHORT_FLOW',1,4,'reached Sepsis +24 h landmark',NULL,COUNT(*),COUNT(*) FILTER(WHERE reached_24h_landmark=1),COUNT(*) FILTER(WHERE reached_24h_landmark=0),ROUND(100.0*COUNT(*) FILTER(WHERE reached_24h_landmark=1)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE reached_24h_landmark=0)/NULLIF(COUNT(*),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.sic_landmark_cohort
 UNION ALL SELECT '01_COHORT_FLOW',1,5,'trajectory pool after landmark + early-KRT exclusion',NULL,(SELECT COUNT(*) FROM study_make30_final.sic_landmark_cohort WHERE reached_24h_landmark=1),COUNT(*),NULL,ROUND(100.0*COUNT(*)/NULLIF((SELECT COUNT(*) FROM study_make30_final.sic_landmark_cohort WHERE reached_24h_landmark=1),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.sic_trajectory_pool
 UNION ALL SELECT '01_COHORT_FLOW',1,6,'PRIMARY mixed-frequency Core-4 eligible: MSI>=9/12 + Na/K/Cl/HCO3 each>=2/4',NULL,COUNT(*),COUNT(*) FILTER(WHERE core4_trajectory_eligible=1),COUNT(*) FILTER(WHERE core4_trajectory_eligible=0),ROUND(100.0*COUNT(*) FILTER(WHERE core4_trajectory_eligible=1)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE core4_trajectory_eligible=0)/NULLIF(COUNT(*),0),2),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.sic_trajectory_coverage
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
 FROM study_make30_final.sic_trajectory_coverage c CROSS JOIN variable_list v GROUP BY v.variable_order,v.variable_name,v.required_bins,v.max_bins
),
all_series AS
(
 SELECT case_id,1 variable_order,'MSI'::text variable_name,time_bin,time_label,msi::double precision value_num FROM study_make30_final.sic_msi_trajectory_long
 UNION ALL SELECT case_id,2,'sodium',time_bin,time_label,sodium::double precision FROM study_make30_final.sic_electrolyte_trajectory_long
 UNION ALL SELECT case_id,3,'potassium',time_bin,time_label,potassium::double precision FROM study_make30_final.sic_electrolyte_trajectory_long
 UNION ALL SELECT case_id,4,'chloride',time_bin,time_label,chloride::double precision FROM study_make30_final.sic_electrolyte_trajectory_long
 UNION ALL SELECT case_id,5,'bicarbonate',time_bin,time_label,bicarbonate::double precision FROM study_make30_final.sic_electrolyte_trajectory_long
),
eligible_series AS
(
 SELECT case_id,1 variable_order,'MSI'::text variable_name,time_bin,time_label,msi::double precision value_num FROM study_make30_final.sic_core4_msi_trajectory_long
 UNION ALL SELECT case_id,2,'sodium',time_bin,time_label,sodium::double precision FROM study_make30_final.sic_core4_electrolyte_trajectory_long
 UNION ALL SELECT case_id,3,'potassium',time_bin,time_label,potassium::double precision FROM study_make30_final.sic_core4_electrolyte_trajectory_long
 UNION ALL SELECT case_id,4,'chloride',time_bin,time_label,chloride::double precision FROM study_make30_final.sic_core4_electrolyte_trajectory_long
 UNION ALL SELECT case_id,5,'bicarbonate',time_bin,time_label,bicarbonate::double precision FROM study_make30_final.sic_core4_electrolyte_trajectory_long
),
all_missing AS (SELECT '03_BIN_LEVEL_MISSINGNESS'::text section,3 sort1,(variable_order*100+time_bin)::integer sort2,variable_name item,time_label time_bin,COUNT(*)::bigint n_total,COUNT(value_num)::bigint n_observed,COUNT(*) FILTER(WHERE value_num IS NULL)::bigint n_missing,ROUND(100.0*COUNT(value_num)/NULLIF(COUNT(*),0),2) observed_pct,ROUND(100.0*COUNT(*) FILTER(WHERE value_num IS NULL)/NULLIF(COUNT(*),0),2) missing_pct,NULL::bigint n_meet_requirement,NULL::bigint n_all_bins,NULL::numeric mean_value,NULL::numeric sd_value,NULL::numeric min_value,NULL::numeric p01,NULL::numeric p05,NULL::numeric p25,NULL::numeric median_value,NULL::numeric p75,NULL::numeric p95,NULL::numeric p99,NULL::numeric max_value,NULL::text detail FROM all_series GROUP BY variable_order,variable_name,time_bin,time_label),
elig_missing AS (SELECT '04_CORE4_ELIGIBLE_BIN_MISSINGNESS'::text section,4 sort1,(variable_order*100+time_bin)::integer sort2,variable_name item,time_label time_bin,COUNT(*)::bigint,COUNT(value_num)::bigint,COUNT(*) FILTER(WHERE value_num IS NULL)::bigint,ROUND(100.0*COUNT(value_num)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE value_num IS NULL)/NULLIF(COUNT(*),0),2),NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::text FROM eligible_series GROUP BY variable_order,variable_name,time_bin,time_label),
value_qc AS (SELECT '05_CORE4_VALUE_DISTRIBUTION_QC'::text section,5 sort1,(variable_order*100+time_bin)::integer sort2,variable_name item,time_label time_bin,COUNT(*)::bigint,COUNT(value_num)::bigint,COUNT(*) FILTER(WHERE value_num IS NULL)::bigint,ROUND(100.0*COUNT(value_num)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE value_num IS NULL)/NULLIF(COUNT(*),0),2),NULL::bigint,NULL::bigint,AVG(value_num)::numeric,STDDEV_SAMP(value_num)::numeric,MIN(value_num)::numeric,PERCENTILE_CONT(.01) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.05) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.25) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.5) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.75) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.95) WITHIN GROUP(ORDER BY value_num)::numeric,PERCENTILE_CONT(.99) WITHIN GROUP(ORDER BY value_num)::numeric,MAX(value_num)::numeric,NULL::text FROM eligible_series GROUP BY variable_order,variable_name,time_bin,time_label),
threshold_qc AS
(
 SELECT '06_MSI_COVERAGE_THRESHOLD_SENSITIVITY'::text section,6 sort1,x.ord sort2,x.label item,NULL::text time_bin,COUNT(*)::bigint,COUNT(*) FILTER(WHERE x.flag=1)::bigint,COUNT(*) FILTER(WHERE x.flag=0)::bigint,ROUND(100.0*COUNT(*) FILTER(WHERE x.flag=1)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE x.flag=0)/NULLIF(COUNT(*),0),2),NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::text
 FROM study_make30_final.sic_trajectory_coverage c CROSS JOIN LATERAL (VALUES (1,'Core-4 if MSI>=6/12',c.core4_eligible_msi6),(2,'Core-4 if MSI>=9/12 [PRIMARY]',c.core4_eligible_msi9),(3,'Core-4 if MSI=12/12',c.core4_eligible_msi12)) x(ord,label,flag) GROUP BY x.ord,x.label
),
anomaly_trace AS
(
 SELECT '07_MSI_GT4_TRACE'::text section,7 sort1,ROW_NUMBER() OVER(ORDER BY hourly_msi DESC,case_id,hour_bin)::integer sort2,'hourly MSI > 4'::text item,msi_2h_label::text time_bin,1::bigint,1::bigint,0::bigint,100::numeric,0::numeric,NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,
 jsonb_build_object('patient_id',patient_id,'case_id',case_id,'hour_bin',hour_bin,'msi_2h_bin',msi_2h_bin,'hr_ecg',hr_ecg,'hr_spo2',hr_spo2,'hr_abp',hr_abp,'hr_used',hr_used,'hr_source',hr_source,'map_invasive',map_invasive,'map_noninvasive',map_noninvasive,'map_abs_diff',map_abs_diff,'map_used',map_used,'map_source',map_source,'hourly_msi',hourly_msi,'msi_2h_max',msi_2h_max)::text FROM study_make30_final.sic_msi_gt4_trace
),
structure_qc AS
(
 SELECT '08_STRUCTURE_QC'::text section,8 sort1,1 sort2,'eligible patients'::text item,NULL::text time_bin,COUNT(*)::bigint,NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::bigint,NULL::bigint,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::numeric,NULL::text FROM study_make30_final.sic_core4_trajectory_cohort
 UNION ALL SELECT '08_STRUCTURE_QC',8,2,'MSI long rows (expected eligible N x 12)',NULL,COUNT(*),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.sic_core4_msi_trajectory_long
 UNION ALL SELECT '08_STRUCTURE_QC',8,3,'electrolyte long rows (expected eligible N x 4)',NULL,COUNT(*),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL FROM study_make30_final.sic_core4_electrolyte_trajectory_long
)
SELECT section,item,time_bin,n_total,n_observed,n_missing,observed_pct,missing_pct,n_meet_requirement,n_all_bins,ROUND(mean_value,4) mean_value,ROUND(sd_value,4) sd_value,ROUND(min_value,4) min_value,ROUND(p01,4) p01,ROUND(p05,4) p05,ROUND(p25,4) p25,ROUND(median_value,4) median_value,ROUND(p75,4) p75,ROUND(p95,4) p95,ROUND(p99,4) p99,ROUND(max_value,4) max_value,detail
FROM (SELECT * FROM implementation_info UNION ALL SELECT * FROM cohort_flow UNION ALL SELECT * FROM patient_coverage UNION ALL SELECT * FROM all_missing UNION ALL SELECT * FROM elig_missing UNION ALL SELECT * FROM value_qc UNION ALL SELECT * FROM threshold_qc UNION ALL SELECT * FROM anomaly_trace UNION ALL SELECT * FROM structure_qc) z ORDER BY sort1,sort2;


INSERT INTO study_make30_final.sic_qc_summary
(section,item,time_bin,n_total,n_observed,n_missing,observed_pct,missing_pct,n_meet_requirement,n_all_bins,mean_value,sd_value,min_value,p01,p05,p25,median_value,p75,p95,p99,max_value,detail)
SELECT '09_MAKE30_QC','MAKE30 evaluable / events',NULL,COUNT(*)::bigint,COUNT(make30)::bigint,COUNT(*) FILTER(WHERE make30 IS NULL)::bigint,
       ROUND(100.0*COUNT(make30)/NULLIF(COUNT(*),0),2),ROUND(100.0*COUNT(*) FILTER(WHERE make30 IS NULL)/NULLIF(COUNT(*),0),2),COUNT(*) FILTER(WHERE make30=1)::bigint,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
       jsonb_build_object('death',COUNT(*) FILTER(WHERE make30_death=1),'new_krt',COUNT(*) FILTER(WHERE make30_new_krt=1),'persistent_kidney_dysfunction',COUNT(*) FILTER(WHERE make30_persistent_kidney_dysfunction=1),'missing_final_creatinine',COUNT(*) FILTER(WHERE final_inpatient_creatinine IS NULL))::text
FROM study_make30_final.sic_make30;

INSERT INTO study_make30_final.sic_qc_summary
(section,item,time_bin,n_total,n_observed,n_missing,observed_pct,missing_pct,n_meet_requirement,n_all_bins,mean_value,sd_value,min_value,p01,p05,p25,median_value,p75,p95,p99,max_value,detail)
SELECT '10_COVARIATE_QC',v.var,NULL,COUNT(*)::bigint,COUNT(v.val)::bigint,(COUNT(*)-COUNT(v.val))::bigint,
       ROUND(100.0*COUNT(v.val)/NULLIF(COUNT(*),0),2),ROUND(100.0*(COUNT(*)-COUNT(v.val))/NULLIF(COUNT(*),0),2),NULL,NULL,
       AVG(v.val),STDDEV_SAMP(v.val),MIN(v.val),PERCENTILE_CONT(.01) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.05) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.25) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.5) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.75) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.95) WITHIN GROUP(ORDER BY v.val),PERCENTILE_CONT(.99) WITHIN GROUP(ORDER BY v.val),MAX(v.val),NULL
FROM study_make30_final.sic_covariates c CROSS JOIN LATERAL (VALUES ('age',c.age),('bmi',c.bmi),('sepsis_sofa_score',c.sepsis_sofa_score),('baseline_creatinine',c.baseline_creatinine),('baseline_egfr_mdrd',c.baseline_egfr_mdrd),('saps3',c.saps3)) v(var,val)
GROUP BY v.var;

SELECT * FROM study_make30_final.sic_analysis_master ORDER BY patient_id,case_id;
