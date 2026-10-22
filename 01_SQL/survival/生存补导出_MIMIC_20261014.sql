-- Read-only supplement. Run in the MIMIC database containing the frozen pipeline tables.
-- Export one row per stay_id; merge onto frozen classifications without changing them.
-- Not executed against the database in this review.
SELECT
    m.subject_id, m.hadm_id, m.stay_id,
    c.icu_intime, c.icu_outtime, c.dischtime,
    c.sepsis_time, c.landmark_time,
    a.hospital_expire_flag,
    a.deathtime AS recorded_death_time,
    a.dischtime AS admission_discharge_time,
    m.death_time AS frozen_death_time,
    CASE WHEN a.hospital_expire_flag = 1 AND a.deathtime IS NULL
         THEN 1 ELSE 0 END AS frozen_death_used_discharge_fallback,
    c.icu_intime + INTERVAL '30 days' AS administrative_end_time,
    m.make30_risk_start, m.make30_censor_time,
    m.make30_evaluable_window, m.make30_death
FROM study_make30_final.mimic_analysis_master m
JOIN study_make30_final.mimic_core4_trajectory_cohort c USING (stay_id)
JOIN mimiciv_hosp.admissions a ON a.hadm_id = m.hadm_id
ORDER BY m.stay_id;
