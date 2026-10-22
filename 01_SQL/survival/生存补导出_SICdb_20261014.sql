-- Read-only supplement using existing normalized source fields.
-- Preserve original day-level discharge components; do not invent finer precision.
-- Offset origin and fallback semantics must be checked before final survival analysis.
-- Not executed against the database in this review.
SELECT
    m.patient_id, m.case_id,
    c.icu_offset_sec, c.time_of_stay_sec,
    c.offset_after_first_admission_sec,
    c.offset_of_death_sec,
    c.hospital_discharge_day, c.hospital_stay_days,
    c.sepsis_offset_sec, c.landmark_offset_sec,
    c.icu_offset_sec + 2592000 AS administrative_end_offset_sec,
    CASE WHEN c.hospital_stay_days IS NOT NULL THEN 'hospital_stay_days'
         WHEN c.hospital_discharge_day IS NOT NULL THEN 'hospital_discharge_day'
         WHEN c.time_of_stay_sec IS NOT NULL THEN 'time_of_stay_sec'
         ELSE 'missing' END AS frozen_discharge_endpoint_source,
    m.hospital_discharge_offset_sec,
    m.make30_risk_start_offset_sec, m.make30_censor_offset_sec,
    m.make30_evaluable_window, m.make30_death
FROM study_make30_final.sic_analysis_master m
JOIN study_make30_final.sic_core4_trajectory_cohort c USING (case_id)
ORDER BY m.case_id;
