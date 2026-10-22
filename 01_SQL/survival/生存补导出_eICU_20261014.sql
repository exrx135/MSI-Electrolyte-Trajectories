-- Read-only supplement. Native offsets below are MINUTES from index ICU admission.
-- Keep discharge status verbatim, including unknown values, for endpoint audit.
-- Historical SQL output has 5867 rows; join to frozen 5839-row cohort by patientunitstayid.
-- Do not reintroduce the 28 excluded records. Not database-executed in this review.
SELECT
    m.uniquepid, m.patienthealthsystemstayid, m.patientunitstayid,
    c.hospitaladmitoffset, c.hospitaldischargeoffset,
    c.hospitaldischargestatus,
    c.unitdischargeoffset, c.unitdischargestatus,
    c.sepsis_offset, c.landmark_offset,
    CASE WHEN LOWER(COALESCE(c.hospitaldischargestatus, '')) = 'expired'
         THEN c.hospitaldischargeoffset END AS death_icu_offset_min_from_status,
    43200::integer AS administrative_end_icu_offset_min,
    m.make30_risk_start_hosp_offset, m.make30_censor_hosp_offset,
    m.make30_evaluable_window, m.make30_death
FROM study_make30_final.eicu_analysis_master m
JOIN study_make30_final.eicu_core4_trajectory_cohort c USING (patientunitstayid)
ORDER BY m.patientunitstayid;
