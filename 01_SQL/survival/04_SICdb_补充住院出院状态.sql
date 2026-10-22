-- Read-only: normalized fields already created by the original extraction SQL.
-- Export the entire result as 04_SICdb_补充住院出院状态.csv in this folder.
-- HospitalStayDays includes pre-ICU hospital time: do NOT use it directly as
-- a MetaVision-admission-origin endpoint. Preserve both day fields for audit.
SELECT m.patient_id, m.case_id,
       e.hospital_discharge_type_ref,
       r.ref_value AS hospital_discharge_type_value,
       r.ref_name AS hospital_discharge_type_name,
       c.icu_offset_sec, c.time_of_stay_sec,
       c.offset_after_first_admission_sec,
       c.offset_of_death_sec,
       c.hospital_discharge_day, c.hospital_stay_days,
       c.sepsis_offset_sec, c.landmark_offset_sec,
       m.make30_death, m.make30_censor_offset_sec
FROM study_make30_final.sic_analysis_master m
JOIN study_make30_final.sic_core4_trajectory_cohort c USING (case_id)
LEFT JOIN study_make30_final.sic_case_extra e USING (case_id)
LEFT JOIN study_make30_final.sic_src_refs r
  ON r.ref_id = e.hospital_discharge_type_ref
ORDER BY m.case_id;
