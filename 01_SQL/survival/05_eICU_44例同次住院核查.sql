-- Read-only. Export ALL rows to 05_eICU_44例同次住院核查.csv.
-- Multiple rows per index stay are expected: all unit stays in the SAME hospital stay.
WITH target AS (
 SELECT * FROM public.patient WHERE patientunitstayid IN (SELECT patientunitstayid FROM local_inputs.eicu_survival_review_ids)
)
SELECT t.patientunitstayid AS index_patientunitstayid,
       t.uniquepid AS index_uniquepid,
       t.patienthealthsystemstayid AS index_hospitalstayid,
       t.hospitaladmitoffset AS index_hospitaladmitoffset,
       p.patientunitstayid AS related_patientunitstayid,
       p.uniquepid, p.patienthealthsystemstayid, p.hospitalid,
       p.unitvisitnumber, p.unitstaytype, p.unittype,
       p.hospitaladmitoffset, p.hospitaldischargeoffset,
       p.hospitaldischargestatus, p.hospitaldischargelocation,
       p.unitdischargeoffset, p.unitdischargestatus, p.unitdischargelocation,
       p.hospitaldischargeoffset - p.hospitaladmitoffset + t.hospitaladmitoffset
          AS hospital_discharge_index_icu_min,
       p.unitdischargeoffset - p.hospitaladmitoffset + t.hospitaladmitoffset
          AS unit_discharge_index_icu_min
FROM target t
JOIN public.patient p
 ON p.patienthealthsystemstayid = t.patienthealthsystemstayid
 AND p.uniquepid = t.uniquepid
ORDER BY t.patientunitstayid, p.unitvisitnumber, p.patientunitstayid;
