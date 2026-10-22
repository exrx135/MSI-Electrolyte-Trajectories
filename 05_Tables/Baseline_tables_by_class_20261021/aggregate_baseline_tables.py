"""Aggregate the existing final Table 1 data without refitting any model.

This script preserves the frozen Table 1 patient-level processing. It adds a
separately audited ICU-admission-to-day-30 in-hospital mortality description.
Only aggregate outputs are written; no patient-level output is exported.
"""
from pathlib import Path
import json
import hashlib
import math
import numpy as np
import pandas as pd

ROOT = Path("Z:/Research/MSI-1")
SOURCE = ROOT / "历史归档/results_整理前_20261019/results"
OUT = ROOT / "analysis_updates/Baseline_tables_by_class_20261021"
OUT.mkdir(parents=True, exist_ok=True)
SURV = SOURCE / "生存分析_补充数据导出"
LAYOUT_PATH = SOURCE / "table1/derived/Table1_rows.json"
LAYOUT = json.loads(LAYOUT_PATH.read_text(encoding="utf-8-sig"))
DBS = {"MIMIC": 9804, "eICU": 5839, "SICdb": 3831}
TWO_DIGITS = {
    "albumin", "creatinine", "baseline_creatinine", "potassium",
    "magnesium", "phosphate", "inr", "msi_m1", "msi_window_mean"
}
INPUTS = [LAYOUT_PATH]
MORTALITY_NOTES = {}
MORTALITY_COMPARISONS = []


def read_csv(path):
    INPUTS.append(path)
    return pd.read_csv(path)


def finite_or_none(v):
    return None if pd.isna(v) or not math.isfinite(float(v)) else float(v)


def audit_summary(df, value_column):
    return {
        "n_total": int(len(df)),
        "n_observed": int(df[value_column].notna().sum()),
        "n_unknown": int(df[value_column].isna().sum()),
        "events": int(df[value_column].sum()),
        "percent_observed": float(100 * df[value_column].mean()),
        "by_class": {
            f"C{int(c)}": {
                "n_total": int(len(g)),
                "n_observed": int(g[value_column].notna().sum()),
                "n_unknown": int(g[value_column].isna().sum()),
                "events": int(g[value_column].sum()),
                "percent_observed": float(100 * g[value_column].mean()),
            } for c, g in df.groupby("assigned_class")
        }
    }


def add_mortality(db, dat):
    """Exact-time result, plus an explicitly separate binary adjudication option."""
    if db == "MIMIC":
        x = read_csv(SURV / "01_MIMIC_生存补导出.csv")
        x = x[x.stay_id.isin(dat.cohort_unit_id)].copy()
        for c in ["icu_intime", "frozen_death_time",
                  "admission_discharge_time", "administrative_end_time"]:
            x[c] = pd.to_datetime(x[c], format="mixed", dayfirst=True, errors="coerce")
        assert (x.administrative_end_time - x.icu_intime).eq(pd.Timedelta(days=30)).all()
        y = pd.Series(np.nan, index=x.index)
        y.loc[x.hospital_expire_flag.eq(0)] = 0
        dead = x.hospital_expire_flag.eq(1)
        valid = dead & x.frozen_death_time.notna() & x.frozen_death_time.ge(x.icu_intime)
        y.loc[valid] = x.loc[valid, "frozen_death_time"].le(
            x.loc[valid, "administrative_end_time"]).astype(int)
        discordant = dead & x.frozen_death_time.gt(x.admission_discharge_time)
        gap_hours = (
            x.loc[discordant, "frozen_death_time"]
            - x.loc[discordant, "admission_discharge_time"]
        ).dt.total_seconds() / 3600
        MORTALITY_NOTES[db] = {
            "endpoint": "In-hospital death from ICU admission through 30 days",
            "clock": "ICU intime to ICU intime + 30 days, inclusive",
            "status": "hospital_expire_flag; hospital-status deaths use frozen_death_time",
            "death_time_after_discharge_n": int(discordant.sum()),
            "death_time_after_discharge_and_day30_event_n": int(y[discordant].sum()),
            "death_time_discharge_conflict_crosses_day30_boundary_n": int(
                (discordant & x.admission_discharge_time.le(x.administrative_end_time)
                 & x.frozen_death_time.gt(x.administrative_end_time)).sum()),
            "death_time_minus_discharge_time_hours": {
                "min": float(gap_hours.min()), "q1": float(gap_hours.quantile(.25)),
                "median": float(gap_hours.median()), "q3": float(gap_hours.quantile(.75)),
                "max": float(gap_hours.max()),
                "within_24_hours_n": int(gap_hours.le(24).sum()),
            },
            "interpretation": (
                "The 97 known hospital deaths with timestamp after the discharge "
                "timestamp are not excluded from this binary endpoint: hospital death "
                "status is known and their day-30 classification is not changed by "
                "either timestamp. These records were excluded from the prior KM "
                "time-to-event analysis. No KM data or results are modified."
            ),
        }
        identifier = "stay_id"
        y_augmented = y.copy()
    elif db == "eICU":
        x = read_csv(SURV / "02_eICU_生存补导出.csv")
        x = x[x.patientunitstayid.isin(dat.cohort_unit_id)].copy()
        q = read_csv(SURV / "05_eICU_44例同次住院核查.csv")
        status = x.hospitaldischargestatus.copy()
        death = pd.to_numeric(x.death_icu_offset_min_from_status, errors="coerce")
        overrides = {"expired": 0, "alive": 0}
        for idx in x.index[~status.isin(["Alive", "Expired"])]:
            z = q[q.index_patientunitstayid.eq(x.at[idx, "patientunitstayid"])]
            d = z[z.unitdischargestatus.eq("Expired") & z.unitdischargelocation.eq("Death")]
            a = z[z.unitdischargestatus.eq("Alive")
                  & z.unitdischargelocation.isin(["Home", "Other Hospital", "OtherHospital"])]
            if len(d):
                status.loc[idx] = "Expired"
                death.loc[idx] = pd.to_numeric(d.unit_discharge_index_icu_min).min()
                overrides["expired"] += 1
            elif len(a):
                status.loc[idx] = "Alive"
                overrides["alive"] += 1
        y = pd.Series(np.nan, index=x.index)
        y.loc[status.eq("Alive")] = 0
        valid = status.eq("Expired") & death.notna() & death.ge(0)
        y.loc[valid] = death.loc[valid].le(43200).astype(int)
        assert overrides == {"expired": 5, "alive": 4}, overrides
        MORTALITY_NOTES[db] = {
            "endpoint": "In-hospital death from ICU admission through 30 days",
            "clock": "Index ICU minute 0 through 43,200, inclusive",
            "status": "Hospital discharge status with the existing same-hospital-stay adjudications",
            "same_stay_resolved_as_dead_n": overrides["expired"],
            "same_stay_resolved_as_alive_n": overrides["alive"],
            "remaining_unknown_hospital_status_n": int((~status.isin(["Alive", "Expired"])).sum()),
            "interpretation": "Unknown hospital discharge status is not assumed to be survival.",
        }
        identifier = "patientunitstayid"
        y_augmented = y.copy()
    else:
        x = read_csv(SURV / "03_SICdb_生存补导出.csv")
        x = x[x.case_id.isin(dat.cohort_unit_id)].copy()
        q = read_csv(SURV / "04_SICdb_补充住院出院状态.csv")
        x = x.merge(q[["case_id", "hospital_discharge_type_ref"]], on="case_id",
                    how="left", validate="one_to_one")
        status = x.hospital_discharge_type_ref.map(
            {3130: "Expired", 3129: "Alive", 3131: "Alive", 3132: "Alive"})
        # Individual adjudications are loaded privately, never distributed with code.
        import os
        adjudication_path = Path(os.environ.get("MSI_SICDB_ADJUDICATIONS", str(ROOT / "private_inputs/sicdb_adjudications.csv")))
        adjudications = pd.read_csv(adjudication_path)
        required = {"case_id", "death_override", "discharge_override", "ignore_negative_discharge_day"}
        assert required.issubset(adjudications.columns) and adjudications.case_id.is_unique
        death_overrides = adjudications.loc[adjudications.death_override.eq(1), "case_id"].tolist()
        alive_overrides = adjudications.loc[adjudications.discharge_override.eq(1), "case_id"].tolist()
        assert len(death_overrides) == 12 and len(alive_overrides) == 6
        # Existing adjudications identified these by case_id.
        status.loc[x.case_id.isin(death_overrides)] = "Expired"
        status.loc[x.case_id.isin(alive_overrides)] = "Alive"
        assert x.case_id.isin(death_overrides).sum() == 12
        assert x.case_id.isin(alive_overrides).sum() == 6
        y = pd.Series(np.nan, index=x.index)
        y.loc[status.eq("Alive")] = 0
        valid = (status.eq("Expired") & x.offset_of_death_sec.notna()
                 & x.offset_of_death_sec.ge(x.icu_offset_sec))
        y.loc[valid] = x.loc[valid, "offset_of_death_sec"].le(
            x.loc[valid, "administrative_end_offset_sec"]).astype(int)
        assert x.administrative_end_offset_sec.sub(x.icu_offset_sec).eq(30 * 86400).all()
        y_augmented = y.copy()
        # Binary adjudication authorized by the main task; no KM event time is assigned.
        # Known ICU observation beyond day 30 proves no death within day 30.
        observed_beyond30 = y_augmented.isna() & x.time_of_stay_sec.gt(
            x.administrative_end_offset_sec)
        y_augmented.loc[observed_beyond30] = 0
        # Two status-confirmed deaths with no exact death timestamp have day-level
        # discharge fields that agree on which side of day 30 death occurred.
        # The task explicitly authorized these two day-level adjudications.
        missing_time = status.eq("Expired") & x.offset_of_death_sec.isna()
        day_a = x.hospital_discharge_day
        day_b = x.hospital_stay_days
        within_a = day_a.mul(86400).le(x.administrative_end_offset_sec)
        within_b = day_b.mul(86400).le(x.administrative_end_offset_sec)
        unambiguous_day_bin = missing_time & day_a.notna() & day_b.notna() & within_a.eq(within_b)
        y_augmented.loc[unambiguous_day_bin] = within_a.loc[unambiguous_day_bin].astype(int)
        MORTALITY_NOTES[db] = {
            "endpoint": "In-hospital death from ICU admission through 30 days",
            "clock": "icu_offset_sec through icu_offset_sec + 2,592,000, inclusive",
            "status": "Supplementary hospital discharge status and the already agreed 12-death / 6-alive adjudications",
            "death_at_or_before_landmark_included_n": int(
                (valid & x.offset_of_death_sec.le(x.landmark_offset_sec) & y.eq(1)).sum()),
            "unknown_hospital_status_n": int(status.isna().sum()),
            "known_dead_missing_exact_death_time_n": int(missing_time.sum()),
            "known_ICU_observation_beyond_day30_n": int(observed_beyond30.sum()),
            "day_level_resolved_death_n": int(unambiguous_day_bin.sum()),
            "day_level_evidence": [
                {
                    "case_id": int(z.case_id),
                    "hospital_discharge_day": float(z.hospital_discharge_day),
                    "hospital_stay_days": float(z.hospital_stay_days),
                    "icu_offset_sec": float(z.icu_offset_sec),
                    "classification": int(y_augmented.loc[i]),
                }
                for i, z in x.loc[unambiguous_day_bin].iterrows()
            ],
            "known_observation_evidence": [
                {
                    "case_id": int(z.case_id),
                    "icu_observation_days": float(
                        (z.time_of_stay_sec - z.icu_offset_sec) / 86400),
                    "classification": int(y_augmented.loc[i]),
                }
                for i, z in x.loc[observed_beyond30].iterrows()
            ],
            "interpretation": (
                "Exact-time endpoint leaves five unknown hospital-status records and "
                "two known hospital deaths with unavailable exact time as unknown. "
                "The final binary adjudication resolves one unknown-status "
                "case by ICU observation beyond day 30 and two deaths by concordant "
                "day-level discharge timing (one within 30 days, one after 30 days). "
                "Neither rule creates or changes any KM death time."
            ),
        }
        identifier = "case_id"
    x["mortality_icu_day30_exact_analysis"] = y
    x["mortality_icu_day30_adjudicated_analysis"] = y_augmented
    x = x.rename(columns={identifier: "cohort_unit_id"})
    ans = dat.merge(x[["cohort_unit_id", "mortality_icu_day30_exact_analysis",
                      "mortality_icu_day30_adjudicated_analysis"]],
                   on="cohort_unit_id", how="left", validate="one_to_one")
    ans["mortality_icu_day30_analysis"] = ans.mortality_icu_day30_adjudicated_analysis
    old = read_csv(SOURCE / f"生存分析_三库图与代码/{db}_逐患者纳排审计.csv")
    check = ans.merge(old[["unit_id", "event", "eligible", "reason"]],
                      left_on="cohort_unit_id", right_on="unit_id", validate="one_to_one")
    kept = check[check.eligible.astype(str).str.lower().eq("true")]
    mismatch = kept.mortality_icu_day30_exact_analysis.ne(kept.event.eq(1).astype(int))
    assert not mismatch.any(), (db, int(mismatch.sum()))
    for key in ["mortality_icu_day30_exact_analysis", "mortality_icu_day30_adjudicated_analysis"]:
        MORTALITY_COMPARISONS.append({
            "database": db, "definition": key, **audit_summary(ans, key)})
    MORTALITY_NOTES[db]["old_KM_eligible_exact_endpoint_mismatches"] = int(mismatch.sum())
    MORTALITY_NOTES[db]["old_KM_n"] = int(len(kept))
    MORTALITY_NOTES[db]["old_KM_death_events"] = int(kept.event.eq(1).sum())
    old_component = ans.make30_death_analysis
    new_endpoint = ans.mortality_icu_day30_analysis
    both_known = old_component.notna() & new_endpoint.notna()
    MORTALITY_NOTES[db]["frozen_MAKE30_death_component_comparison"] = {
        "old_component_observed_n": int(old_component.notna().sum()),
        "old_component_events": int(old_component.sum()),
        "old_component_missing_n": int(old_component.isna().sum()),
        "both_endpoints_known_n": int(both_known.sum()),
        "new_death_old_component_zero_n": int(
            (both_known & new_endpoint.eq(1) & old_component.eq(0)).sum()),
        "new_no_death_old_component_one_n": int(
            (both_known & new_endpoint.eq(0) & old_component.eq(1)).sum()),
        "new_death_old_component_missing_n": int(
            (old_component.isna() & new_endpoint.eq(1)).sum()),
        "interpretation": (
            "The new descriptive endpoint and the frozen MAKE30 death component "
            "are retained as different variables. No frozen MAKE30 values or "
            "regression inputs are changed."
        ),
    }
    return ans


def summarize(db, data, row):
    variable = row["variable"]
    field = variable + "_analysis"
    structural_na = row.get(db) == "N.A."
    v = pd.to_numeric(data[field], errors="coerce")
    if structural_na:
        return {"display": "N.A.", "n_total": int(len(data)), "n_observed": 0,
                "n_missing": int(len(data)), "structural_unavailability": True,
                "events": None, "percent": None, "q1": None, "median": None, "q3": None}
    n = int(v.notna().sum())
    info = {"n_total": int(len(data)), "n_observed": n,
            "n_missing": int(v.isna().sum()), "structural_unavailability": False}
    binary = (variable == "male" or row["section"].startswith("Comorbidities")
              or row["section"].startswith("Treatment")
              or row["section"].startswith("Outcomes"))
    if binary:
        events = int(v.eq(1).sum())
        pct = 100 * events / n if n else np.nan
        if variable.startswith("mortality_icu_day30"):
            display = f"{events:,}/{n:,} ({pct:.1f})" if n else "N.A."
        else:
            display = f"{events} ({pct:.1f})" if n else "N.A."
        info.update(display=display, events=events, percent=finite_or_none(pct),
                    q1=None, median=None, q3=None)
    else:
        q1, median, q3 = (v.quantile([.25, .5, .75]).to_list() if n else [np.nan] * 3)
        digits = 2 if variable in TWO_DIGITS else 1
        display = (f"{median:.{digits}f} [{q1:.{digits}f}, {q3:.{digits}f}]"
                   if n else "N.A.")
        info.update(display=display, events=None, percent=None,
                    q1=finite_or_none(q1), median=finite_or_none(median), q3=finite_or_none(q3))
    return info


def main():
    output = {"description": "Final Table 1 processing; full frozen trajectory cohorts",
              "cohorts": {}, "overall_rows": [],
              "mortality_definition": "adjudicated_binary_in_hospital_death_within_ICU_day30",
              "notes": [
                  "Continuous data are median [Q1, Q3]; categorical data are n (%).",
                  "Class percentages use observed denominators for each variable.",
                  "No new cleaning, imputation, model fitting or SQL was performed.",
                  "All original Table 1 variables, labels and structural N.A. were preserved.",
                  "MSI 0–2 h maximum is distinct from mean of observed two-hour maxima.",
                  "Mortality is in-hospital mortality, not all-cause post-discharge mortality.",
              ]}
    rows = []
    for r in LAYOUT:
        r = r.copy()
        if r["section"].startswith("Outcomes"):
            r["section"] = "Outcomes"
        rows.append(r)
    rows.append({"section": "Outcomes", "variable": "mortality_icu_day30",
                 "label": "30-day in-hospital mortality after ICU admission"})
    aggregate = []
    discrepancies = []
    for db, n in DBS.items():
        dat = read_csv(SOURCE / f"table1/derived/{db}_Table1_analysis.csv")
        assert len(dat) == n and dat.cohort_unit_id.is_unique
        dat = add_mortality(db, dat)
        block = {"n_total": n,
                 "class_n": {f"C{int(k)}": int(v) for k, v in dat.assigned_class.value_counts().sort_index().items()},
                 "rows": []}
        for r in rows:
            vals = {"Overall": summarize(db, dat, r)}
            vals.update({f"C{int(c)}": summarize(db, g, r)
                         for c, g in dat.groupby("assigned_class")})
            block["rows"].append({"section": r["section"], "variable": r["variable"],
                                  "label": r["label"], "values": vals})
            for stratum, details in vals.items():
                aggregate.append({"database": db, "stratum": stratum, "section": r["section"],
                                  "variable": r["variable"], "label": r["label"], **details})
            if r.get(db) is not None and vals["Overall"]["display"] != r[db]:
                discrepancies.append({"database": db, "variable": r["variable"],
                                      "original": r[db], "reproduced": vals["Overall"]["display"]})
        # Exact timestamp audit explicitly separate; final table uses final adjudication.
        alt_row = {"section": "Outcomes", "variable": "mortality_icu_day30_exact",
                   "label": "30-day in-hospital mortality after ICU admission"}
        block["mortality_exact_timestamp_audit"] = {
            "Overall": summarize(db, dat, alt_row),
            **{f"C{int(c)}": summarize(db, g, alt_row) for c, g in dat.groupby("assigned_class")}}
        output["cohorts"][db] = block
        pd.DataFrame([
            {"section": a["section"], "variable": a["variable"], "Variable": a["label"],
             "Overall": a["values"]["Overall"]["display"],
             **{c: a["values"][c]["display"] for c in ["C1", "C2", "C3"]}}
            for a in block["rows"]
        ]).to_csv(OUT / f"{db}_baseline_by_class_display.csv", index=False, encoding="utf-8-sig")
    # Numeric quantile formatting is required to match the actual final Table 1.
    if discrepancies:
        raise AssertionError({"original_Table1_display_mismatches": discrepancies})
    for i, row in enumerate(rows):
        output["overall_rows"].append({
            "section": row["section"], "variable": row["variable"], "label": row["label"],
            **{db: output["cohorts"][db]["rows"][i]["values"]["Overall"] for db in DBS}})
    (OUT / "baseline_tables_by_class.json").write_text(
        json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    audit = pd.DataFrame(aggregate)
    audit.to_csv(OUT / "baseline_summaries_numeric_long.csv", index=False, encoding="utf-8-sig")
    msi = audit[audit.variable.isin(
        ["msi_m1", "msi_window_mean", "sodium", "potassium", "chloride", "bicarbonate"])]
    msi.to_csv(OUT / "MSI_and_electrolytes_overall_by_class.csv", index=False, encoding="utf-8-sig")
    (OUT / "mortality_endpoint_audit.json").write_text(
        json.dumps({"notes": MORTALITY_NOTES, "summaries": MORTALITY_COMPARISONS},
                   ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    mortality_cells = {
        db: output["cohorts"][db]["rows"][-1]["values"]["Overall"] for db in DBS
    }
    p_reference = mortality_cells["MIMIC"]["events"] / mortality_cells["MIMIC"]["n_observed"]
    s1_rows = []
    for db, cell in mortality_cells.items():
        p = cell["events"] / cell["n_observed"]
        denom = math.sqrt((p * (1 - p) + p_reference * (1 - p_reference)) / 2)
        smd = abs(p - p_reference) / denom if denom else 0.0
        s1_rows.append({
            "database": db,
            "variable": "mortality_icu_day30",
            "label": "30-day in-hospital mortality after ICU admission",
            "events": cell["events"], "n_observed": cell["n_observed"],
            "display": cell["display"], "n_total": cell["n_total"],
            "n_unknown": cell["n_missing"],
            "unknown_pct": 100 * cell["n_missing"] / cell["n_total"],
            "unknown_display": f'{cell["n_missing"]:,}/{cell["n_total"]:,} ({100*cell["n_missing"]/cell["n_total"]:.1f}%)',
            "binary_smd_vs_MIMIC": None if db == "MIMIC" else smd,
            "binary_smd_display": "—" if db == "MIMIC" else f"{smd:.3f}",
        })
    pd.DataFrame(s1_rows).to_csv(OUT / "new_mortality_row_for_Table1_and_S1.csv",
                               index=False, encoding="utf-8-sig")
    (OUT / "new_mortality_row_for_Table1_and_S1.json").write_text(
        json.dumps(s1_rows, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    provenance = [{"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                  for p in dict.fromkeys(INPUTS)]
    (OUT / "input_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    readme = """Baseline tables by trajectory class: aggregate data handoff
==========================================================

Scope: MIMIC-IV 9,804; eICU-CRD 5,839; SICdb 3,831.
Patient-level processing is exactly the final Table 1 processing. No SQL,
trajectory fitting, imputation, regression, MAKE30 or KM result was rerun or changed.
Recomputed original Table 1 overall display cells matched all original cells.

Primary file: baseline_tables_by_class.json
- cohorts[database].class_n: frozen C1/C2/C3 counts.
- cohorts[database].rows: ordered section, variable, label, values.
- values[Overall/C1/C2/C3]: display, n_total, n_observed, n_missing, events,
  percent, q1, median, q3, structural_unavailability.
- overall_rows: same rows with MIMIC/eICU/SICdb overall cells.
- The formal new mortality variable is mortality_icu_day30.
- Exact-timestamp-only mortality is an internal audit, not the formal new row.

Display CSVs: MIMIC_baseline_by_class_display.csv, eICU_baseline_by_class_display.csv,
SICdb_baseline_by_class_display.csv.
New Table 1 / Table S1 row: new_mortality_row_for_Table1_and_S1.json and .csv.
Numeric audit: baseline_summaries_numeric_long.csv.
MSI and electrolyte writing values: MSI_and_electrolytes_overall_by_class.csv.

Mortality description
---------------------
30-day in-hospital mortality after ICU admission uses adjudicated hospital
discharge status and recorded timing, regardless of the old KM landmark inclusion.
MIMIC 1,632/9,804 (16.6%); eICU 1,188/5,804 (20.5%);
SICdb 630/3,827 (16.5%).
Unknown outcomes: 0, 35 and 4, respectively. The table displays events/observed N (%).
These are descriptive proportions, not day-30 Kaplan-Meier survival estimates.
SICdb two missing exact death times are adjudicated only as within/after day 30
from known hospital death status and concordant discharge-day evidence.
No exact death timestamps are created. One unknown-discharge-status patient
has ICU observation extending beyond day 30 and therefore no death by day 30.

MIMIC 97 deaths have death time later than the discharge time, all by less than
24 hours: minimum 1 minute, median 20.42 h, IQR 16.75-21.78 h, maximum 23 h 59 min.
The two timestamps agree on day-30 classification for all 97: 93 inside, 4 outside.
These add 93 events relative to the existing KM eligible-event count. The old KM
analysis remains unchanged. Details and frozen-component comparison are in
mortality_endpoint_audit.json.

MSI interpretation
------------------
msi_m1 is the maximum in the first observed 0-2 h window.
msi_window_mean is each patient's mean across observed two-hour window maxima
within the first 24 h; it is not a mean of all raw bedside MSI readings.
The cohort values are medians [Q1, Q3] across patients. The four electrolyte
rows use the existing final Table 1 first-24-hour mean values.

Suggested short new mortality footnote
-------------------------------------
30-day in-hospital mortality was assessed from ICU admission using separately
adjudicated hospital status and timing and is reported as events/evaluable patients (%).

Input paths and SHA-256 hashes are recorded in input_provenance.json.
"""
    (OUT / "README_evidence_and_schema.txt").write_text(readme, encoding="utf-8")
    print(json.dumps({"output": str(OUT), "Table1_display_mismatches": discrepancies,
                      "mortality": MORTALITY_COMPARISONS}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
