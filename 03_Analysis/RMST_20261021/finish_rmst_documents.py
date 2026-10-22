from pathlib import Path
import importlib.util, csv, json

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("table_builder", ROOT.parent / "Incremental_prediction_20261021/build_prediction_tables.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
mod.OUT = ROOT / "word_staging"
mod.OUT.mkdir(exist_ok=True)
with (ROOT/"RMST_table_ready.csv").open(encoding="utf-8-sig", newline="") as f:
    data = list(csv.DictReader(f))
def ci(s):
    return s.replace(" (", "\n(").replace(" to ", ", ").replace("-", "−")
values = [
    [r["cohort"], r["tau_days"],
     ci(r["C1_RMST_95CI"]),ci(r["C2_RMST_95CI"]),ci(r["C3_RMST_95CI"]),
     ci(r["C2_minus_C1_95CI"]),mod.pval(r["C2_vs_C1_P"]),
     ci(r["C3_minus_C1_95CI"]),mod.pval(r["C3_vs_C1_P"])]
    for r in data
]
mod.make_doc(
    "Table S8. Restricted mean in-hospital survival by trajectory class",
    "Unadjusted Kaplan–Meier estimates at 14 and 21 days after the 24-hour landmark",
    ["Cohort","Horizon\n(days)","C1 RMST\n(95% CI)","C2 RMST\n(95% CI)","C3 RMST\n(95% CI)",
     "C2 − C1\nDifference (95% CI)","P","C3 − C1\nDifference (95% CI)","P"],
    [2.5,1.5,3.3,3.3,3.3,4.4,1.2,4.4,1.2],
    values,
    [
        "RMST and differences are expressed in days. Follow-up starts 24 hours after sepsis onset, with live discharge or ICU day 30 treated as censoring, matching Figure S1. Fourteen days is the primary horizon; 21 days is the extended horizon. The same 9,706 MIMIC-IV, 5,804 eICU-CRD and 3,811 SICdb patients were used.",
        "RMST is the integral of the Kaplan–Meier in-hospital survival curve; variance uses the integrated Greenwood covariance. Differences are C2 or C3 minus C1, with normal-approximation 95% CIs and two-sided P values. Estimates are unadjusted and use the KM independent-censoring assumption. Negative differences indicate shorter restricted mean survival. RMST, restricted mean survival time; CI, confidence interval; ICU, intensive care unit.",
    ],
    "Table S8. Restricted mean in-hospital survival by trajectory class.docx",
)
wording = {
    "methods": "Unadjusted restricted mean survival time (RMST) was estimated by integrating the same Kaplan–Meier curves to 14 days after the landmark, with 21 days as an extended horizon. Class differences versus C1 were expressed in days with Greenwood-based 95% confidence intervals",
    "results": "RMST quantified the survival gradient in days. At 14 days after the landmark, C2 minus C1 differences were −1.73 days (95% CI −1.93 to −1.54) in MIMIC-IV, −1.97 (−2.23 to −1.70) in eICU-CRD and −1.52 (−1.82 to −1.22) in SICdb; C3 differences were −0.71, −0.55 and −0.60 days. At 21 days, C2 differences reached −2.96, −3.07 and −2.66 days, with all C2 and C3 intervals below zero (Figure S4; Table S8).",
    "conclusion": "RMST quantified the shorter survival associated with C2 and C3, while incremental analyses demonstrated prognostic information beyond conventional severity scores.",
    "legend": "Unadjusted restricted mean survival time (RMST) differences relative to C1, estimated from the same Kaplan–Meier data as Figure S1. The panels use horizons of 14 and 21 days after the 24-hour sepsis landmark. Points and bars show differences in days and 95% confidence intervals; negative values indicate shorter RMST. MIMIC-IV, eICU-CRD and SICdb are represented by circles, diamonds and triangles, respectively. Live discharge and the ICU-day-30 administrative cutoff are treated as censoring.",
}
(ROOT/"manuscript_wording.json").write_text(json.dumps(wording,ensure_ascii=False,indent=2),encoding="utf-8")
(ROOT/"Methods_and_interpretation.txt").write_text(
    "Restricted mean survival time extension\n\n" + wording["methods"] + ".\n\n" + wording["results"] +
    "\n\nAll analyses use the frozen KM cohorts and endpoints. Both horizons have observed risk-set support in every cohort and class. "
    "The estimated quantity is the integral of the KM in-hospital survival function under its live-discharge censoring assumption. "
    "It is not observed hospital length of stay or postdischarge all-cause survival. No adjustment or MAKE30 imputation was reused. "
    "The existing MAKE30, KM and trajectory results were not changed. Absolute RMSTs and integrated Greenwood variances were independently "
    "checked against survival::summary(survfit, rmean=...) for all 18 class/cohort/horizon combinations.\n",
    encoding="utf-8",
)
print("RMST Word table and manuscript wording saved.")
