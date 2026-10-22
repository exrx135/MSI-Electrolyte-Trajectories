from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT = Path(r"Z:/Research/MSI-1")
OUT = ROOT / "analysis_updates/Incremental_prediction_20261021"
PRIVATE = OUT / "private_patient_level"
PRIVATE.mkdir(parents=True, exist_ok=True)
SOURCE = ROOT / "历史归档/results_整理前_20261019/results/table1/derived"
frozen = {
    "MIMIC": ("MIMIC_预后关联/V319_MIMIC_analysis_cohort_G3.csv", "stay_id", 9758, 2442),
    "eICU": ("eICU_外部验证/V321_analysis_cohort_G3.csv", "patientunitstayid", 5709, 2060),
    "SICdb": ("SICdb_外部验证/V322_analysis_cohort_G3.csv", "case_id", 3635, 949),
}
coverage, manifests, distributions = [], [], []
for db, (rel, idcol, expected_n, expected_events) in frozen.items():
    source = SOURCE / f"{db}_Table1_analysis.csv"
    ref = ROOT / "STEP5_冻结分型与预后外部验证" / rel
    d = pd.read_csv(source)
    d = d.loc[d.make30_analysis.notna()].copy()
    f = pd.read_csv(ref, low_memory=False)
    f = f.loc[f.make30.notna()].copy()
    assert len(d) == expected_n and int(d.make30_analysis.sum()) == expected_events
    assert d.cohort_unit_id.is_unique and f[idcol].is_unique
    z = d[["cohort_unit_id", "assigned_class", "make30_analysis"]].merge(
        f[[idcol, "assigned_class", "make30"]],
        left_on="cohort_unit_id", right_on=idcol, suffixes=("_new", "_frozen"), validate="one_to_one")
    assert len(z) == expected_n
    assert (z.assigned_class_new == z.assigned_class_frozen).all()
    assert (z.make30_analysis == z.make30).all()
    score = "saps3_analysis" if db == "SICdb" else "aps_score_analysis"
    p = d[["cohort_unit_id", "assigned_class", "make30_analysis", "sepsis_sofa_score_analysis", score]].copy()
    p.columns = ["id", "class", "y", "sofa", "additional_score"]
    for col in ["sofa", "additional_score"]:
        assert (p[col].dropna() >= 0).all()
        distributions.append(dict(dataset=db, score=col, available=int(p[col].notna().sum()),
                                  minimum=p[col].min(), maximum=p[col].max(),
                                  median=p[col].median()))
    assert p.sofa.notna().all()
    cc = p[["sofa", "additional_score"]].notna().all(axis=1)
    n, e = int(cc.sum()), int(p.loc[cc, "y"].sum())
    p.loc[cc].sort_values("id").to_csv(PRIVATE / f"{db}_complete_score_input.csv", index=False)
    if db == "eICU":
        p.sort_values("id").to_csv(PRIVATE / "eICU_full_SOFA_input.csv", index=False)
    coverage.append(dict(dataset=db, make30_cohort_n=len(p), make30_events=int(p.y.sum()),
                         sofa_available=int(p.sofa.notna().sum()), additional_score_available=n,
                         additional_score_label="SAPS 3" if db == "SICdb" else "APS",
                         complete_score_n=n, complete_score_events=e, coverage_percent=100 * cc.mean(),
                         excluded_missing_additional_score=int((~cc).sum()),
                         excluded_events=int(p.loc[~cc, "y"].sum()),
                         ids_classes_outcomes_verified=True))
    for path in [source, ref]:
        manifests.append({"dataset": db, "source": str(path),
                          "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
pd.DataFrame(coverage).to_csv(OUT / "Score_coverage_and_analysis_denominators.csv", index=False)
pd.DataFrame(distributions).to_csv(OUT / "Observed_score_ranges.csv", index=False)
(OUT / "Input_source_manifest.json").write_text(json.dumps(manifests, ensure_ascii=False, indent=2), encoding="utf-8")
print(pd.DataFrame(coverage).to_string(index=False))
