"""Copy frozen KM inputs unchanged into an analysis-only ASCII-path workspace."""
from pathlib import Path
import hashlib
import json
import shutil
import pandas as pd

OUT = Path(__file__).resolve().parent
SOURCE = Path("Z:/Research/MSI-1/历史归档/results_整理前_20261019/results/生存分析_三库图与代码")
PRIVATE = OUT / "private_patient_level"
PRIVATE.mkdir(parents=True, exist_ok=True)
records = []
for cohort in ("MIMIC", "eICU", "SICdb"):
    source = SOURCE / f"{cohort}_KM分析数据.csv"
    target = PRIVATE / f"{cohort}_frozen_KM.csv"
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(source, target)
    assert hashlib.sha256(target.read_bytes()).hexdigest() == source_hash
    data = pd.read_csv(target)
    assert data.unit_id.is_unique
    assert data.eligible.all()
    assert (data.time_days > 0).all()
    assert set(data.event.unique()).issubset({0, 1, 2})
    assert (data.death == (data.event == 1).astype(int)).all()
    assert set(data.phenotype.unique()) == {"C1", "C2", "C3"}
    records.append({
        "cohort": cohort,
        "source": str(source),
        "working_copy": str(target),
        "sha256": source_hash,
        "n": len(data),
        "deaths": int(data.death.sum()),
        "live_discharge_censored": int((data.event == 2).sum()),
        "administratively_censored": int((data.event == 0).sum())
    })
(OUT / "Input_source_manifest.json").write_text(
    json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(records, ensure_ascii=False, indent=2))
