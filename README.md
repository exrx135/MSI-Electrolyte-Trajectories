# Joint MSI–electrolyte trajectories in sepsis: analysis code

SQL, R and Python source code supporting the MSI–electrolyte trajectory study in MIMIC-IV, eICU-CRD and SICdb. This code snapshot corresponds to the JCC submission materials prepared on 22 October 2026.

## Contents

| Directory | Contents |
| --- | --- |
| `01_SQL` | Cohort, trajectory-variable, outcome, covariate and survival-data extraction |
| `02_Trajectory` | Data preparation; G2–G20 screening; structural comparison; multistart refinement; stability assessment; frozen scoring and parameter transport |
| `03_Analysis` | Three-model modified Poisson analyses, MICE, subgroup interactions, survival, RMST, incremental prediction, calibration, decision curves, organ-support analyses and missing-data sensitivity analyses |
| `04_Figures` | Final flowchart, trajectory, Kaplan–Meier and state-transition drawing scripts; related plots also reside beside their analysis code |
| `05_Tables` | Baseline, class-specific and supplementary-table builders |
| `docs` | Workflow, current figure/table mapping, input requirements and packaging verification |
| `tools` | Preparation of a separate private execution workspace |

See [CODE_INDEX.md](docs/CODE_INDEX.md) for every source file and [FIGURE_TABLE_MAP.md](docs/FIGURE_TABLE_MAP.md) for the current manuscript numbering.

## Data access

No original database extracts, patient-level derived data, imputed datasets, individual class assignments, fitted model objects or database credentials are distributed. Data must be obtained through the original databases under their access conditions. Requests concerning study-specific data or adjudication records should be directed to the author team; access remains subject to the original databases' rules.

Some historical standalone SQL files embedded private patient rows. In this public edition those rows are replaced by local input-table interfaces. Patient-specific adjudications in two analysis scripts are likewise loaded from a local CSV. Extraction and statistical rules are retained. See [PRIVATE_INPUTS.md](docs/PRIVATE_INPUTS.md) and [PUBLIC_COPY_CHANGES.json](docs/PUBLIC_COPY_CHANGES.json).

## Running the code

This is a curated source archive, not a data-inclusive one-command pipeline. It retains the executed analysis stages, versions, seeds, checks and necessary upstream dependencies. The `Vxxx` prefixes and some original figure numbers identify historical execution stages; the current numbering is in the figure/table map.

1. Obtain database access and prepare an authorised private working directory **outside this repository**.
2. Read [WORKFLOW.md](docs/WORKFLOW.md). Select the required stage; do not execute all scripts indiscriminately. Some historical wrappers install or resume specific checkpoints and some presentation helpers expect locally held templates.
3. Optionally run `python tools/prepare_private_workspace.py /absolute/path/to/a/new/private/workspace`. This copies source files to their recorded execution layout and translates the main project paths. It does not supply data or run an analysis. Review module-specific paths, R executable paths and font requirements before use.
4. Supply local inputs and run the stages in dependency order. Original study identity/hash checks are retained; they are checks against the study's frozen inputs, not universal checks for an arbitrary new cohort.

The statistical scripts use R (recorded execution version 4.5.2). Python is used for preparation, validation, plotting and Word-table generation. See [SOFTWARE.md](docs/SOFTWARE.md). Database extraction is PostgreSQL-specific.

The repository preparation involved source selection, removal of embedded patient data and static syntax checks. No statistical models or clinical database queries were rerun during packaging.

## 中文使用说明

本文件夹是供 GitHub 上传的代码包，按提取、建模、分析、绘图和制表整理。只上传这个文件夹，不要把整个 MSI-1 项目作为仓库。

- 不含原始数据、患者级处理数据、患者编号清单、插补数据或 RDS 模型对象。
- SQL 内嵌的患者数据已移除，改从授权环境中的本地表读取；原项目内的 SQL 和数据没有改动。
- 最新 RMST、增量预测、Model 3 口径敏感性分析、器官支持分析及新版 Figure 1 均已收录。
- 代码路径、需要的输入和主图/补充图表的对应关系在 `docs` 中。历史版本号用于追溯，不代表当前投稿图表编号。
- 本次未上传 GitHub、未重新拟合模型，也未改动稿件和现有分析结果。
