# Execution order and provenance

The source filenames preserve historical versions. `manifest.json` records the source's relative location, its execution location, the public-copy checksum and any privacy-related changes. `tools/prepare_private_workspace.py` prepares this execution layout in a new directory and normalises the archived-results path to `results`.

## 1. SQL and data preparation

Use the three `*_FINAL_*ELECTROLYTEQC.sql` files to extract cohort variables, trajectory measurements and MAKE30 components. Use the survival supplements and covariate queries when their upstream cohort records are ready. The eICU raw extraction precedes the existing time-alignment stage; do not substitute the raw export for the later, aligned analysis cohort.

`GBMTM_stage1_MIMIC_fix1.R` prepares/scales the development data. `GBTM_V267_stage1_full.R` provides the single-indicator trajectory dictionaries. `V316_STEP5_update_G3.R` is an original update/installation wrapper that expects its private update kit; it is retained for provenance, not as a universal first execution step. Its replay/transport dependencies are included.

## 2. Trajectory model development and selection

`GBMTM_V251_G5G20.R`, `GBMTM_V258_merge_scheduler.R` and `V258_lowG_screen.R` implement the joint-model screening and low-class-count extension. `V292_step2_visual_G2G20.R` and `V294_stage3_G2G20.R` assemble the screening/structural comparisons. `V297_joint_portrait_G3G6.R` and `build_level_shape_support.py` provide the level-versus-shape comparisons.

`V295_step4_2lane_G3G5.R` implements the multistart refinement. `V306_STEP4B_parallel_G3.R` and `V309_STEP4B_resume_G3.R` implement/resume the resampling stability work. These are sequential or alternative stage tools, not interchangeable model specifications.

## 3. Frozen scoring and external application

`V318_STEP5_score_G3.R` and the V320/V321/V322 scripts implement frozen classification and original cohort-specific evaluation. The V325 scoring/MIMIC/eICU/SICdb reproduction scripts are path-adjusted reproductions of this upstream work. Their earlier association models are **not** the final three-model MICE-adjusted analysis. Use the next stage for the reported Model 1–3 results.

## 4. Baseline covariates and final Model 1–3 analyses

Run the privately populated standalone covariate SQL, SICdb vitals supplement and eICU diagnosis/history extraction. `calculate_eicu_cci.py` computes the eICU Charlson score. `table1/code/build_table1.py` creates the analysis-ready baseline variables and table summaries.

In `关联分析_三模型_MICE_RR/code`, run `prepare_inputs.py`, then `run_models_final.R` (optionally one cohort argument: `MIMIC`, `eICU`, or `SICdb`), then `collect_and_validate.py` and `draw_forest.py`. The data-free `model_specification.json` defines the final covariates and imputation settings.

`补充表_亚组与交互作用/code/run_subgroups_safe.R` runs the final subgroup analyses by cohort, reading the main analysis inputs. `build_table.py` and `export_word.py` create the subgroup outputs. The approved concise table note is supplied as `05_Tables/table_legend.txt`.

## 5. Survival and RMST

`Figure3_三库院内死亡累积发生.R` prepares the adjudicated endpoint and related summaries; it requires the private SICdb adjudication file and eICU status-review export. `Figure3_KM_logrank.R` computes the Kaplan–Meier results. The final `figure3/draw_figure3.py` renders them.

In `RMST_20261021`, `prepare_rmst_inputs.py` stages the frozen eligible survival data, `run_rmst.R` runs RMST, and `finish_rmst_documents.py` formats its tables using the neighbouring incremental-prediction table helper. RMST and KM use the same survival cohort.

## 6. Incremental prediction

In `Incremental_prediction_20261021`, run `prepare_prediction_inputs.py`, `run_incremental_prediction.R`, `finish_prediction_outputs.R` and `build_prediction_tables.py`. These implement the score comparisons, repeated cross-validation, calibration, decision curves and tables. Their score-complete samples differ from the full MAKE30 association samples; the existing coverage checks are retained.

## 7. Organ support and state transitions

`FigureS5_organ_support/analyse_and_plot.R` generates the organ-support analysis and the aligned alluvial figure. It accepts a private project root, output directory and optional `--plot-only`; its current main-text number is Figure 6. `build_outputs.py` is its historical presentation helper and expects local Word templates.

The upstream V327/V330 scripts generate trajectory/state summary data. The final Figure 2, Figure 5 and external-state drawing scripts apply the selected layout and palette. Their original output names are mapped to current submission numbers in `FIGURE_TABLE_MAP.md`.

## 8. Missing-data sensitivity analysis

The V331 scripts construct and diagnose the trajectory imputations and MAR/MNAR scenarios. Their earlier adjusted association outputs are retained only as upstream context. **Use `S4_Model3_20261019/run_s4_model3.R` and `verify_s4.R` for the current sensitivity associations using the main Model 3 adjustment set.** The current supplementary table is Table S12.

## 9. Presentation

The latest Figure 1 script excludes neither missing MAKE30 outcomes nor survival-specific endpoint issues from the trajectory cohort. The baseline aggregation/table scripts add the class-specific summaries and 30-day in-hospital mortality.

Historical document builders preserve their original output numbers and may require existing private presentation templates. Use `FIGURE_TABLE_MAP.md` for final numbering. Document-only manuscript revisions, reference lookups, redline generation and directory-cleanup scripts are not analysis code and are not included.

The package has been statically checked. Data-dependent stages were not rerun during repository preparation; source packaging does not establish an end-to-end reproducibility run on a new system.
