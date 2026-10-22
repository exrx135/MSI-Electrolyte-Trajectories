# Current JCC figure and table mapping

These numbers refer to the current JCC submission, not the historical numbers embedded in some source filenames.

| Submission item | Main source / generation chain |
| --- | --- |
| Figure 1. Study flowchart | `04_Figures/figure1/draw_figure1.py` — latest full-trajectory cohort version |
| Figure 2. MSI and electrolyte trajectories across three cohorts | Upstream `V330_STEP6_publication_figures_G3.R`; `04_Figures/figure2/draw_figure2.py` |
| Figure 3. Kaplan–Meier survival curves across three cohorts | `03_Analysis/upstream/Figure3_三库院内死亡累积发生.R`; `Figure3_KM_logrank.R`; `04_Figures/figure3/draw_figure3.py` |
| Figure 4. Associations of trajectory classes with MAKE30 | `03_Analysis/关联分析_三模型_MICE_RR/code/run_models_final.R`, `collect_and_validate.py`, `draw_forest.py` |
| Figure 5. Within-class changes in MSI and electrolyte states in MIMIC-IV | Upstream state summaries; `04_Figures/figure5/draw_figure5.py` |
| Figure 6. Trajectory classes, early organ support and MAKE30 | `03_Analysis/FigureS5_organ_support/analyse_and_plot.R` — aligned final alluvial layout |
| Figure S1. Internally validated calibration for MAKE30 | `03_Analysis/Incremental_prediction_20261021/run_incremental_prediction.R`, `finish_prediction_outputs.R` |
| Figure S2. Internally validated decision curves for MAKE30 | Same incremental-prediction scripts |
| Figure S3. Restricted mean in-hospital survival differences | `03_Analysis/RMST_20261021/run_rmst.R` |
| Figure S4. Within-class changes in MSI and electrolyte states in eICU-CRD | `05_Tables/补充材料_20261019/代码/draw_external_states.py eICU` |
| Figure S5. Within-class changes in MSI and electrolyte states in SICdb | Same script with `SICdb` argument |
| Table 1. Baseline characteristics of the three study cohorts | `05_Tables/table1/code/build_table1.py`, `render_table1.py`; baseline aggregation/table builders below |
| Table 2. Subgroup associations of trajectory classes with MAKE30 in MIMIC-IV | `03_Analysis/补充表_亚组与交互作用/code/run_subgroups_safe.R`, `build_table.py`, `export_word.py` |
| Table S1. Incremental discrimination and reclassification for MAKE30 | Incremental-prediction analysis and `build_prediction_tables.py` |
| Table S2. Internally validated performance of MAKE30 risk models | Same incremental-prediction chain |
| Table S3. Between-cohort standardised differences and missingness of baseline characteristics | `05_Tables/table1/code/build_table1.py`; baseline table builders |
| Tables S4–S6. Baseline characteristics by trajectory class in MIMIC-IV, eICU-CRD and SICdb | `05_Tables/Baseline_tables_by_class_20261021/aggregate_baseline_tables.py`, `build_baseline_word_tables.py` |
| Table S7. MAKE30 and component outcomes by trajectory class | `05_Tables/候选补充材料_20261018/代码/build_candidate_tables.py` and subsequent supplementary-table formatting |
| Table S8. Restricted mean in-hospital survival by trajectory class | RMST analysis and `finish_rmst_documents.py` |
| Table S9. Early organ support and MAKE30 by trajectory class | `03_Analysis/FigureS5_organ_support/analyse_and_plot.R`, `build_outputs.py` |
| Tables S10–S11. Subgroup associations in eICU-CRD and SICdb | Same subgroup analysis/Word-export chain as Table 2 |
| Table S12. Sensitivity to missing trajectory measurements | Upstream V331 trajectory imputations; **`03_Analysis/S4_Model3_20261019/run_s4_model3.R`**; `05_Tables/update_s4_document.py` |

The table formatter names retain their historical numbering. For example, `S4_Model3_20261019` supports the current Table S12, and `FigureS5_organ_support` supports current Figure 6. This mapping records presentation renumbering; no statistical result is changed.

## Denominator checks for presentation

| Analysis | MIMIC-IV | eICU-CRD | SICdb |
| --- | ---: | ---: | ---: |
| Full trajectory cohort | 9,804 | 5,839 | 3,831 |
| MAKE30 association cohort | 9,758 | 5,709 | 3,635 |
| KM and RMST survival cohort | 9,706 | 5,804 | 3,811 |
| Complete-score incremental-prediction comparison | 9,758 | 5,209 | 3,635 |

MAKE30 exclusions from the trajectory cohort are 46, 130 and 196; survival exclusions from the trajectory cohort are 98, 35 and 20. Use the relevant denominator in each table or legend. These are aggregate study counts; no patient records accompany them.
