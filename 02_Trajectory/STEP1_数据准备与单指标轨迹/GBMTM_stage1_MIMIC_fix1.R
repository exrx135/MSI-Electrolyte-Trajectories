# ============================================================
# GBMTM Stage 1: final data preparation + QC + 2-group smoke test
# Fix1: corrected data.table time-point coverage aggregation and hardened posterior ID mapping
# Development cohort: MIMIC-IV
# External cohorts prepared identically: eICU, SICdb
# Model: group-based multivariate trajectory model (gbmt)
# Indicators: MSI, sodium, potassium, chloride, bicarbonate
# MSI time points: 1,3,5,...,23 h (12 x 2-h bins)
# Electrolyte time points: 3,9,15,21 h (4 x 6-h bins)
#
# IMPORTANT:
# 1) Missing trajectory values are NOT imputed here.
# 2) Global z-standardization is estimated ONLY from MIMIC.
# 3) eICU/SICdb use the MIMIC means/SDs; they are NOT re-standardized.
# 4) gbmt is therefore called with scaling = 0.
# 5) This script fits ONLY a 2-group quadratic smoke-test model.
#    Do not use the 2-group result as the final class solution.
# ============================================================

options(stringsAsFactors = FALSE)
set.seed(20260831)

# -----------------------------
# 0. Install/load packages
# -----------------------------
required_packages <- c("data.table", "gbmt")

installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(gbmt)
})

cat("\n============================================================\n")
cat("GBMTM Stage 1 started\n")
cat("R version:", R.version.string, "\n")
cat("gbmt version:", as.character(packageVersion("gbmt")), "\n")
cat("============================================================\n\n")

# -----------------------------
# 1. Select the three final CSV files
# -----------------------------
cat("Please select the FINAL MIMIC CSV:\n")
mimic_path <- file.choose()
cat("Please select the FINAL eICU CSV:\n")
eicu_path <- file.choose()
cat("Please select the FINAL SICdb CSV:\n")
sicdb_path <- file.choose()

# Output folder: created beside the selected MIMIC CSV
output_dir <- file.path(dirname(mimic_path), "GBMTM_stage1_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\nOutput directory:\n", output_dir, "\n\n")

# -----------------------------
# 2. Read files
# -----------------------------
cat("Reading CSV files...\n")
mimic <- fread(mimic_path, na.strings = c("", "NA", "NULL"))
eicu  <- fread(eicu_path,  na.strings = c("", "NA", "NULL"))
sicdb <- fread(sicdb_path, na.strings = c("", "NA", "NULL"))

cat("MIMIC dimensions:", nrow(mimic), "x", ncol(mimic), "\n")
cat("eICU dimensions :", nrow(eicu),  "x", ncol(eicu),  "\n")
cat("SICdb dimensions:", nrow(sicdb), "x", ncol(sicdb), "\n\n")

# -----------------------------
# 3. Exact trajectory-column validation
# -----------------------------
msi_cols <- paste0("msi_m", 1:12)
sodium_cols <- paste0("sodium_t", 1:4)
potassium_cols <- paste0("potassium_t", 1:4)
chloride_cols <- paste0("chloride_t", 1:4)
bicarbonate_cols <- paste0("bicarbonate_t", 1:4)
trajectory_cols <- c(
  msi_cols,
  sodium_cols,
  potassium_cols,
  chloride_cols,
  bicarbonate_cols
)

check_required_columns <- function(dt, db_name, id_col) {
  required <- c(id_col, trajectory_cols, "make30")
  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        db_name,
        " is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  if (anyDuplicated(dt[[id_col]]) > 0) {
    stop(paste0(db_name, " has duplicated values in ", id_col, "."))
  }
}

check_required_columns(mimic, "MIMIC", "stay_id")
check_required_columns(eicu,  "eICU",  "patientunitstayid")
check_required_columns(sicdb, "SICdb", "case_id")

if (nrow(mimic) != 9804) {
  warning("MIMIC row count is not 9804. Check whether the correct final file was selected.")
}
if (nrow(eicu) != 5867) {
  warning("eICU row count is not 5867. Check whether the correct final file was selected.")
}
if (nrow(sicdb) != 3831) {
  warning("SICdb row count is not 3831. Check whether the correct final file was selected.")
}

cat("Required columns and patient-level uniqueness: PASSED\n\n")

# -----------------------------
# 4. Wide -> long conversion
# -----------------------------
# One row per patient per 2-hour MSI midpoint.
# Electrolytes are inserted only at the corresponding 6-hour-bin midpoints:
# T1 -> 3 h, T2 -> 9 h, T3 -> 15 h, T4 -> 21 h.
# No interpolation, LOCF, NOCB, or other filling is performed.

make_long <- function(dt, id_col, db_name) {
  id_values <- as.character(dt[[id_col]])
  n <- nrow(dt)

  # 12 common time rows per patient, defined by MSI bin midpoints
  time_hour <- seq(1, 23, by = 2)

  out <- data.table(
    unit_id = rep(id_values, each = 12),
    database = db_name,
    time_hour = rep(time_hour, times = n)
  )

  # MSI: 12 time points
  msi_matrix <- as.matrix(dt[, ..msi_cols])
  storage.mode(msi_matrix) <- "double"
  out[, msi := as.vector(t(msi_matrix))]

  # Initialize electrolyte markers as missing
  out[, `:=`(
    sodium = NA_real_,
    potassium = NA_real_,
    chloride = NA_real_,
    bicarbonate = NA_real_
  )]

  electrolyte_hours <- c(3, 9, 15, 21)

  sodium_matrix <- as.matrix(dt[, ..sodium_cols])
  potassium_matrix <- as.matrix(dt[, ..potassium_cols])
  chloride_matrix <- as.matrix(dt[, ..chloride_cols])
  bicarbonate_matrix <- as.matrix(dt[, ..bicarbonate_cols])

  storage.mode(sodium_matrix) <- "double"
  storage.mode(potassium_matrix) <- "double"
  storage.mode(chloride_matrix) <- "double"
  storage.mode(bicarbonate_matrix) <- "double"

  for (j in seq_along(electrolyte_hours)) {
    h <- electrolyte_hours[j]
    idx <- which(out$time_hour == h)
    out[idx, sodium := sodium_matrix[, j]]
    out[idx, potassium := potassium_matrix[, j]]
    out[idx, chloride := chloride_matrix[, j]]
    out[idx, bicarbonate := bicarbonate_matrix[, j]]
  }

  # Center/rescale time only for numerical stability in polynomial fitting.
  # This is a linear transformation of the same 0-24 h time axis.
  # The original human-readable hour is retained in time_hour.
  out[, time_model := (time_hour - 12) / 12]

  setcolorder(
    out,
    c(
      "unit_id", "database", "time_hour", "time_model",
      "msi", "sodium", "potassium", "chloride", "bicarbonate"
    )
  )

  out
}

cat("Converting final wide files to longitudinal GBMTM format...\n")
mimic_long <- make_long(mimic, "stay_id", "MIMIC")
eicu_long  <- make_long(eicu,  "patientunitstayid", "eICU")
sicdb_long <- make_long(sicdb, "case_id", "SICdb")

# Validate exactly 12 unique time rows per patient
validate_long <- function(dt, expected_n, db_name) {
  if (nrow(dt) != expected_n * 12) {
    stop(paste0(db_name, ": unexpected number of long-format rows."))
  }

  duplicate_id_time <- dt[, .N, by = .(unit_id, time_hour)][N > 1]
  if (nrow(duplicate_id_time) > 0) {
    stop(paste0(db_name, ": duplicate unit_id + time_hour rows detected."))
  }

  n_time <- dt[, uniqueN(time_hour), by = unit_id]
  if (any(n_time$V1 != 12)) {
    stop(paste0(db_name, ": at least one patient does not have exactly 12 time rows."))
  }
}

validate_long(mimic_long, nrow(mimic), "MIMIC")
validate_long(eicu_long,  nrow(eicu),  "eICU")
validate_long(sicdb_long, nrow(sicdb), "SICdb")

cat("Long-format construction: PASSED\n\n")

# -----------------------------
# 5. Pre-model QC on ORIGINAL scale
# -----------------------------
markers <- c("msi", "sodium", "potassium", "chloride", "bicarbonate")

marker_qc <- function(dt, db_name) {
  rbindlist(lapply(markers, function(v) {
    x <- dt[[v]]
    x_obs <- x[is.finite(x)]
    data.table(
      database = db_name,
      marker = v,
      n_total_rows = length(x),
      n_observed = length(x_obs),
      missing_pct = 100 * mean(is.na(x)),
      mean = mean(x_obs),
      sd = sd(x_obs),
      p01 = as.numeric(quantile(x_obs, 0.01, names = FALSE)),
      p05 = as.numeric(quantile(x_obs, 0.05, names = FALSE)),
      median = median(x_obs),
      p95 = as.numeric(quantile(x_obs, 0.95, names = FALSE)),
      p99 = as.numeric(quantile(x_obs, 0.99, names = FALSE)),
      min = min(x_obs),
      max = max(x_obs)
    )
  }))
}

qc_original <- rbindlist(list(
  marker_qc(mimic_long, "MIMIC"),
  marker_qc(eicu_long,  "eICU"),
  marker_qc(sicdb_long, "SICdb")
))

fwrite(qc_original, file.path(output_dir, "01_original_marker_QC.csv"))

# Per-time-point missingness on original scale
coverage_qc <- function(dt, db_name) {
  # data.table requires every expression inside `by` to have either one
  # value per row or a compatible grouping length. `db_name` is a single
  # scalar string, so we first aggregate by time_hour and then add the
  # database label explicitly.
  out <- dt[,
    lapply(.SD, function(x) 100 * mean(!is.na(x))),
    by = time_hour,
    .SDcols = markers
  ]
  out[, database := db_name]
  setcolorder(out, c("database", "time_hour", markers))
  out
}

coverage_all <- rbindlist(list(
  coverage_qc(mimic_long, "MIMIC"),
  coverage_qc(eicu_long,  "eICU"),
  coverage_qc(sicdb_long, "SICdb")
))

fwrite(coverage_all, file.path(output_dir, "02_timepoint_coverage_pct.csv"))

# -----------------------------
# 6. GLOBAL standardization using MIMIC only
# -----------------------------
# Why not gbmt(scaling=2)?
# gbmt scaling=2 standardizes each indicator WITHIN EACH UNIT/PATIENT.
# For this clinical phenotype study, we want to retain between-patient level
# differences (e.g. one patient being systematically more hypernatremic than
# another), while putting the five indicators on comparable scales.
# Therefore we estimate one global mean and SD per marker from MIMIC and then
# call gbmt with scaling=0.

scaling_stats <- rbindlist(lapply(markers, function(v) {
  x <- mimic_long[[v]]
  x <- x[is.finite(x)]
  data.table(
    marker = v,
    mimic_mean = mean(x),
    mimic_sd = sd(x)
  )
}))

if (any(!is.finite(scaling_stats$mimic_sd)) || any(scaling_stats$mimic_sd <= 0)) {
  stop("Invalid MIMIC scaling SD detected.")
}

fwrite(scaling_stats, file.path(output_dir, "03_MIMIC_global_scaling_parameters.csv"))

apply_mimic_scaling <- function(dt, stats) {
  out <- copy(dt)
  for (v in markers) {
    mu <- stats[marker == v, mimic_mean]
    sig <- stats[marker == v, mimic_sd]
    new_name <- paste0(v, "_z")
    out[, (new_name) := (get(v) - mu) / sig]
  }
  out
}

mimic_long_z <- apply_mimic_scaling(mimic_long, scaling_stats)
eicu_long_z  <- apply_mimic_scaling(eicu_long,  scaling_stats)
sicdb_long_z <- apply_mimic_scaling(sicdb_long, scaling_stats)

z_markers <- paste0(markers, "_z")

# Verify development-cohort global standardization
z_qc <- rbindlist(lapply(z_markers, function(v) {
  x <- mimic_long_z[[v]]
  x <- x[is.finite(x)]
  data.table(
    marker_z = v,
    mean_after_scaling = mean(x),
    sd_after_scaling = sd(x)
  )
}))

fwrite(z_qc, file.path(output_dir, "04_MIMIC_scaled_marker_QC.csv"))

# Save prepared long datasets.
# RDS preserves NA and numeric types efficiently.
saveRDS(mimic_long_z, file.path(output_dir, "MIMIC_GBMTM_long_scaled.rds"))
saveRDS(eicu_long_z,  file.path(output_dir, "eICU_GBMTM_long_scaled.rds"))
saveRDS(sicdb_long_z, file.path(output_dir, "SICdb_GBMTM_long_scaled.rds"))

# Also save the MIMIC long dataset as CSV for easy visual inspection.
fwrite(
  mimic_long_z,
  file.path(output_dir, "MIMIC_GBMTM_long_scaled_preview.csv")
)

cat("Global MIMIC scaling completed.\n")
print(scaling_stats)
cat("\n")

# -----------------------------
# 7. 2-group quadratic SMOKE-TEST model
# -----------------------------
# This is ONLY to verify that the data structure and gbmt model run correctly.
# It is NOT the final class-number selection.
#
# nstart=3 is intentionally modest for the smoke test.
# Final 2-10 group enumeration will use more robust random restarts.

model_data <- as.data.frame(
  mimic_long_z[, c("unit_id", "time_model", z_markers), with = FALSE]
)

# gbmt requires unit to be character/factor, time numeric, indicators numeric.
model_data$unit_id <- as.character(model_data$unit_id)
model_data$time_model <- as.numeric(model_data$time_model)
for (v in z_markers) {
  model_data[[v]] <- as.numeric(model_data[[v]])
}

if (anyNA(model_data$unit_id) || anyNA(model_data$time_model)) {
  stop("unit_id or time_model contains missing values, which gbmt does not allow.")
}

if (anyDuplicated(model_data[, c("unit_id", "time_model")]) > 0) {
  stop("Duplicate unit_id + time_model rows exist before gbmt fitting.")
}

cat("============================================================\n")
cat("Starting 2-group quadratic GBMTM smoke test...\n")
cat("This may take some time. Do not close RStudio while it runs.\n")
cat("============================================================\n\n")

fit_g2_d2 <- gbmt(
  x.names = z_markers,
  unit = "unit_id",
  time = "time_model",
  ng = 2,
  d = 2,
  data = model_data,
  scaling = 0,
  pruning = TRUE,
  nstart = 3,
  tol = 1e-4,
  maxit = 1000,
  quiet = FALSE
)

# -----------------------------
# 8. Save the smoke-test model and diagnostics
# -----------------------------
saveRDS(fit_g2_d2, file.path(output_dir, "fit_G2_D2_smoketest.rds"))

# Basic fit statistics
fit_stats <- data.table(
  model = "G2_D2_smoketest",
  groups = 2,
  degree = 2,
  logLik = as.numeric(logLik(fit_g2_d2)),
  AIC = as.numeric(AIC(fit_g2_d2)),
  BIC = as.numeric(BIC(fit_g2_d2)),
  npar = fit_g2_d2$npar
)
fwrite(fit_stats, file.path(output_dir, "05_G2_D2_fit_statistics.csv"))

# Raw information criteria object from gbmt
capture.output(
  print(fit_g2_d2$ic),
  file = file.path(output_dir, "06_G2_D2_gbmt_ic.txt")
)

# EM runs and convergence
em_table <- as.data.table(fit_g2_d2$em, keep.rownames = "run")
fwrite(em_table, file.path(output_dir, "07_G2_D2_EM_runs.csv"))

# Posterior probabilities, class assignment, APPA and group size
post_mat <- as.matrix(fit_g2_d2$posterior)
assign_vec <- as.integer(fit_g2_d2$assign)

# gbmt stores the unit IDs as row names of the posterior matrix. Using those
# row names is safer than reconstructing an order from the original data.
unit_order <- rownames(post_mat)
if (is.null(unit_order) || length(unit_order) != nrow(post_mat)) {
  unit_order <- unique(as.character(fit_g2_d2$data.orig$unit_id))
}

if (nrow(post_mat) != length(unit_order)) {
  stop("Posterior row count does not match the number of model units.")
}
if (length(assign_vec) != length(unit_order)) {
  stop("Assignment count does not match the number of model units.")
}
if (ncol(post_mat) != 2) {
  stop("Unexpected posterior matrix column count for a 2-group model.")
}

posterior_dt <- as.data.table(post_mat)
setnames(posterior_dt, paste0("posterior_g", seq_len(ncol(post_mat))))
posterior_dt[, unit_id := unit_order]
posterior_dt[, assigned_group := assign_vec]
setcolorder(
  posterior_dt,
  c("unit_id", "assigned_group", paste0("posterior_g", 1:2))
)

fwrite(posterior_dt, file.path(output_dir, "08_G2_D2_patient_posterior.csv"))

group_summary <- posterior_dt[, {
  g <- unique(assigned_group)
  selected_post <- get(paste0("posterior_g", g))
  .(
    n = .N,
    proportion = .N / nrow(posterior_dt),
    APPA_manual = mean(selected_post, na.rm = TRUE),
    median_assigned_posterior = median(selected_post, na.rm = TRUE),
    p10_assigned_posterior = as.numeric(quantile(selected_post, 0.10, na.rm = TRUE))
  )
}, by = assigned_group][order(assigned_group)]

fwrite(group_summary, file.path(output_dir, "09_G2_D2_group_summary.csv"))

# Save package-provided APPA separately for comparison
capture.output(
  print(fit_g2_d2$appa),
  file = file.path(output_dir, "10_G2_D2_package_APPA.txt")
)

# Model regression summaries
capture.output(
  summary(fit_g2_d2),
  file = file.path(output_dir, "11_G2_D2_model_summary.txt")
)

# Built-in trajectory plot on standardized scale
pdf(
  file.path(output_dir, "12_G2_D2_trajectory_plot_standardized.pdf"),
  width = 12,
  height = 9
)
plot(
  fit_g2_d2,
  bands = TRUE,
  observed = FALSE,
  equal.scale = TRUE,
  add.legend = TRUE,
  mar = c(4.0, 4.0, 3.0, 1.5)
)
dev.off()

# -----------------------------
# 9. Console summary
# -----------------------------
cat("\n============================================================\n")
cat("STAGE 1 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")
cat("Development cohort (MIMIC):", nrow(mimic), "patients\n")
cat("MIMIC long rows:", nrow(mimic_long_z), "\n")
cat("External eICU prepared:", nrow(eicu), "patients\n")
cat("External SICdb prepared:", nrow(sicdb), "patients\n\n")
cat("Smoke-test fit statistics:\n")
print(fit_stats)
cat("\nSmoke-test group summary:\n")
print(group_summary)
cat("\nEM run diagnostics:\n")
print(em_table)
cat("\nAll outputs were written to:\n", output_dir, "\n")
cat("\nDO NOT choose the final number of trajectory groups from this 2-group smoke test.\n")
cat("The next stage is systematic 2-10 group model enumeration after this test passes.\n")
cat("============================================================\n")
