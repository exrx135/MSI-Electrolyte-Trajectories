#!/usr/bin/env Rscript

# =============================================================================
# V294_stage3_G2G20.R
# FOURTH STUDY / MSI-Electrolyte GBMTM
#
# PURPOSE
#   Local revised Stage 3 structural adjudication for G2-G20.
#
# INPUT
#   Exact compact outputs created by V292_step2_visual_G2G20.R:
#
#   <current Windows profile>/Desktop/Research/MSI/提取轨迹文件/
#   GBMTM_stage2_output/V292_step2_visual_G2G20/tables
#
# SCOPE
#   - NO GBMTM refitting.
#   - NO MAKE30, mortality, KRT, or any outcome data.
#   - Compare G2 -> G3 -> ... -> G20 continuously.
#   - Separate persistent absolute LEVEL differences from temporal SHAPE.
#   - Retain full fitted-point shape for structural comparison.
#   - Restrict electrolyte clinical shape evidence to 3/9/15/21 h anchors.
#   - Use adjacent-G lineage, ARI/NMI, posterior quality, small-class burden,
#     level-vs-shape decomposition, and multivariate coupling.
#   - No single fixed RMSE cutoff is allowed to define a new phenotype.
#
# OUTPUT
#   <current Windows profile>/Desktop/Research/MSI/提取轨迹文件/
#   GBMTM_stage3_output/
#
#   The script intentionally DOES NOT automatically freeze the final G or the
#   Step-4A shortlist. It creates the complete evidence package for final
#   statistical/clinical adjudication after upload.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(zip)
})

VERSION <- "V294"
SOURCE_VERSION <- "V292"

G_VALUES <- 2:20
LOW_G_FOCUS <- 2:6

INDICATORS <- c("MSI", "Na", "K", "Cl", "HCO3")
ELECTROLYTES <- c("Na", "K", "Cl", "HCO3")

MSI_GRID <- seq(1, 23, by = 2)
ELECTROLYTE_ANCHORS <- c(3, 9, 15, 21)
COMMON_COUPLING_ANCHORS <- c(3, 9, 15, 21)

# These proportions are lineage-descriptive signals only.
# They are NOT trajectory-distance thresholds.
MINOR_SPLIT_PCT <- 5
SUBSTANTIVE_SPLIT_PCT <- 10

# =============================================================================
# 1. Paths
# =============================================================================

home_r <- path.expand("~")
home_env <- Sys.getenv("USERPROFILE", unset = "")

user_home <- if (nzchar(home_env) && dir.exists(home_env)) {
  home_env
} else {
  home_r
}

if (!dir.exists(user_home)) {
  stop(
    paste0(
      "Could not resolve current Windows profile root.\n",
      "path.expand(\"~\") = ", home_r, "\n",
      "USERPROFILE = ", ifelse(nzchar(home_env), home_env, "<empty>")
    ),
    call. = FALSE
  )
}

user_home <- normalizePath(
  user_home,
  winslash = "/",
  mustWork = TRUE
)

project_dir <- file.path(
  user_home,
  "Desktop", "Research", "MSI", "提取轨迹文件"
)

stage2_base <- file.path(
  project_dir,
  "GBMTM_stage2_output"
)

stage3_base <- file.path(
  project_dir,
  "GBMTM_stage3_output"
)

source_root <- file.path(
  stage2_base,
  "V292_step2_visual_G2G20"
)

source_table_dir <- file.path(
  source_root,
  "tables"
)

source_figure_dir <- file.path(
  source_root,
  "figures"
)

# The user explicitly requested Stage-3 results under GBMTM_stage3_output.
# No Windows username is hard-coded; only the frozen relative hierarchy is used.
output_dir <- stage3_base

table_dir <- file.path(
  output_dir,
  "tables"
)

figure_dir <- file.path(
  output_dir,
  "source_figures"
)

resume_dir <- file.path(
  output_dir,
  "_resume"
)

required_dirs <- c(
  project_dir,
  stage2_base,
  source_root,
  source_table_dir
)

missing_dirs <- required_dirs[
  !dir.exists(required_dirs)
]

if (length(missing_dirs) > 0L) {
  stop(
    paste0(
      "Required directory does not exist. Checked exact path(s):\n",
      paste(
        normalizePath(
          missing_dirs,
          winslash = "/",
          mustWork = FALSE
        ),
        collapse = "\n"
      ),
      "\nRun V292_step2_visual_G2G20.R successfully before V294."
    ),
    call. = FALSE
  )
}

for (d in c(
  output_dir,
  table_dir,
  figure_dir,
  resume_dir
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(d)) {
    stop(
      paste0(
        "Could not create output directory:\n",
        normalizePath(
          d,
          winslash = "/",
          mustWork = FALSE
        )
      ),
      call. = FALSE
    )
  }
}

log_path <- file.path(
  output_dir,
  "V294_stage3_log_G2G20.txt"
)

report_path <- file.path(
  output_dir,
  "V294_stage3_report_G2G20.txt"
)

zip_path <- file.path(
  output_dir,
  "V294_stage3_bundle_G2G20.zip"
)

# =============================================================================
# 2. Exact required V292 inputs
# =============================================================================

required_input_names <- c(
  "V292_trajectory_all_G2G20.csv",
  "V292_class_quality_all_G2G20.csv",
  "V292_model_summary_all_G2G20.csv",
  "V292_assignment_all_G2G20.csv",
  "V292_trajectory_summary_G2G20.csv",
  "V292_level_shape_by_indicator_G2G20.csv",
  "V292_level_shape_joint_G2G20.csv",
  "V292_nearest_shape_pair_G2G20.csv",
  "V292_transition_G2G20.csv",
  "V292_transition_metrics_G2G20.csv"
)

required_input_paths <- file.path(
  source_table_dir,
  required_input_names
)

missing_inputs <- required_input_paths[
  !file.exists(required_input_paths)
]

if (length(missing_inputs) > 0L) {
  stop(
    paste0(
      "Required V292 compact input file(s) are missing.\n",
      "Checked exact path(s):\n",
      paste(
        normalizePath(
          missing_inputs,
          winslash = "/",
          mustWork = FALSE
        ),
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}

# =============================================================================
# 3. Helpers
# =============================================================================

log_message <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0(
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    " | ",
    msg
  )

  cat(line, "\n")
  cat(
    line,
    "\n",
    file = log_path,
    append = TRUE,
    sep = ""
  )

  flush.console()
}

required_columns <- function(x, cols, label) {
  missing_cols <- setdiff(
    cols,
    names(x)
  )

  if (length(missing_cols) > 0L) {
    stop(
      paste0(
        label,
        " is missing exact required field(s): ",
        paste(
          missing_cols,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::median(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 3L) {
    return(NA_real_)
  }

  xx <- x[ok]
  yy <- y[ok]

  sx <- safe_sd(xx)
  sy <- safe_sd(yy)

  if (
    !is.finite(sx) ||
    !is.finite(sy) ||
    sx == 0 ||
    sy == 0
  ) {
    return(NA_real_)
  }

  stats::cor(xx, yy)
}

write_csv_exact <- function(x, path) {
  fwrite(
    x,
    path,
    na = ""
  )

  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size <= 0
  ) {
    stop(
      "Failed to write output file: ",
      path,
      call. = FALSE
    )
  }
}

md5_one <- function(path) {
  value <- tolower(
    unname(
      tools::md5sum(path)
    )
  )

  if (
    length(value) != 1L ||
    is.na(value) ||
    !nzchar(value)
  ) {
    stop(
      "Failed to calculate MD5: ",
      path,
      call. = FALSE
    )
  }

  value
}

rank_fraction <- function(x) {
  out <- rep(
    NA_real_,
    length(x)
  )

  ok <- is.finite(x)
  n_ok <- sum(ok)

  if (n_ok == 0L) {
    return(out)
  }

  if (n_ok == 1L) {
    out[ok] <- 0.5
    return(out)
  }

  rr <- frank(
    x[ok],
    ties.method = "average"
  )

  out[ok] <- (
    rr - 1
  ) / (
    n_ok - 1
  )

  out
}

normalised_entropy <- function(p) {
  p <- p[
    is.finite(p) &
      p > 0
  ]

  if (length(p) <= 1L) {
    return(0)
  }

  p <- p / sum(p)

  h <- -sum(
    p * log(p)
  )

  h / log(length(p))
}

turn_count <- function(x) {
  if (length(x) < 3L) {
    return(0L)
  }

  d <- diff(x)
  d <- d[
    is.finite(d) &
      d != 0
  ]

  if (length(d) < 2L) {
    return(0L)
  }

  as.integer(
    sum(
      sign(d[-1L]) !=
        sign(d[-length(d)])
    )
  )
}

curvature_rms <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 3L) {
    return(NA_real_)
  }

  d2 <- diff(
    x,
    differences = 2L
  )

  sqrt(
    mean(d2^2)
  )
}

stage_marker_path <- function(stage_name) {
  file.path(
    resume_dir,
    paste0(
      "V294_",
      stage_name,
      "_done.rds"
    )
  )
}

stage_is_valid <- function(
  stage_name,
  source_fingerprint,
  output_paths
) {
  marker_path <- stage_marker_path(
    stage_name
  )

  if (!file.exists(marker_path)) {
    return(FALSE)
  }

  if (
    length(output_paths) == 0L ||
    any(!file.exists(output_paths))
  ) {
    return(FALSE)
  }

  marker <- tryCatch(
    readRDS(marker_path),
    error = function(e) NULL
  )

  if (is.null(marker)) {
    return(FALSE)
  }

  if (!identical(
    marker$version,
    VERSION
  )) {
    return(FALSE)
  }

  if (!identical(
    marker$source_fingerprint,
    source_fingerprint
  )) {
    return(FALSE)
  }

  current_paths <- normalizePath(
    output_paths,
    winslash = "/",
    mustWork = TRUE
  )

  if (!identical(
    marker$output_paths,
    current_paths
  )) {
    return(FALSE)
  }

  current_md5 <- vapply(
    output_paths,
    md5_one,
    character(1)
  )

  identical(
    unname(marker$output_md5),
    unname(current_md5)
  )
}

save_stage_marker <- function(
  stage_name,
  source_fingerprint,
  output_paths
) {
  marker <- list(
    version = VERSION,
    stage = stage_name,
    completed_at = format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    source_fingerprint =
      source_fingerprint,
    output_paths = normalizePath(
      output_paths,
      winslash = "/",
      mustWork = TRUE
    ),
    output_md5 = vapply(
      output_paths,
      md5_one,
      character(1)
    )
  )

  saveRDS(
    marker,
    stage_marker_path(stage_name)
  )
}

get_pair_row <- function(
  pair_dt,
  g,
  class_a,
  class_b
) {
  lo <- min(
    class_a,
    class_b
  )

  hi <- max(
    class_a,
    class_b
  )

  out <- pair_dt[
    requested_g == g &
      class_a == lo &
      class_b == hi
  ]

  if (nrow(out) == 0L) {
    return(NULL)
  }

  out[1L]
}

fmt <- function(x, digits = 4L) {
  ifelse(
    is.na(x),
    "NA",
    format(
      round(x, digits),
      trim = TRUE,
      scientific = FALSE
    )
  )
}

# =============================================================================
# 4. Source fingerprint
# =============================================================================

source_manifest <- data.table(
  file_name = required_input_names,
  full_path = normalizePath(
    required_input_paths,
    winslash = "/",
    mustWork = TRUE
  ),
  size_bytes = as.numeric(
    file.info(required_input_paths)$size
  ),
  modified_time = format(
    file.info(required_input_paths)$mtime,
    "%Y-%m-%d %H:%M:%S"
  ),
  md5 = vapply(
    required_input_paths,
    md5_one,
    character(1)
  )
)

source_fingerprint <- paste(
  paste0(
    source_manifest$file_name,
    "=",
    source_manifest$md5
  ),
  collapse = "|"
)

source_manifest_path <- file.path(
  table_dir,
  "V294_source_manifest_G2G20.csv"
)

write_csv_exact(
  source_manifest,
  source_manifest_path
)

# =============================================================================
# 5. Read and validate V292 compact inputs
# =============================================================================

log_message(
  VERSION,
  " started"
)

log_message(
  "V292 source directory: ",
  normalizePath(
    source_root,
    winslash = "/",
    mustWork = TRUE
  )
)

trajectory <- fread(
  file.path(
    source_table_dir,
    "V292_trajectory_all_G2G20.csv"
  )
)

class_quality <- fread(
  file.path(
    source_table_dir,
    "V292_class_quality_all_G2G20.csv"
  )
)

model_summary <- fread(
  file.path(
    source_table_dir,
    "V292_model_summary_all_G2G20.csv"
  )
)

assignment <- fread(
  file.path(
    source_table_dir,
    "V292_assignment_all_G2G20.csv"
  )
)

trajectory_summary_v292 <- fread(
  file.path(
    source_table_dir,
    "V292_trajectory_summary_G2G20.csv"
  )
)

level_shape_indicator <- fread(
  file.path(
    source_table_dir,
    "V292_level_shape_by_indicator_G2G20.csv"
  )
)

level_shape_joint <- fread(
  file.path(
    source_table_dir,
    "V292_level_shape_joint_G2G20.csv"
  )
)

nearest_shape_v292 <- fread(
  file.path(
    source_table_dir,
    "V292_nearest_shape_pair_G2G20.csv"
  )
)

transition <- fread(
  file.path(
    source_table_dir,
    "V292_transition_G2G20.csv"
  )
)

transition_metrics <- fread(
  file.path(
    source_table_dir,
    "V292_transition_metrics_G2G20.csv"
  )
)

required_columns(
  trajectory,
  c(
    "requested_g",
    "class",
    "class_label",
    "n",
    "class_pct",
    "APPA_model",
    "marker",
    "indicator",
    "model_time",
    "time_hour_center",
    "fitted_z",
    "fitted_native",
    "centered_z"
  ),
  "V292_trajectory_all_G2G20.csv"
)

required_columns(
  class_quality,
  c(
    "requested_g",
    "class",
    "n",
    "class_pct",
    "APPA_model",
    "posterior_soft_mass",
    "posterior_soft_mass_pct",
    "mean_max_posterior_in_assigned_class",
    "median_max_posterior_in_assigned_class",
    "class_label"
  ),
  "V292_class_quality_all_G2G20.csv"
)

required_columns(
  model_summary,
  c(
    "requested_g",
    "n_patients",
    "status",
    "final_logLik",
    "npar",
    "BIC",
    "min_class_pct",
    "min_APPA",
    "posterior_rows",
    "posterior_columns",
    "assignment_modal_mismatch_n",
    "mean_max_posterior",
    "median_max_posterior",
    "mean_posterior_entropy_normalised",
    "em_iterations",
    "package_converged",
    "source_model_file",
    "extraction_seconds"
  ),
  "V292_model_summary_all_G2G20.csv"
)

required_columns(
  assignment,
  c(
    "requested_g",
    "unit_id",
    "assigned_class",
    "modal_class",
    "assignment_modal_mismatch",
    "max_posterior",
    "posterior_entropy_normalised"
  ),
  "V292_assignment_all_G2G20.csv"
)

required_columns(
  level_shape_indicator,
  c(
    "requested_g",
    "class_a",
    "class_b",
    "indicator",
    "mean_level_offset_z",
    "absolute_level_offset_z",
    "mean_level_offset_native",
    "absolute_level_offset_native",
    "total_RMS_difference_z",
    "baseline_centered_shape_RMS_z",
    "baseline_centered_shape_max_abs_z",
    "baseline_centered_shape_correlation"
  ),
  "V292_level_shape_by_indicator_G2G20.csv"
)

required_columns(
  level_shape_joint,
  c(
    "requested_g",
    "class_a",
    "class_b",
    "indicator_count",
    "joint_level_RMS_z",
    "joint_total_RMS_z",
    "joint_centered_shape_RMS_z",
    "mean_centered_shape_correlation",
    "shape_share_of_level_plus_shape"
  ),
  "V292_level_shape_joint_G2G20.csv"
)

required_columns(
  transition,
  c(
    "from_g",
    "to_g",
    "class_from",
    "class_to",
    "N",
    "row_total",
    "row_pct",
    "column_total",
    "column_pct"
  ),
  "V292_transition_G2G20.csv"
)

required_columns(
  transition_metrics,
  c(
    "from_g",
    "to_g",
    "n_units",
    "adjusted_rand_index",
    "normalised_mutual_information"
  ),
  "V292_transition_metrics_G2G20.csv"
)

# Explicit numeric conversion of exact fields.
for (col in c(
  "requested_g",
  "class",
  "n",
  "class_pct",
  "APPA_model",
  "model_time",
  "time_hour_center",
  "fitted_z",
  "fitted_native",
  "centered_z"
)) {
  set(
    trajectory,
    j = col,
    value = safe_numeric(
      trajectory[[col]]
    )
  )
}

for (col in c(
  "requested_g",
  "class",
  "n",
  "class_pct",
  "APPA_model",
  "posterior_soft_mass",
  "posterior_soft_mass_pct",
  "mean_max_posterior_in_assigned_class",
  "median_max_posterior_in_assigned_class"
)) {
  set(
    class_quality,
    j = col,
    value = safe_numeric(
      class_quality[[col]]
    )
  )
}

for (col in c(
  "requested_g",
  "n_patients",
  "final_logLik",
  "npar",
  "BIC",
  "min_class_pct",
  "min_APPA",
  "posterior_rows",
  "posterior_columns",
  "assignment_modal_mismatch_n",
  "mean_max_posterior",
  "median_max_posterior",
  "mean_posterior_entropy_normalised",
  "em_iterations",
  "package_converged",
  "extraction_seconds"
)) {
  set(
    model_summary,
    j = col,
    value = safe_numeric(
      model_summary[[col]]
    )
  )
}

for (col in c(
  "requested_g",
  "assigned_class",
  "modal_class",
  "assignment_modal_mismatch",
  "max_posterior",
  "posterior_entropy_normalised"
)) {
  set(
    assignment,
    j = col,
    value = safe_numeric(
      assignment[[col]]
    )
  )
}

assignment[
  ,
  unit_id := as.character(unit_id)
]

trajectory[
  ,
  indicator := as.character(indicator)
]

observed_g <- sort(
  unique(
    as.integer(
      model_summary$requested_g
    )
  )
)

if (!identical(
  observed_g,
  G_VALUES
)) {
  stop(
    paste0(
      "V292 model summary does not contain exactly G2-G20.\n",
      "Expected: ",
      paste(G_VALUES, collapse = ","),
      "\nObserved: ",
      paste(observed_g, collapse = ",")
    ),
    call. = FALSE
  )
}

for (g in G_VALUES) {
  q <- class_quality[
    requested_g == g
  ]

  if (
    nrow(q) != g ||
    uniqueN(q$class) != g
  ) {
    stop(
      "Class-quality structure invalid for G",
      g,
      ".",
      call. = FALSE
    )
  }

  a <- assignment[
    requested_g == g
  ]

  if (
    nrow(a) == 0L ||
    anyDuplicated(a$unit_id) > 0L
  ) {
    stop(
      "Assignment structure invalid for G",
      g,
      ".",
      call. = FALSE
    )
  }

  if (
    uniqueN(a$assigned_class) != g
  ) {
    stop(
      "G",
      g,
      " does not contain exactly ",
      g,
      " occupied classes.",
      call. = FALSE
    )
  }

  t <- trajectory[
    requested_g == g
  ]

  if (
    uniqueN(t$class) != g
  ) {
    stop(
      "Trajectory class count invalid for G",
      g,
      ".",
      call. = FALSE
    )
  }

  if (
    !setequal(
      unique(t$indicator),
      INDICATORS
    )
  ) {
    stop(
      paste0(
        "Trajectory indicator set invalid for G",
        g,
        ".\nExpected: ",
        paste(INDICATORS, collapse = ","),
        "\nObserved: ",
        paste(
          sort(unique(t$indicator)),
          collapse = ","
        )
      ),
      call. = FALSE
    )
  }
}

reference_units <- sort(
  assignment[
    requested_g == G_VALUES[1L],
    unit_id
  ]
)

for (g in G_VALUES[-1L]) {
  current_units <- sort(
    assignment[
      requested_g == g,
      unit_id
    ]
  )

  if (!identical(
    reference_units,
    current_units
  )) {
    stop(
      "Patient unit_id set differs between G2 and G",
      g,
      ".",
      call. = FALSE
    )
  }
}

# =============================================================================
# 6. Stage A: class-level clinical structure
# =============================================================================

anchor_path <- file.path(
  table_dir,
  "V294_class_anchor_G2G20.csv"
)

class_structure_path <- file.path(
  table_dir,
  "V294_class_structure_G2G20.csv"
)

class_coupling_path <- file.path(
  table_dir,
  "V294_class_coupling_G2G20.csv"
)

stage_a_outputs <- c(
  anchor_path,
  class_structure_path,
  class_coupling_path
)

if (
  stage_is_valid(
    "stageA_class_structure",
    source_fingerprint,
    stage_a_outputs
  )
) {
  log_message(
    "Stage A valid checkpoint found; reusing outputs"
  )

  class_anchor <- fread(anchor_path)
  class_structure <- fread(
    class_structure_path
  )
  class_coupling <- fread(
    class_coupling_path
  )
} else {
  log_message(
    "Stage A computing class-level clinical structure"
  )

  anchor_parts <- list()
  structure_parts <- list()
  coupling_parts <- list()

  anchor_idx <- 0L
  structure_idx <- 0L
  coupling_idx <- 0L

  for (g in G_VALUES) {
    for (cl in seq_len(g)) {
      for (indicator_now in INDICATORS) {
        d <- trajectory[
          requested_g == g &
            class == cl &
            indicator == indicator_now
        ][order(time_hour_center)]

        if (nrow(d) == 0L) {
          stop(
            "No trajectory rows for G",
            g,
            " C",
            cl,
            " ",
            indicator_now,
            ".",
            call. = FALSE
          )
        }

        clinical_times <- if (
          indicator_now == "MSI"
        ) {
          MSI_GRID
        } else {
          ELECTROLYTE_ANCHORS
        }

        d_clinical <- d[
          time_hour_center %in%
            clinical_times
        ]

        setorder(
          d_clinical,
          time_hour_center
        )

        if (
          nrow(d_clinical) !=
            length(clinical_times) ||
          !isTRUE(
            all.equal(
              d_clinical$time_hour_center,
              clinical_times,
              tolerance = 1e-12,
              check.attributes = FALSE
            )
          )
        ) {
          stop(
            paste0(
              "Clinical anchor validation failed for G",
              g,
              " C",
              cl,
              " ",
              indicator_now,
              ". Expected exact times: ",
              paste(
                clinical_times,
                collapse = ","
              )
            ),
            call. = FALSE
          )
        }

        baseline_z <-
          d_clinical$fitted_z[1L]

        baseline_native <-
          d_clinical$fitted_native[1L]

        d_clinical[
          ,
          `:=`(
            clinical_centered_z =
              fitted_z - baseline_z,
            clinical_centered_native =
              fitted_native -
              baseline_native,
            clinical_anchor_type =
              if (
                indicator_now == "MSI"
              ) {
                "MSI_complete_2h_grid"
              } else {
                "electrolyte_real_3_9_15_21h"
              }
          )
        ]

        anchor_idx <- anchor_idx + 1L

        anchor_parts[[anchor_idx]] <-
          d_clinical[
            ,
            .(
              requested_g,
              class,
              class_label,
              n,
              class_pct,
              APPA_model,
              indicator,
              time_hour_center,
              fitted_z,
              fitted_native,
              clinical_centered_z,
              clinical_centered_native,
              clinical_anchor_type
            )
          ]

        z <- d_clinical$fitted_z
        native <- d_clinical$fitted_native
        centered <- z - z[1L]
        dz <- diff(z)

        structure_idx <-
          structure_idx + 1L

        structure_parts[[structure_idx]] <- data.table(
          requested_g = g,
          class = cl,
          class_label =
            as.character(
              d_clinical$class_label[1L]
            ),
          n =
            as.numeric(
              d_clinical$n[1L]
            ),
          class_pct =
            as.numeric(
              d_clinical$class_pct[1L]
            ),
          APPA_model =
            as.numeric(
              d_clinical$APPA_model[1L]
            ),
          indicator = indicator_now,
          clinical_time_definition =
            if (
              indicator_now == "MSI"
            ) {
              "1,3,5,...,23 h"
            } else {
              "3,9,15,21 h only"
            },
          start_time_h =
            d_clinical$time_hour_center[1L],
          end_time_h =
            d_clinical$time_hour_center[
              nrow(d_clinical)
            ],
          start_native = native[1L],
          end_native = native[length(native)],
          delta_native =
            native[length(native)] -
            native[1L],
          mean_native = mean(native),
          min_native = min(native),
          max_native = max(native),
          start_z = z[1L],
          end_z = z[length(z)],
          delta_z =
            z[length(z)] -
            z[1L],
          mean_z = mean(z),
          clinical_centered_shape_RMS_z =
            sqrt(
              mean(centered^2)
            ),
          clinical_centered_shape_range_z =
            max(centered) -
            min(centered),
          total_variation_z =
            if (length(dz) > 0L) {
              sum(abs(dz))
            } else {
              0
            },
          absolute_net_change_z =
            abs(
              z[length(z)] -
                z[1L]
            ),
          path_excess_z =
            if (length(dz) > 0L) {
              sum(abs(dz)) -
                abs(
                  z[length(z)] -
                    z[1L]
                )
            } else {
              0
            },
          turn_count =
            turn_count(z),
          curvature_RMS_z =
            curvature_rms(z)
        )
      }

      # Descriptive multivariate coupling at common real anchor times.
      common <- trajectory[
        requested_g == g &
          class == cl &
          indicator %in% INDICATORS &
          time_hour_center %in%
            COMMON_COUPLING_ANCHORS,
        .(
          indicator,
          time_hour_center,
          fitted_z,
          fitted_native
        )
      ]

      if (
        nrow(common) !=
          length(INDICATORS) *
          length(COMMON_COUPLING_ANCHORS)
      ) {
        stop(
          "Common-anchor coupling grid incomplete for G",
          g,
          " C",
          cl,
          ".",
          call. = FALSE
        )
      }

      wide_z <- dcast(
        common,
        time_hour_center ~ indicator,
        value.var = "fitted_z"
      )

      setorder(
        wide_z,
        time_hour_center
      )

      for (indicator_now in INDICATORS) {
        if (
          !indicator_now %in%
            names(wide_z)
        ) {
          stop(
            "Missing coupling indicator ",
            indicator_now,
            " for G",
            g,
            " C",
            cl,
            ".",
            call. = FALSE
          )
        }
      }

      coupling_row <- data.table(
        requested_g = g,
        class = cl
      )

      for (indicator_now in INDICATORS) {
        values <- wide_z[[indicator_now]]

        coupling_row[
          ,
          (
            paste0(
              "delta_z_",
              indicator_now,
              "_3to21h"
            )
          ) :=
            values[length(values)] -
            values[1L]
        ]
      }

      indicator_pairs <- combn(
        INDICATORS,
        2L
      )

      for (
        j in seq_len(
          ncol(indicator_pairs)
        )
      ) {
        a_name <-
          indicator_pairs[1L, j]

        b_name <-
          indicator_pairs[2L, j]

        coupling_row[
          ,
          (
            paste0(
              "corr_",
              a_name,
              "_",
              b_name,
              "_3_9_15_21h"
            )
          ) :=
            safe_cor(
              wide_z[[a_name]],
              wide_z[[b_name]]
            )
        ]
      }

      coupling_idx <-
        coupling_idx + 1L

      coupling_parts[[coupling_idx]] <- coupling_row
    }
  }

  class_anchor <- rbindlist(
    anchor_parts,
    use.names = TRUE,
    fill = TRUE
  )

  class_structure <- rbindlist(
    structure_parts,
    use.names = TRUE,
    fill = TRUE
  )

  class_coupling <- rbindlist(
    coupling_parts,
    use.names = TRUE,
    fill = TRUE
  )

  class_anchor[
    ,
    indicator_order :=
      match(
        indicator,
        INDICATORS
      )
  ]

  setorder(
    class_anchor,
    requested_g,
    class,
    indicator_order,
    time_hour_center
  )

  class_anchor[
    ,
    indicator_order := NULL
  ]

  class_structure[
    ,
    indicator_order :=
      match(
        indicator,
        INDICATORS
      )
  ]

  setorder(
    class_structure,
    requested_g,
    class,
    indicator_order
  )

  class_structure[
    ,
    indicator_order := NULL
  ]

  setorder(
    class_coupling,
    requested_g,
    class
  )

  write_csv_exact(
    class_anchor,
    anchor_path
  )

  write_csv_exact(
    class_structure,
    class_structure_path
  )

  write_csv_exact(
    class_coupling,
    class_coupling_path
  )

  save_stage_marker(
    "stageA_class_structure",
    source_fingerprint,
    stage_a_outputs
  )
}

# =============================================================================
# 7. Stage B: clinical-anchor pairwise shape + pair context
# =============================================================================

clinical_pair_path <- file.path(
  table_dir,
  "V294_pair_clinical_shape_G2G20.csv"
)

pair_context_path <- file.path(
  table_dir,
  "V294_pair_context_G2G20.csv"
)

nearest_pair_path <- file.path(
  table_dir,
  "V294_nearest_pair_G2G20.csv"
)

class_nearest_path <- file.path(
  table_dir,
  "V294_class_nearest_G2G20.csv"
)

stage_b_outputs <- c(
  clinical_pair_path,
  pair_context_path,
  nearest_pair_path,
  class_nearest_path
)

if (
  stage_is_valid(
    "stageB_pair_context",
    source_fingerprint,
    stage_b_outputs
  )
) {
  log_message(
    "Stage B valid checkpoint found; reusing outputs"
  )

  clinical_pair <- fread(
    clinical_pair_path
  )

  pair_context <- fread(
    pair_context_path
  )

  nearest_pair <- fread(
    nearest_pair_path
  )

  class_nearest <- fread(
    class_nearest_path
  )
} else {
  log_message(
    "Stage B computing clinical-anchor pairwise shape context"
  )

  clinical_pair_parts <- list()
  cp_idx <- 0L

  for (g in G_VALUES) {
    classes <- seq_len(g)

    pair_matrix <- combn(
      classes,
      2L
    )

    for (
      j in seq_len(
        ncol(pair_matrix)
      )
    ) {
      class_a <-
        pair_matrix[1L, j]

      class_b <-
        pair_matrix[2L, j]

      indicator_rows <- list()
      ii <- 0L

      for (indicator_now in INDICATORS) {
        a <- class_anchor[
          requested_g == g &
            class == class_a &
            indicator == indicator_now
        ][order(time_hour_center)]

        b <- class_anchor[
          requested_g == g &
            class == class_b &
            indicator == indicator_now
        ][order(time_hour_center)]

        if (
          nrow(a) == 0L ||
          nrow(a) != nrow(b) ||
          !identical(
            a$time_hour_center,
            b$time_hour_center
          )
        ) {
          stop(
            paste0(
              "Clinical pair alignment failed for G",
              g,
              " C",
              class_a,
              " vs C",
              class_b,
              " ",
              indicator_now,
              "."
            ),
            call. = FALSE
          )
        }

        diff_z <-
          a$fitted_z -
          b$fitted_z

        centered_a <-
          a$fitted_z -
          a$fitted_z[1L]

        centered_b <-
          b$fitted_z -
          b$fitted_z[1L]

        centered_diff <-
          centered_a -
          centered_b

        ii <- ii + 1L

        indicator_rows[[ii]] <-
          data.table(
            requested_g = g,
            class_a = class_a,
            class_b = class_b,
            indicator = indicator_now,
            clinical_point_count =
              nrow(a),
            clinical_level_offset_z =
              mean(diff_z),
            clinical_absolute_level_offset_z =
              abs(mean(diff_z)),
            clinical_total_RMS_z =
              sqrt(
                mean(diff_z^2)
              ),
            clinical_centered_shape_RMS_z =
              sqrt(
                mean(centered_diff^2)
              ),
            clinical_centered_shape_max_abs_z =
              max(abs(centered_diff)),
            clinical_centered_shape_correlation =
              safe_cor(
                centered_a,
                centered_b
              )
          )
      }

      indicator_dt <- rbindlist(
        indicator_rows,
        use.names = TRUE,
        fill = TRUE
      )

      cp_idx <- cp_idx + 1L

      clinical_pair_parts[[cp_idx]] <-
        indicator_dt[
          ,
          .(
            clinical_indicator_count = .N,
            clinical_joint_level_RMS_z =
              sqrt(
                mean(
                  clinical_level_offset_z^2
                )
              ),
            clinical_joint_total_RMS_z =
              sqrt(
                mean(
                  clinical_total_RMS_z^2
                )
              ),
            clinical_joint_centered_shape_RMS_z =
              sqrt(
                mean(
                  clinical_centered_shape_RMS_z^2
                )
              ),
            clinical_mean_centered_shape_correlation =
              safe_mean(
                clinical_centered_shape_correlation
              ),
            clinical_shape_share_of_level_plus_shape = {
              shape_value <- sqrt(
                mean(
                  clinical_centered_shape_RMS_z^2
                )
              )

              level_value <- sqrt(
                mean(
                  clinical_level_offset_z^2
                )
              )

              if (
                shape_value +
                  level_value >
                  0
              ) {
                shape_value /
                  (
                    shape_value +
                      level_value
                  )
              } else {
                NA_real_
              }
            }
          ),
          by = .(
            requested_g,
            class_a,
            class_b
          )
        ]
    }
  }

  clinical_pair <- rbindlist(
    clinical_pair_parts,
    use.names = TRUE,
    fill = TRUE
  )

  setorder(
    clinical_pair,
    requested_g,
    class_a,
    class_b
  )

  pair_context <- merge(
    copy(level_shape_joint),
    clinical_pair,
    by = c(
      "requested_g",
      "class_a",
      "class_b"
    ),
    all = TRUE,
    sort = TRUE
  )

  if (
    nrow(pair_context) !=
      nrow(level_shape_joint)
  ) {
    stop(
      paste0(
        "V292 fitted-point pair table and V294 clinical-anchor pair table ",
        "failed one-to-one merge."
      ),
      call. = FALSE
    )
  }

  pair_context[
    ,
    `:=`(
      total_distance_rank_fraction =
        rank_fraction(
          joint_total_RMS_z
        ),
      level_distance_rank_fraction =
        rank_fraction(
          joint_level_RMS_z
        ),
      fitted_shape_rank_fraction =
        rank_fraction(
          joint_centered_shape_RMS_z
        ),
      clinical_shape_rank_fraction =
        rank_fraction(
          clinical_joint_centered_shape_RMS_z
        ),
      clinical_total_rank_fraction =
        rank_fraction(
          clinical_joint_total_RMS_z
        ),
      shape_share_rank_fraction =
        rank_fraction(
          shape_share_of_level_plus_shape
        )
    ),
    by = requested_g
  ]

  pair_context[
    ,
    level_minus_fitted_shape_z :=
      joint_level_RMS_z -
      joint_centered_shape_RMS_z
  ]

  pair_context[
    ,
    level_minus_clinical_shape_z :=
      clinical_joint_level_RMS_z -
      clinical_joint_centered_shape_RMS_z
  ]

  nearest_parts <- list()
  np_idx <- 0L

  for (g in G_VALUES) {
    d <- pair_context[
      requested_g == g
    ]

    total_row <- d[
      order(
        joint_total_RMS_z,
        class_a,
        class_b
      )
    ][1L]

    shape_row <- d[
      order(
        joint_centered_shape_RMS_z,
        -mean_centered_shape_correlation,
        class_a,
        class_b
      )
    ][1L]

    clinical_shape_row <- d[
      order(
        clinical_joint_centered_shape_RMS_z,
        -clinical_mean_centered_shape_correlation,
        class_a,
        class_b
      )
    ][1L]

    np_idx <- np_idx + 1L

    nearest_parts[[np_idx]] <- data.table(
      requested_g = g,
      nearest_total_class_a =
        total_row$class_a,
      nearest_total_class_b =
        total_row$class_b,
      nearest_total_RMS_z =
        total_row$joint_total_RMS_z,
      nearest_total_level_RMS_z =
        total_row$joint_level_RMS_z,
      nearest_total_fitted_shape_RMS_z =
        total_row$joint_centered_shape_RMS_z,
      nearest_shape_class_a =
        shape_row$class_a,
      nearest_shape_class_b =
        shape_row$class_b,
      nearest_fitted_shape_RMS_z =
        shape_row$joint_centered_shape_RMS_z,
      nearest_shape_level_RMS_z =
        shape_row$joint_level_RMS_z,
      nearest_shape_shape_share =
        shape_row$shape_share_of_level_plus_shape,
      nearest_clinical_shape_class_a =
        clinical_shape_row$class_a,
      nearest_clinical_shape_class_b =
        clinical_shape_row$class_b,
      nearest_clinical_shape_RMS_z =
        clinical_shape_row$clinical_joint_centered_shape_RMS_z,
      nearest_clinical_shape_level_RMS_z =
        clinical_shape_row$clinical_joint_level_RMS_z,
      nearest_clinical_shape_share =
        clinical_shape_row$clinical_shape_share_of_level_plus_shape
    )
  }

  nearest_pair <- rbindlist(
    nearest_parts,
    use.names = TRUE,
    fill = TRUE
  )

  setorder(
    nearest_pair,
    requested_g
  )

  pair_long <- rbindlist(
    list(
      pair_context[
        ,
        .(
          requested_g,
          class = class_a,
          neighbour = class_b,
          joint_total_RMS_z,
          joint_level_RMS_z,
          joint_centered_shape_RMS_z,
          shape_share_of_level_plus_shape,
          clinical_joint_total_RMS_z,
          clinical_joint_level_RMS_z,
          clinical_joint_centered_shape_RMS_z,
          clinical_shape_share_of_level_plus_shape
        )
      ],
      pair_context[
        ,
        .(
          requested_g,
          class = class_b,
          neighbour = class_a,
          joint_total_RMS_z,
          joint_level_RMS_z,
          joint_centered_shape_RMS_z,
          shape_share_of_level_plus_shape,
          clinical_joint_total_RMS_z,
          clinical_joint_level_RMS_z,
          clinical_joint_centered_shape_RMS_z,
          clinical_shape_share_of_level_plus_shape
        )
      ]
    ),
    use.names = TRUE,
    fill = TRUE
  )

  nearest_total_by_class <-
    pair_long[
      order(
        requested_g,
        class,
        joint_total_RMS_z,
        neighbour
      ),
      .SD[1L],
      by = .(
        requested_g,
        class
      )
    ][
      ,
      .(
        requested_g,
        class,
        nearest_total_class =
          neighbour,
        nearest_total_RMS_z =
          joint_total_RMS_z,
        nearest_total_level_RMS_z =
          joint_level_RMS_z,
        nearest_total_shape_RMS_z =
          joint_centered_shape_RMS_z
      )
    ]

  nearest_shape_by_class <-
    pair_long[
      order(
        requested_g,
        class,
        joint_centered_shape_RMS_z,
        clinical_joint_centered_shape_RMS_z,
        neighbour
      ),
      .SD[1L],
      by = .(
        requested_g,
        class
      )
    ][
      ,
      .(
        requested_g,
        class,
        nearest_shape_class =
          neighbour,
        nearest_shape_RMS_z =
          joint_centered_shape_RMS_z,
        nearest_shape_level_RMS_z =
          joint_level_RMS_z,
        nearest_shape_shape_share =
          shape_share_of_level_plus_shape,
        nearest_shape_clinical_RMS_z =
          clinical_joint_centered_shape_RMS_z
      )
    ]

  class_nearest <- merge(
    nearest_total_by_class,
    nearest_shape_by_class,
    by = c(
      "requested_g",
      "class"
    ),
    all = TRUE,
    sort = TRUE
  )

  write_csv_exact(
    clinical_pair,
    clinical_pair_path
  )

  write_csv_exact(
    pair_context,
    pair_context_path
  )

  write_csv_exact(
    nearest_pair,
    nearest_pair_path
  )

  write_csv_exact(
    class_nearest,
    class_nearest_path
  )

  save_stage_marker(
    "stageB_pair_context",
    source_fingerprint,
    stage_b_outputs
  )
}

# =============================================================================
# 8. Stage C: lineage, parent split, child origin, refinement
# =============================================================================

parent_split_path <- file.path(
  table_dir,
  "V294_parent_split_G2G20.csv"
)

child_origin_path <- file.path(
  table_dir,
  "V294_child_origin_G2G20.csv"
)

refinement_path <- file.path(
  table_dir,
  "V294_refinement_G2G20.csv"
)

stage_c_outputs <- c(
  parent_split_path,
  child_origin_path,
  refinement_path
)

if (
  stage_is_valid(
    "stageC_refinement",
    source_fingerprint,
    stage_c_outputs
  )
) {
  log_message(
    "Stage C valid checkpoint found; reusing outputs"
  )

  parent_split <- fread(
    parent_split_path
  )

  child_origin <- fread(
    child_origin_path
  )

  refinement <- fread(
    refinement_path
  )
} else {
  log_message(
    "Stage C computing adjacent-G lineage and refinement"
  )

  parent_parts <- list()
  child_parts <- list()

  p_idx <- 0L
  c_idx <- 0L

  for (
    g in G_VALUES[
      -length(G_VALUES)
    ]
  ) {
    child_g <- g + 1L

    t <- transition[
      from_g == g &
        to_g == child_g
    ]

    if (nrow(t) == 0L) {
      stop(
        "No transition rows for G",
        g,
        " -> G",
        child_g,
        ".",
        call. = FALSE
      )
    }

    for (
      parent_class in
      sort(
        unique(
          t$class_from
        )
      )
    ) {
      d <- t[
        class_from == parent_class
      ][
        order(
          -row_pct,
          -N,
          class_to
        )
      ]

      p <- d$row_pct / 100

      top <- d[1L]

      second <- if (
        nrow(d) >= 2L
      ) {
        d[2L]
      } else {
        NULL
      }

      p_idx <- p_idx + 1L

      parent_parts[[p_idx]] <-
        data.table(
          from_g = g,
          to_g = child_g,
          parent_class =
            parent_class,
          parent_n =
            as.integer(
              top$row_total
            ),
          n_child_destinations =
            nrow(d),
          top_child_class =
            as.integer(
              top$class_to
            ),
          top_child_pct_of_parent =
            as.numeric(
              top$row_pct
            ),
          top_child_pct_of_child =
            as.numeric(
              top$column_pct
            ),
          second_child_class =
            if (!is.null(second)) {
              as.integer(
                second$class_to
              )
            } else {
              NA_integer_
            },
          second_child_pct_of_parent =
            if (!is.null(second)) {
              as.numeric(
                second$row_pct
              )
            } else {
              0
            },
          second_child_pct_of_child =
            if (!is.null(second)) {
              as.numeric(
                second$column_pct
              )
            } else {
              NA_real_
            },
          top2_combined_pct_of_parent =
            if (!is.null(second)) {
              as.numeric(
                top$row_pct +
                  second$row_pct
              )
            } else {
              as.numeric(
                top$row_pct
              )
            },
          split_entropy_normalised =
            normalised_entropy(p),
          split_effective_children =
            exp(
              -sum(
                p[p > 0] *
                  log(p[p > 0])
              )
            ),
          split_HHI =
            sum(p^2),
          n_children_ge_5pct =
            sum(
              d$row_pct >= 5
            ),
          n_children_ge_10pct =
            sum(
              d$row_pct >= 10
            ),
          n_children_ge_20pct =
            sum(
              d$row_pct >= 20
            )
        )
    }

    for (
      child_class in
      sort(
        unique(
          t$class_to
        )
      )
    ) {
      d <- t[
        class_to == child_class
      ][
        order(
          -column_pct,
          -N,
          class_from
        )
      ]

      p <- d$column_pct / 100

      top <- d[1L]

      second <- if (
        nrow(d) >= 2L
      ) {
        d[2L]
      } else {
        NULL
      }

      c_idx <- c_idx + 1L

      child_parts[[c_idx]] <-
        data.table(
          from_g = g,
          to_g = child_g,
          child_class =
            child_class,
          child_n =
            as.integer(
              top$column_total
            ),
          n_parent_sources =
            nrow(d),
          top_parent_class =
            as.integer(
              top$class_from
            ),
          top_parent_pct_of_child =
            as.numeric(
              top$column_pct
            ),
          top_parent_pct_of_parent =
            as.numeric(
              top$row_pct
            ),
          second_parent_class =
            if (!is.null(second)) {
              as.integer(
                second$class_from
              )
            } else {
              NA_integer_
            },
          second_parent_pct_of_child =
            if (!is.null(second)) {
              as.numeric(
                second$column_pct
              )
            } else {
              0
            },
          origin_entropy_normalised =
            normalised_entropy(p),
          origin_effective_parents =
            exp(
              -sum(
                p[p > 0] *
                  log(p[p > 0])
              )
            ),
          origin_HHI =
            sum(p^2)
        )
    }
  }

  parent_split <- rbindlist(
    parent_parts,
    use.names = TRUE,
    fill = TRUE
  )

  child_origin <- rbindlist(
    child_parts,
    use.names = TRUE,
    fill = TRUE
  )

  setorder(
    parent_split,
    from_g,
    parent_class
  )

  setorder(
    child_origin,
    from_g,
    child_class
  )

  refinement_parts <- list()
  r_idx <- 0L

  for (
    i in seq_len(
      nrow(parent_split)
    )
  ) {
    row <- parent_split[i]

    child_g <-
      as.integer(
        row$to_g
      )

    top_child <-
      as.integer(
        row$top_child_class
      )

    second_child <-
      as.integer(
        row$second_child_class
      )

    sibling <- NULL

    if (is.finite(second_child)) {
      sibling <- get_pair_row(
        pair_context,
        child_g,
        top_child,
        second_child
      )
    }

    sibling_total <-
      if (!is.null(sibling)) {
        sibling$joint_total_RMS_z
      } else {
        NA_real_
      }

    sibling_level <-
      if (!is.null(sibling)) {
        sibling$joint_level_RMS_z
      } else {
        NA_real_
      }

    sibling_fitted_shape <-
      if (!is.null(sibling)) {
        sibling$joint_centered_shape_RMS_z
      } else {
        NA_real_
      }

    sibling_clinical_shape <-
      if (!is.null(sibling)) {
        sibling$clinical_joint_centered_shape_RMS_z
      } else {
        NA_real_
      }

    sibling_shape_share <-
      if (!is.null(sibling)) {
        sibling$shape_share_of_level_plus_shape
      } else {
        NA_real_
      }

    sibling_total_rank <-
      if (!is.null(sibling)) {
        sibling$total_distance_rank_fraction
      } else {
        NA_real_
      }

    sibling_level_rank <-
      if (!is.null(sibling)) {
        sibling$level_distance_rank_fraction
      } else {
        NA_real_
      }

    sibling_fitted_shape_rank <-
      if (!is.null(sibling)) {
        sibling$fitted_shape_rank_fraction
      } else {
        NA_real_
      }

    sibling_clinical_shape_rank <-
      if (!is.null(sibling)) {
        sibling$clinical_shape_rank_fraction
      } else {
        NA_real_
      }

    sibling_shape_share_rank <-
      if (!is.null(sibling)) {
        sibling$shape_share_rank_fraction
      } else {
        NA_real_
      }

    substantive_split <-
      row$second_child_pct_of_parent >=
      SUBSTANTIVE_SPLIT_PCT

    minor_split <-
      row$second_child_pct_of_parent <
      MINOR_SPLIT_PCT

    complex_repartition_signal <-
      (
        row$top2_combined_pct_of_parent <
          70
      ) ||
      (
        row$split_effective_children >=
          2.5
      )

    near_shape_signal <-
      is.finite(
        sibling_fitted_shape_rank
      ) &&
      is.finite(
        sibling_clinical_shape_rank
      ) &&
      sibling_fitted_shape_rank <= 0.25 &&
      sibling_clinical_shape_rank <= 0.25

    level_dominant_signal <-
      near_shape_signal &&
      is.finite(
        sibling_level_rank
      ) &&
      sibling_level_rank >= 0.50 &&
      is.finite(
        sibling_shape_share_rank
      ) &&
      sibling_shape_share_rank <= 0.50

    near_total_signal <-
      is.finite(
        sibling_total_rank
      ) &&
      sibling_total_rank <= 0.25

    distinct_shape_signal <-
      is.finite(
        sibling_fitted_shape_rank
      ) &&
      is.finite(
        sibling_clinical_shape_rank
      ) &&
      is.finite(
        sibling_shape_share_rank
      ) &&
      sibling_fitted_shape_rank >= 0.50 &&
      sibling_clinical_shape_rank >= 0.50 &&
      sibling_shape_share_rank >= 0.50

    label <- "manual_review"

    if (minor_split) {
      label <-
        "minor_split_or_stable_parent"
    } else if (
      complex_repartition_signal
    ) {
      label <-
        "complex_repartition"
    } else if (
      substantive_split &&
      level_dominant_signal
    ) {
      label <-
        "mainly_level_stratification"
    } else if (
      substantive_split &&
      near_shape_signal &&
      near_total_signal &&
      !level_dominant_signal
    ) {
      label <-
        "probable_fragmentation"
    } else if (
      substantive_split &&
      distinct_shape_signal
    ) {
      label <-
        "meaningful_new_phenotype"
    } else if (
      substantive_split &&
      is.finite(
        sibling_fitted_shape_rank
      ) &&
      is.finite(
        sibling_clinical_shape_rank
      ) &&
      sibling_fitted_shape_rank >= 0.25 &&
      sibling_clinical_shape_rank >= 0.25
    ) {
      label <-
        "meaningful_subphenotype"
    }

    top_q <- class_quality[
      requested_g == child_g &
        class == top_child
    ]

    second_q <- if (
      is.finite(second_child)
    ) {
      class_quality[
        requested_g == child_g &
          class == second_child
      ]
    } else {
      data.table(
        class_pct = NA_real_,
        APPA_model = NA_real_
      )
    }

    r_idx <- r_idx + 1L

    refinement_parts[[r_idx]] <-
      data.table(
        from_g =
          row$from_g,
        to_g =
          row$to_g,
        parent_class =
          row$parent_class,
        parent_n =
          row$parent_n,
        top_child_class =
          top_child,
        second_child_class =
          second_child,
        top_child_pct_of_parent =
          row$top_child_pct_of_parent,
        second_child_pct_of_parent =
          row$second_child_pct_of_parent,
        top2_combined_pct_of_parent =
          row$top2_combined_pct_of_parent,
        split_entropy_normalised =
          row$split_entropy_normalised,
        split_effective_children =
          row$split_effective_children,
        split_HHI =
          row$split_HHI,
        sibling_joint_total_RMS_z =
          sibling_total,
        sibling_joint_level_RMS_z =
          sibling_level,
        sibling_fitted_shape_RMS_z =
          sibling_fitted_shape,
        sibling_clinical_shape_RMS_z =
          sibling_clinical_shape,
        sibling_shape_share =
          sibling_shape_share,
        sibling_total_rank_fraction =
          sibling_total_rank,
        sibling_level_rank_fraction =
          sibling_level_rank,
        sibling_fitted_shape_rank_fraction =
          sibling_fitted_shape_rank,
        sibling_clinical_shape_rank_fraction =
          sibling_clinical_shape_rank,
        sibling_shape_share_rank_fraction =
          sibling_shape_share_rank,
        top_child_class_pct =
          top_q$class_pct[1L],
        top_child_APPA =
          top_q$APPA_model[1L],
        second_child_class_pct =
          second_q$class_pct[1L],
        second_child_APPA =
          second_q$APPA_model[1L],
        complex_repartition_signal =
          as.integer(
            complex_repartition_signal
          ),
        near_shape_signal =
          as.integer(
            near_shape_signal
          ),
        level_dominant_signal =
          as.integer(
            level_dominant_signal
          ),
        near_total_signal =
          as.integer(
            near_total_signal
          ),
        distinct_shape_signal =
          as.integer(
            distinct_shape_signal
          ),
        refinement_review_label =
          label,
        label_is_multisignal_screening_only =
          1L
      )
  }

  refinement <- rbindlist(
    refinement_parts,
    use.names = TRUE,
    fill = TRUE
  )

  setorder(
    refinement,
    from_g,
    parent_class
  )

  write_csv_exact(
    parent_split,
    parent_split_path
  )

  write_csv_exact(
    child_origin,
    child_origin_path
  )

  write_csv_exact(
    refinement,
    refinement_path
  )

  save_stage_marker(
    "stageC_refinement",
    source_fingerprint,
    stage_c_outputs
  )
}

# =============================================================================
# 9. Stage D: G-level synthesis
# =============================================================================

g_summary_path <- file.path(
  table_dir,
  "V294_g_summary_G2G20.csv"
)

review_matrix_path <- file.path(
  table_dir,
  "V294_review_matrix_G2G20.csv"
)

lowg_path <- file.path(
  table_dir,
  "V294_lowG_focus_G2G6.csv"
)

stage_d_outputs <- c(
  g_summary_path,
  review_matrix_path,
  lowg_path
)

if (
  stage_is_valid(
    "stageD_g_summary",
    source_fingerprint,
    stage_d_outputs
  )
) {
  log_message(
    "Stage D valid checkpoint found; reusing outputs"
  )

  g_summary <- fread(
    g_summary_path
  )

  review_matrix <- fread(
    review_matrix_path
  )

  lowg_focus <- fread(
    lowg_path
  )
} else {
  log_message(
    "Stage D computing G-level structural synthesis"
  )

  quality_summary <- class_quality[
    ,
    .(
      class_count = .N,
      min_class_pct_recomputed =
        min(class_pct),
      median_class_pct =
        median(class_pct),
      n_class_lt_5pct =
        sum(class_pct < 5),
      pct_patients_in_classes_lt_5pct =
        sum(
          class_pct[
            class_pct < 5
          ]
        ),
      min_APPA_recomputed =
        min(APPA_model),
      median_APPA =
        median(APPA_model),
      mean_APPA =
        mean(APPA_model),
      min_mean_max_posterior_in_class =
        min(
          mean_max_posterior_in_assigned_class
        ),
      median_mean_max_posterior_in_class =
        median(
          mean_max_posterior_in_assigned_class
        )
    ),
    by = requested_g
  ]

  pair_summary <- pair_context[
    ,
    .(
      median_joint_total_RMS_z =
        median(joint_total_RMS_z),
      median_joint_level_RMS_z =
        median(joint_level_RMS_z),
      median_fitted_shape_RMS_z =
        median(
          joint_centered_shape_RMS_z
        ),
      median_clinical_shape_RMS_z =
        median(
          clinical_joint_centered_shape_RMS_z
        ),
      median_shape_share =
        median(
          shape_share_of_level_plus_shape
        ),
      q25_fitted_shape_RMS_z =
        as.numeric(
          quantile(
            joint_centered_shape_RMS_z,
            0.25,
            names = FALSE
          )
        ),
      q25_clinical_shape_RMS_z =
        as.numeric(
          quantile(
            clinical_joint_centered_shape_RMS_z,
            0.25,
            names = FALSE
          )
        ),
      proportion_pairs_level_gt_fitted_shape =
        mean(
          joint_level_RMS_z >
            joint_centered_shape_RMS_z
        ),
      proportion_pairs_level_gt_clinical_shape =
        mean(
          clinical_joint_level_RMS_z >
            clinical_joint_centered_shape_RMS_z
        )
    ),
    by = requested_g
  ]

  refinement_counts <- refinement[
    ,
    .(
      n_parent_rows = .N,
      n_meaningful_new_phenotype =
        sum(
          refinement_review_label ==
            "meaningful_new_phenotype"
        ),
      n_meaningful_subphenotype =
        sum(
          refinement_review_label ==
            "meaningful_subphenotype"
        ),
      n_mainly_level_stratification =
        sum(
          refinement_review_label ==
            "mainly_level_stratification"
        ),
      n_probable_fragmentation =
        sum(
          refinement_review_label ==
            "probable_fragmentation"
        ),
      n_complex_repartition =
        sum(
          refinement_review_label ==
            "complex_repartition"
        ),
      n_minor_split_or_stable_parent =
        sum(
          refinement_review_label ==
            "minor_split_or_stable_parent"
        ),
      n_manual_review =
        sum(
          refinement_review_label ==
            "manual_review"
        ),
      mean_parent_split_entropy =
        mean(
          split_entropy_normalised
        ),
      mean_parent_effective_children =
        mean(
          split_effective_children
        )
    ),
    by = .(
      requested_g = to_g
    )
  ]

  inbound_metrics <- transition_metrics[
    ,
    .(
      requested_g = to_g,
      inbound_from_g = from_g,
      inbound_ARI =
        adjusted_rand_index,
      inbound_NMI =
        normalised_mutual_information
    )
  ]

  outbound_metrics <- transition_metrics[
    ,
    .(
      requested_g = from_g,
      outbound_to_g = to_g,
      outbound_ARI =
        adjusted_rand_index,
      outbound_NMI =
        normalised_mutual_information
    )
  ]

  g_summary <- Reduce(
    function(x, y) {
      merge(
        x,
        y,
        by = "requested_g",
        all = TRUE,
        sort = TRUE
      )
    },
    list(
      copy(model_summary),
      quality_summary,
      nearest_pair,
      pair_summary,
      inbound_metrics,
      outbound_metrics,
      refinement_counts
    )
  )

  setorder(
    g_summary,
    requested_g
  )

  g_summary[
    ,
    `:=`(
      delta_logLik_vs_previous_G =
        final_logLik -
        shift(
          final_logLik,
          1L
        ),
      delta_BIC_vs_previous_G =
        BIC -
        shift(
          BIC,
          1L
        ),
      relative_BIC_change_vs_previous_G =
        (
          BIC -
            shift(
              BIC,
              1L
            )
        ) /
        abs(
          shift(
            BIC,
            1L
          )
        )
    )
  ]

  # Descriptive fragmentation evidence columns.
  # These are intentionally NOT collapsed into one automated G score.
  g_summary[
    ,
    `:=`(
      fragmentation_small_class_signal =
        as.integer(
          n_class_lt_5pct > 0
        ),
      fragmentation_near_shape_signal =
        as.integer(
          rank_fraction(
            nearest_clinical_shape_RMS_z
          ) <= 0.25
        ),
      fragmentation_repartition_signal =
        as.integer(
          fifelse(
            is.na(n_complex_repartition),
            0,
            n_complex_repartition
          ) > 0
        ),
      level_stratification_signal =
        as.integer(
          fifelse(
            is.na(
              n_mainly_level_stratification
            ),
            0,
            n_mainly_level_stratification
          ) > 0
        )
    )
  ]

  review_matrix <- g_summary[
    ,
    .(
      requested_g,
      final_logLik,
      BIC,
      delta_logLik_vs_previous_G,
      delta_BIC_vs_previous_G,
      min_class_pct,
      n_class_lt_5pct,
      pct_patients_in_classes_lt_5pct,
      min_APPA,
      median_APPA,
      mean_max_posterior,
      median_max_posterior,
      mean_posterior_entropy_normalised,
      nearest_total_RMS_z,
      nearest_fitted_shape_RMS_z,
      nearest_clinical_shape_RMS_z,
      median_joint_level_RMS_z,
      median_fitted_shape_RMS_z,
      median_clinical_shape_RMS_z,
      median_shape_share,
      proportion_pairs_level_gt_fitted_shape,
      proportion_pairs_level_gt_clinical_shape,
      inbound_ARI,
      inbound_NMI,
      outbound_ARI,
      outbound_NMI,
      n_meaningful_new_phenotype,
      n_meaningful_subphenotype,
      n_mainly_level_stratification,
      n_probable_fragmentation,
      n_complex_repartition,
      n_manual_review,
      fragmentation_small_class_signal,
      fragmentation_near_shape_signal,
      fragmentation_repartition_signal,
      level_stratification_signal
    )
  ]

  lowg_focus <- review_matrix[
    requested_g %in%
      LOW_G_FOCUS
  ]

  write_csv_exact(
    g_summary,
    g_summary_path
  )

  write_csv_exact(
    review_matrix,
    review_matrix_path
  )

  write_csv_exact(
    lowg_focus,
    lowg_path
  )

  save_stage_marker(
    "stageD_g_summary",
    source_fingerprint,
    stage_d_outputs
  )
}

# =============================================================================
# 10. Class-level review table
# =============================================================================

class_review_path <- file.path(
  table_dir,
  "V294_class_review_G2G20.csv"
)

class_level_base <- class_quality[
  ,
  .(
    requested_g,
    class,
    class_label,
    n,
    class_pct,
    APPA_model,
    posterior_soft_mass,
    posterior_soft_mass_pct,
    mean_max_posterior_in_assigned_class,
    median_max_posterior_in_assigned_class
  )
]

class_review <- merge(
  class_level_base,
  class_nearest,
  by = c(
    "requested_g",
    "class"
  ),
  all = TRUE,
  sort = TRUE
)

class_structure_wide <- dcast(
  class_structure,
  requested_g + class ~ indicator,
  value.var = c(
    "start_native",
    "end_native",
    "delta_native",
    "mean_native",
    "delta_z",
    "clinical_centered_shape_RMS_z",
    "total_variation_z",
    "path_excess_z",
    "turn_count",
    "curvature_RMS_z"
  ),
  sep = "_"
)

class_review <- merge(
  class_review,
  class_structure_wide,
  by = c(
    "requested_g",
    "class"
  ),
  all = TRUE,
  sort = TRUE
)

class_review <- merge(
  class_review,
  class_coupling,
  by = c(
    "requested_g",
    "class"
  ),
  all = TRUE,
  sort = TRUE
)

setorder(
  class_review,
  requested_g,
  class
)

write_csv_exact(
  class_review,
  class_review_path
)

# =============================================================================
# 11. Copy key V292 figures for the adjudication bundle
# =============================================================================

if (dir.exists(source_figure_dir)) {
  expected_figure_names <- c(
    "V292_native_all_G2G20.pdf",
    "V292_centered_all_G2G20.pdf",
    unlist(
      lapply(
        INDICATORS,
        function(indicator_now) {
          c(
            paste0(
              "V292_compare_native_",
              indicator_now,
              "_G2G20.png"
            ),
            paste0(
              "V292_compare_centered_",
              indicator_now,
              "_G2G20.png"
            ),
            paste0(
              "V292_compare_native_",
              indicator_now,
              "_G2G6.png"
            ),
            paste0(
              "V292_compare_centered_",
              indicator_now,
              "_G2G6.png"
            )
          )
        }
      ),
      use.names = FALSE
    )
  )

  existing_copies <- list.files(
    figure_dir,
    full.names = TRUE,
    recursive = TRUE,
    all.files = TRUE,
    include.dirs = FALSE
  )

  if (length(existing_copies) > 0L) {
    unlink(
      existing_copies,
      force = TRUE
    )
  }

  for (
    file_name in
    expected_figure_names
  ) {
    src <- file.path(
      source_figure_dir,
      file_name
    )

    if (file.exists(src)) {
      dst <- file.path(
        figure_dir,
        file_name
      )

      ok <- file.copy(
        src,
        dst,
        overwrite = TRUE,
        copy.mode = TRUE
      )

      if (!isTRUE(ok)) {
        stop(
          "Failed to copy V292 figure: ",
          src,
          call. = FALSE
        )
      }

      if (!identical(
        md5_one(src),
        md5_one(dst)
      )) {
        stop(
          "MD5 mismatch after copying V292 figure: ",
          src,
          call. = FALSE
        )
      }
    }
  }
}

# =============================================================================
# 12. Human-readable report
# =============================================================================

lowg_display <- lowg_focus[
  ,
  .(
    requested_g,
    BIC,
    delta_BIC_vs_previous_G,
    min_class_pct,
    min_APPA,
    mean_max_posterior,
    mean_posterior_entropy_normalised,
    nearest_clinical_shape_RMS_z,
    median_joint_level_RMS_z,
    median_clinical_shape_RMS_z,
    inbound_ARI,
    inbound_NMI,
    n_meaningful_new_phenotype,
    n_meaningful_subphenotype,
    n_mainly_level_stratification,
    n_probable_fragmentation,
    n_complex_repartition
  )
]

refinement_counts_display <- refinement[
  ,
  .N,
  by = .(
    to_g,
    refinement_review_label
  )
][
  order(
    to_g,
    refinement_review_label
  )
]

report_lines <- c(
  "FOURTH STUDY / MSI-Electrolyte GBMTM",
  "V294 revised local Stage 3 structural adjudication: G2-G20",
  "",
  paste0(
    "Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste0(
    "V292 source root: ",
    normalizePath(
      source_root,
      winslash = "/",
      mustWork = TRUE
    )
  ),
  "",
  "SCIENTIFIC SCOPE",
  "- No GBMTM refitting.",
  "- No MAKE30, mortality, KRT, or any outcome information.",
  "- G2-G20 are compared continuously.",
  "- Persistent level and temporal shape are explicitly separated.",
  "- Full fitted-point centered shape is retained for model-structure comparison.",
  "- Electrolyte clinical shape evidence uses only 3/9/15/21 h anchors.",
  "- MSI clinical structure uses 1,3,5,...,23 h.",
  "- Multivariate coupling is described at the common 3/9/15/21 h anchors.",
  "- Refinement labels are multi-signal screening labels only.",
  "- No single fixed RMSE cutoff determines phenotype novelty.",
  "- No automatic final G or Step-4A shortlist is frozen by this script.",
  "",
  "PRIMARY OUTPUTS",
  "- V294_g_summary_G2G20.csv",
  "- V294_review_matrix_G2G20.csv",
  "- V294_lowG_focus_G2G6.csv",
  "- V294_class_review_G2G20.csv",
  "- V294_class_anchor_G2G20.csv",
  "- V294_class_structure_G2G20.csv",
  "- V294_class_coupling_G2G20.csv",
  "- V294_pair_context_G2G20.csv",
  "- V294_nearest_pair_G2G20.csv",
  "- V294_class_nearest_G2G20.csv",
  "- V294_parent_split_G2G20.csv",
  "- V294_child_origin_G2G20.csv",
  "- V294_refinement_G2G20.csv",
  "",
  "LOW-G FOCUS: G2-G6",
  paste(
    capture.output(
      print(
        lowg_display
      )
    ),
    collapse = "\n"
  ),
  "",
  "REFINEMENT LABEL COUNTS BY CHILD G",
  paste(
    capture.output(
      print(
        refinement_counts_display
      )
    ),
    collapse = "\n"
  ),
  "",
  "FINAL INTERPRETATION RULE",
  paste0(
    "The final shortlist must be selected after joint review of fit quality, ",
    "posterior quality, class size, native level, centered shape, clinical-anchor ",
    "shape, multivariate coupling, lineage/reclassification, and fragmentation. ",
    "A lower BIC alone is insufficient, and a persistent level offset alone is ",
    "not automatically a new dynamic phenotype."
  ),
  "",
  "Upload V294_stage3_bundle_G2G20.zip back to ChatGPT for final G2-G20",
  "structural adjudication and selection of 2-3 G values for Step 4A."
)

writeLines(
  report_lines,
  report_path,
  useBytes = TRUE
)

# =============================================================================
# 13. Session info, output manifest, ZIP
# =============================================================================

session_path <- file.path(
  output_dir,
  "V294_sessionInfo_G2G20.txt"
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  session_path,
  useBytes = TRUE
)

readme_path <- file.path(
  output_dir,
  "V294_README_G2G20.txt"
)

readme_lines <- c(
  "V294 Stage 3 G2-G20",
  "",
  "Run V292_step2_visual_G2G20.R first.",
  "Then run this V294 script.",
  "",
  "V294 reads only the V292 compact CSV outputs.",
  "All Stage-3 results are written under GBMTM_stage3_output.",
  "It does not reopen the 12+ GB archived GBMTM RDS set.",
  "It does not refit any model and does not read any outcome.",
  "",
  "Upload:",
  "V294_stage3_bundle_G2G20.zip",
  "",
  "for final structural adjudication and Step-4A shortlist selection."
)

writeLines(
  readme_lines,
  readme_path,
  useBytes = TRUE
)

manifest_path <- file.path(
  table_dir,
  "V294_output_manifest_G2G20.csv"
)

manifest_files <- list.files(
  output_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE,
  include.dirs = FALSE
)

manifest_files <- manifest_files[
  !grepl(
    "[\\\\/]_resume[\\\\/]",
    manifest_files
  ) &
    basename(manifest_files) !=
      basename(zip_path) &
    normalizePath(
      manifest_files,
      winslash = "/",
      mustWork = TRUE
    ) !=
      normalizePath(
        manifest_path,
        winslash = "/",
        mustWork = FALSE
      )
]

root_norm <- normalizePath(
  output_dir,
  winslash = "/",
  mustWork = TRUE
)

output_manifest <- data.table(
  relative_path = substring(
    normalizePath(
      manifest_files,
      winslash = "/",
      mustWork = TRUE
    ),
    nchar(root_norm) + 2L
  ),
  size_bytes = as.numeric(
    file.info(
      manifest_files
    )$size
  ),
  md5 = vapply(
    manifest_files,
    md5_one,
    character(1)
  )
)

setorder(
  output_manifest,
  relative_path
)

write_csv_exact(
  output_manifest,
  manifest_path
)

if (file.exists(zip_path)) {
  unlink(
    zip_path,
    force = TRUE
  )
}

old_wd <- getwd()

tryCatch(
  {
    setwd(output_dir)

    upload_items <- list.files(
      output_dir,
      all.files = FALSE,
      full.names = FALSE
    )

    upload_items <- upload_items[
      upload_items != basename(zip_path) &
        upload_items != basename(resume_dir)
    ]

    zip::zipr(
      zipfile = basename(zip_path),
      files = upload_items,
      recurse = TRUE
    )
  },
  finally = {
    setwd(old_wd)
  }
)

if (
  !file.exists(zip_path) ||
  is.na(file.info(zip_path)$size) ||
  file.info(zip_path)$size <= 0
) {
  stop(
    "Final Stage-3 ZIP creation failed: ",
    zip_path,
    call. = FALSE
  )
}

log_message(
  "V294 Stage 3 completed successfully."
)

log_message(
  "Final ZIP: ",
  normalizePath(
    zip_path,
    winslash = "/",
    mustWork = TRUE
  )
)

cat("\n")
cat("============================================================\n")
cat("V294 Stage 3 completed successfully.\n")
cat("Output directory:\n")
cat(
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)
cat("Upload this ZIP back to ChatGPT:\n")
cat(
  normalizePath(
    zip_path,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)
cat("============================================================\n")
