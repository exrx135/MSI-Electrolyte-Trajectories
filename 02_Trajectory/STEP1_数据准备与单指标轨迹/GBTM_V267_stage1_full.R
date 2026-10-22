# ============================================================
# V267 - STAGE 1 SINGLE-INDICATOR GBTM SCREEN
# ============================================================
#
#
# V267 corrections relative to V263:
#   - Uses stored fit$fitted trajectories directly for trajectory export.
#   - Preserves exact original model time-level labels.
#   - Audits trajectory row counts and finite fitted values.
#   - Correctly recognizes G1 unnamed V1/V2/V3 EM columns.
#   - Keeps the validated V263 PSOCK direct-argument architecture.
#
# Research role:
#   Stage 1 exploratory single-indicator GBTM.
#
# Indicators:
#   MSI, Na, K, Cl, HCO3
#
# Models:
#   Each indicator is modeled SEPARATELY.
#   G1-G6 are fitted for each indicator.
#
# Purpose:
#   Identify the major single-indicator trajectory shapes that will
#   later help interpret Stage-2 five-indicator joint GBMTM solutions.
#
# IMPORTANT:
#   - Stage 1 classes are NOT final joint phenotypes.
#   - Do NOT Cartesian-product single-indicator classes.
#   - MAKE30 or any outcome is NOT used here.
#
# Modeling settings:
#   - G1-G6 for each indicator
#   - degree d = 2
#   - scaling = 0 because input is already globally standardized
#   - pruning = TRUE
#   - tol = 1 for Stage-1 exploratory screen
#   - G1: nstart = 1
#   - G2-G6: nstart = 3 random starts; best LL retained by gbmt
#   - delete.empty = FALSE
#   - no automatic rerun to force non-empty classes
#   - no wall-clock timeout
#   - maxit uses .Machine$integer.max only because gbmt requires
#     a finite positive integer; it is intended to be non-binding.
#
# Missingness handling for Stage 1:
#   The prepared joint input has 12 rows/patient because MSI uses
#   12 x 2-hour bins, while electrolytes use 4 x 6-hour bins.
#
#   For each SINGLE-indicator model, rows where that indicator is NA
#   are removed BEFORE gbmt(). Therefore:
#     MSI  -> actual observed 2-hour bins only
#     Na   -> actual observed 6-hour bins only
#     K    -> actual observed 6-hour bins only
#     Cl   -> actual observed 6-hour bins only
#     HCO3 -> actual observed 6-hour bins only
#
#   No deterministic interpolation, LOCF, NOCB, or external imputation
#   is added by this script.
#
# Windows path rule:
#   Never hard-code a Windows username.
#   USERPROFILE is used on Windows.
#
# Exact base hierarchy:
#   <USERPROFILE>/Desktop/Research/MSI/提取轨迹文件/GBMTM_stage1_output
#
# Output naming:
#   Version_Content_Indicator_GX
#
# Examples:
#   V267_task_summary_MSI_G3.csv
#   V267_class_quality_Na_G4.csv
#   V267_trajectory_K_G5.csv
#   V267_screen_fit_HCO3_G6.rds
#
# ============================================================


options(stringsAsFactors = FALSE)


# ============================================================
# 0. Fixed configuration
# ============================================================

VERSION_TAG <- "V267"

GROUPS_TO_SCREEN <- 1:6

MODEL_DEGREE <- 2L

SCREEN_TOL <- 1

NSTART_G1 <- 1L
NSTART_G2PLUS <- 3L

N_WORKERS <- 8L

START_STAGGER_SECONDS <- 5L

NONBINDING_MAXIT <- .Machine$integer.max

MASTER_SEED <- 20260905L

EXPECTED_PATIENTS <- 9804L
EXPECTED_LONG_ROWS <- 117648L

PREFERRED_MIN_CLASS_PCT <- 5
PREFERRED_MIN_APPA <- 0.70

INDICATOR_MAP <- data.frame(
  indicator = c(
    "MSI",
    "Na",
    "K",
    "Cl",
    "HCO3"
  ),
  variable = c(
    "msi_z",
    "sodium_z",
    "potassium_z",
    "chloride_z",
    "bicarbonate_z"
  ),
  expected_min_points = c(
    9L,
    2L,
    2L,
    2L,
    2L
  ),
  expected_max_points = c(
    12L,
    4L,
    4L,
    4L,
    4L
  ),
  stringsAsFactors = FALSE
)

REQUIRED_COLS <- c(
  "unit_id",
  "time_model",
  INDICATOR_MAP$variable
)


# ============================================================
# 1. Required packages
# ============================================================

required_packages <- c(
  "data.table",
  "gbmt"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before running V267."
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(gbmt)
  library(parallel)
})


# ============================================================
# 2. Prevent nested threading inside PSOCK workers
# ============================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

data.table::setDTthreads(1L)


# ============================================================
# 3. Verify gbmt interface and version
# ============================================================

installed_gbmt_version <- as.character(
  packageVersion("gbmt")
)

cat(
  "Installed gbmt version: ",
  installed_gbmt_version,
  "\n",
  sep = ""
)

gbmt_formals <- names(
  formals(
    gbmt::gbmt
  )
)

required_gbmt_args <- c(
  "x.names",
  "unit",
  "time",
  "ng",
  "d",
  "data",
  "scaling",
  "pruning",
  "nstart",
  "tol",
  "maxit",
  "delete.empty"
)

missing_gbmt_args <- setdiff(
  required_gbmt_args,
  gbmt_formals
)

if (length(missing_gbmt_args) > 0L) {
  stop(
    paste0(
      "Installed gbmt::gbmt() does not expose required argument(s): ",
      paste(missing_gbmt_args, collapse = ", "),
      "\nInstalled gbmt version: ",
      installed_gbmt_version,
      "\nStop here. Do not silently run with a different interface."
    )
  )
}


# ============================================================
# 4. Resolve portable path
# ============================================================

if (.Platform$OS.type == "windows") {

  user_home <- Sys.getenv(
    "USERPROFILE",
    unset = ""
  )

  if (!nzchar(user_home)) {
    stop(
      "Windows USERPROFILE is empty. The script will not guess another user directory."
    )
  }

} else {

  user_home <- path.expand("~")

}

if (!dir.exists(user_home)) {
  stop(
    paste0(
      "Resolved user profile directory does not exist:\n",
      normalizePath(
        user_home,
        winslash = "/",
        mustWork = FALSE
      )
    )
  )
}

base_dir <- file.path(
  user_home,
  "Desktop",
  "Research",
  "MSI",
  "提取轨迹文件",
  "GBMTM_stage1_output"
)

input_rds <- file.path(
  base_dir,
  "MIMIC_GBMTM_long_scaled.rds"
)

output_dir <- file.path(
  base_dir,
  "GBTM_V267_single_screen"
)

if (!dir.exists(base_dir)) {
  stop(
    paste0(
      "Required base directory does not exist.\n",
      "Checked exactly:\n",
      normalizePath(
        base_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      "\nThe script will not guess OneDrive, another drive, ",
      "another username, or another folder hierarchy."
    )
  )
}

if (!file.exists(input_rds)) {
  stop(
    paste0(
      "Required input file does not exist.\n",
      "Checked exactly:\n",
      normalizePath(
        input_rds,
        winslash = "/",
        mustWork = FALSE
      )
    )
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(output_dir)) {
  stop(
    paste0(
      "Could not create output directory:\n",
      normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
      )
    )
  )
}


# ============================================================
# 5. Console header
# ============================================================

cat("\n")
cat("============================================================\n")
cat("V267 STAGE 1 SINGLE-INDICATOR GBTM SCREEN\n")
cat("============================================================\n")
cat(
  "Input : ",
  normalizePath(
    input_rds,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)
cat(
  "Output: ",
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)
cat("Indicators              : MSI, Na, K, Cl, HCO3\n")
cat("G range                 : G1-G6 each indicator\n")
cat("Workers                 : 8\n")
cat("Polynomial degree       : 2\n")
cat("Screening tolerance     : 1\n")
cat("G1 nstart               : 1\n")
cat("G2-G6 nstart            : 3\n")
cat("delete.empty            : FALSE\n")
cat("Automatic rerun         : NONE\n")
cat("Wall-clock timeout      : NONE\n")
cat("Practical maxit budget  : NONE\n")
cat(
  "Technical maxit        : ",
  NONBINDING_MAXIT,
  "\n",
  sep = ""
)
cat("============================================================\n\n")


# ============================================================
# 6. Read and validate prepared joint input
# ============================================================

cat("Reading prepared MIMIC input...\n")

full_dt <- as.data.table(
  readRDS(
    input_rds
  )
)

missing_cols <- setdiff(
  REQUIRED_COLS,
  names(full_dt)
)

if (length(missing_cols) > 0L) {
  stop(
    paste0(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

full_dt[
  ,
  unit_id := as.character(unit_id)
]

full_dt[
  ,
  time_model := as.numeric(time_model)
]

for (v in INDICATOR_MAP$variable) {
  full_dt[
    ,
    (v) := as.numeric(get(v))
  ]
}

n_patients <- uniqueN(
  full_dt$unit_id
)

n_rows <- nrow(
  full_dt
)

if (n_patients != EXPECTED_PATIENTS) {
  stop(
    paste0(
      "Unexpected patient count: ",
      n_patients,
      "; expected ",
      EXPECTED_PATIENTS,
      "."
    )
  )
}

if (n_rows != EXPECTED_LONG_ROWS) {
  stop(
    paste0(
      "Unexpected row count: ",
      n_rows,
      "; expected ",
      EXPECTED_LONG_ROWS,
      "."
    )
  )
}

if (anyNA(full_dt$unit_id)) {
  stop("unit_id contains missing values.")
}

if (anyNA(full_dt$time_model)) {
  stop("time_model contains missing values.")
}

duplicate_time_rows <- full_dt[
  ,
  .N,
  by = .(
    unit_id,
    time_model
  )
][
  N > 1L
]

if (nrow(duplicate_time_rows) > 0L) {
  stop(
    "Duplicate unit_id + time_model rows detected."
  )
}

cat(
  "Prepared joint input validated: ",
  n_patients,
  " patients; ",
  n_rows,
  " rows.\n\n",
  sep = ""
)


# ============================================================
# 7. Create and validate one observed-only input per indicator
# ============================================================

indicator_input_info <- list()

for (
  i in seq_len(
    nrow(INDICATOR_MAP)
  )
) {

  indicator_name <- INDICATOR_MAP$indicator[i]
  variable_name <- INDICATOR_MAP$variable[i]

  expected_min_points <- INDICATOR_MAP$expected_min_points[i]
  expected_max_points <- INDICATOR_MAP$expected_max_points[i]

  cat(
    "Preparing observed-only input for ",
    indicator_name,
    "...\n",
    sep = ""
  )

  indicator_dt <- full_dt[
    !is.na(
      get(variable_name)
    ),
    .(
      unit_id,
      time_model,
      value = get(variable_name)
    )
  ]

  setnames(
    indicator_dt,
    "value",
    variable_name
  )

  if (nrow(indicator_dt) == 0L) {
    stop(
      paste0(
        indicator_name,
        " has no observed rows."
      )
    )
  }

  if (
    anyNA(
      indicator_dt[[variable_name]]
    )
  ) {
    stop(
      paste0(
        indicator_name,
        " observed-only input still contains NA values."
      )
    )
  }

  indicator_patients <- uniqueN(
    indicator_dt$unit_id
  )

  if (indicator_patients != EXPECTED_PATIENTS) {
    stop(
      paste0(
        indicator_name,
        " observed-only input contains ",
        indicator_patients,
        " patients; expected ",
        EXPECTED_PATIENTS,
        ".\n",
        "The prepared cohort is expected to satisfy this indicator's ",
        "trajectory eligibility rule for every included patient."
      )
    )
  }

  point_count_dt <- indicator_dt[
    ,
    .(
      n_points = .N
    ),
    by = unit_id
  ]

  actual_min_points <- min(
    point_count_dt$n_points
  )

  actual_max_points <- max(
    point_count_dt$n_points
  )

  if (
    actual_min_points <
      expected_min_points
  ) {
    stop(
      paste0(
        indicator_name,
        " has a patient with only ",
        actual_min_points,
        " observed point(s); expected minimum ",
        expected_min_points,
        "."
      )
    )
  }

  if (
    actual_max_points >
      expected_max_points
  ) {
    stop(
      paste0(
        indicator_name,
        " has a patient with ",
        actual_max_points,
        " observed points; expected maximum ",
        expected_max_points,
        "."
      )
    )
  }

  duplicate_indicator_time <- indicator_dt[
    ,
    .N,
    by = .(
      unit_id,
      time_model
    )
  ][
    N > 1L
  ]

  if (
    nrow(
      duplicate_indicator_time
    ) > 0L
  ) {
    stop(
      paste0(
        indicator_name,
        " observed-only input contains duplicate unit_id + time_model rows."
      )
    )
  }

  indicator_file <- file.path(
    output_dir,
    paste0(
      "V267_model_input_",
      indicator_name,
      ".rds"
    )
  )

  saveRDS(
    as.data.frame(
      indicator_dt
    ),
    indicator_file,
    compress = FALSE
  )

  if (!file.exists(indicator_file)) {
    stop(
      paste0(
        "Failed to write observed-only input for ",
        indicator_name,
        "."
      )
    )
  }

  indicator_input_info[[indicator_name]] <- data.table(
    indicator = indicator_name,
    variable = variable_name,
    input_file = indicator_file,
    n_patients = indicator_patients,
    n_rows = nrow(indicator_dt),
    min_points_per_patient = actual_min_points,
    max_points_per_patient = actual_max_points,
    expected_min_points = expected_min_points,
    expected_max_points = expected_max_points,
    unique_time_points = uniqueN(
      indicator_dt$time_model
    ),
    min_time_model = min(
      indicator_dt$time_model
    ),
    max_time_model = max(
      indicator_dt$time_model
    ),
    file_size_mb = file.info(
      indicator_file
    )$size / 1024^2
  )

  cat(
    "  ",
    indicator_name,
    ": ",
    nrow(indicator_dt),
    " observed rows; ",
    actual_min_points,
    "-",
    actual_max_points,
    " points/patient; ",
    uniqueN(indicator_dt$time_model),
    " unique time points.\n",
    sep = ""
  )

}

indicator_input_dt <- rbindlist(
  indicator_input_info,
  fill = TRUE
)

fwrite(
  indicator_input_dt,
  file.path(
    output_dir,
    "V267_indicator_input_summary.csv"
  )
)

rm(
  full_dt,
  indicator_dt,
  point_count_dt,
  duplicate_time_rows,
  duplicate_indicator_time
)

gc(full = TRUE)

cat("\nAll five observed-only indicator inputs validated.\n\n")


# ============================================================
# 8. Save run configuration
# ============================================================

logical_cores <- parallel::detectCores(
  logical = TRUE
)

physical_cores <- parallel::detectCores(
  logical = FALSE
)

config_dt <- data.table(
  parameter = c(
    "VERSION_TAG",
    "input_rds",
    "output_dir",
    "GROUPS_TO_SCREEN",
    "MODEL_DEGREE",
    "SCREEN_TOL",
    "NSTART_G1",
    "NSTART_G2PLUS",
    "N_WORKERS",
    "START_STAGGER_SECONDS",
    "NONBINDING_MAXIT",
    "DELETE_EMPTY",
    "PRUNING",
    "MASTER_SEED",
    "EXPECTED_PATIENTS",
    "EXPECTED_LONG_ROWS",
    "PREFERRED_MIN_CLASS_PCT",
    "PREFERRED_MIN_APPA",
    "logical_cores",
    "physical_cores",
    "R_version",
    "gbmt_version",
    "data.table_version"
  ),
  value = c(
    VERSION_TAG,
    input_rds,
    output_dir,
    paste(
      GROUPS_TO_SCREEN,
      collapse = ","
    ),
    as.character(MODEL_DEGREE),
    as.character(SCREEN_TOL),
    as.character(NSTART_G1),
    as.character(NSTART_G2PLUS),
    as.character(N_WORKERS),
    as.character(START_STAGGER_SECONDS),
    as.character(NONBINDING_MAXIT),
    "FALSE",
    "TRUE",
    as.character(MASTER_SEED),
    as.character(EXPECTED_PATIENTS),
    as.character(EXPECTED_LONG_ROWS),
    as.character(PREFERRED_MIN_CLASS_PCT),
    as.character(PREFERRED_MIN_APPA),
    as.character(logical_cores),
    as.character(physical_cores),
    R.version.string,
    installed_gbmt_version,
    as.character(
      packageVersion("data.table")
    )
  )
)

fwrite(
  config_dt,
  file.path(
    output_dir,
    "V267_config.csv"
  )
)


# ============================================================
# 9. Helper functions
# ============================================================

append_stage <- function(
  path,
  text
) {

  cat(
    paste0(
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      " | ",
      text,
      "\n"
    ),
    file = path,
    append = TRUE
  )

}


safe_ic_value <- function(
  fit,
  target
) {

  x <- fit$ic

  if (
    is.null(x) ||
    is.null(names(x))
  ) {
    return(NA_real_)
  }

  idx <- which(
    tolower(
      names(x)
    ) ==
      tolower(target)
  )

  if (length(idx) == 0L) {
    return(NA_real_)
  }

  as.numeric(
    x[
      idx[1L]
    ]
  )

}


extract_em_table <- function(
  fit
) {

  x <- fit$em

  if (is.null(x)) {
    return(
      data.table()
    )
  }

  if (is.null(dim(x))) {
    x <- matrix(
      x,
      nrow = 1L
    )
  }

  out <- as.data.table(
    as.data.frame(
      x,
      stringsAsFactors = FALSE
    )
  )

  out[
    ,
    restart := seq_len(.N)
  ]

  setcolorder(
    out,
    c(
      "restart",
      setdiff(
        names(out),
        "restart"
      )
    )
  )

  out

}


extract_em_status <- function(
  fit
) {

  em_dt <- extract_em_table(
    fit
  )

  if (nrow(em_dt) == 0L) {
    return(
      list(
        any_converged = FALSE,
        all_converged = FALSE,
        best_restart = NA_integer_,
        best_logLik = NA_real_,
        best_iterations = NA_real_,
        n_restarts = 0L,
        n_converged = 0L
      )
    )
  }

  nm <- names(em_dt)

  payload_cols <- setdiff(
    nm,
    "restart"
  )

  ll_col <- grep(
    "^log.?lik|likelihood|^ll$",
    nm,
    ignore.case = TRUE,
    value = TRUE
  )

  iter_col <- grep(
    "iter",
    nm,
    ignore.case = TRUE,
    value = TRUE
  )

  conv_col <- grep(
    "converg",
    nm,
    ignore.case = TRUE,
    value = TRUE
  )

  # gbmt 0.1.4 special case:
  # G1 can return the EM record without semantic column names,
  # which becomes V1/V2/V3 after as.data.frame().
  # Positional fallback:
  #   first payload column  = logLik
  #   second payload column = n.iter
  #   third payload column  = converged
  if (
    length(ll_col) == 0L &&
    length(payload_cols) >= 1L
  ) {
    ll_col <- payload_cols[1L]
  }

  if (
    length(iter_col) == 0L &&
    length(payload_cols) >= 2L
  ) {
    iter_col <- payload_cols[2L]
  }

  if (
    length(conv_col) == 0L &&
    length(payload_cols) >= 3L
  ) {
    conv_col <- payload_cols[3L]
  }

  if (
    length(ll_col) == 0L ||
    length(iter_col) == 0L ||
    length(conv_col) == 0L
  ) {
    stop(
      "Unable to identify logLik, n.iter, and converged columns in fit$em."
    )
  }

  ll_vec <- suppressWarnings(
    as.numeric(
      em_dt[[ll_col[1L]]]
    )
  )

  iter_vec <- suppressWarnings(
    as.numeric(
      em_dt[[iter_col[1L]]]
    )
  )

  conv_vec <- suppressWarnings(
    as.numeric(
      em_dt[[conv_col[1L]]]
    )
  )

  if (
    any(
      is.finite(ll_vec)
    )
  ) {

    best_restart <- which.max(
      ifelse(
        is.finite(ll_vec),
        ll_vec,
        -Inf
      )
    )

    best_ll <- ll_vec[
      best_restart
    ]

    best_iterations <- iter_vec[
      best_restart
    ]

  } else {

    best_restart <- NA_integer_
    best_ll <- NA_real_
    best_iterations <- NA_real_

  }

  finite_conv <- conv_vec[
    is.finite(
      conv_vec
    )
  ]

  list(
    any_converged = any(
      finite_conv == 1
    ),
    all_converged = length(finite_conv) > 0L &&
      all(
        finite_conv == 1
      ),
    best_restart = as.integer(
      best_restart
    ),
    best_logLik = as.numeric(
      best_ll
    ),
    best_iterations = as.numeric(
      best_iterations
    ),
    n_restarts = as.integer(
      nrow(em_dt)
    ),
    n_converged = as.integer(
      sum(
        finite_conv == 1
      )
    )
  )

}


extract_class_quality <- function(
  fit,
  requested_g,
  indicator
) {

  post <- as.matrix(
    fit$posterior
  )

  assign_vec <- as.integer(
    fit$assign
  )

  if (nrow(post) != length(assign_vec)) {
    stop(
      "Posterior row count differs from assignment length."
    )
  }

  posterior_ncol <- ncol(
    post
  )

  stored_ng <- suppressWarnings(
    as.integer(
      fit$call$ng
    )
  )

  if (
    length(stored_ng) == 0L ||
    is.na(stored_ng)
  ) {
    stored_ng <- posterior_ncol
  }

  if (posterior_ncol < 1L) {
    stop(
      "Returned posterior matrix has zero columns."
    )
  }

  modal_post <- max.col(
    post,
    ties.method = "first"
  )

  mismatch_n <- sum(
    assign_vec != modal_post,
    na.rm = TRUE
  )

  class_dt <- rbindlist(
    lapply(
      seq_len(
        posterior_ncol
      ),
      function(k) {

        idx <- which(
          assign_vec == k
        )

        n_k <- length(
          idx
        )

        appa <- if (n_k > 0L) {

          mean(
            post[
              idx,
              k
            ],
            na.rm = TRUE
          )

        } else {

          NA_real_

        }

        soft_mass <- sum(
          post[
            ,
            k
          ],
          na.rm = TRUE
        )

        data.table(
          indicator = indicator,
          requested_g = as.integer(
            requested_g
          ),
          stored_ng = as.integer(
            stored_ng
          ),
          class = as.integer(k),
          n = as.integer(n_k),
          pct = as.numeric(
            100 *
              n_k /
              length(assign_vec)
          ),
          APPA = as.numeric(
            appa
          ),
          posterior_soft_mass = as.numeric(
            soft_mass
          ),
          posterior_soft_mass_pct = as.numeric(
            100 *
              soft_mass /
              nrow(post)
          ),
          max_posterior_for_class = as.numeric(
            max(
              post[
                ,
                k
              ],
              na.rm = TRUE
            )
          )
        )

      }
    )
  )

  occupied_groups <- sum(
    class_dt$n > 0L
  )

  empty_classes <- class_dt[
    n == 0L,
    class
  ]

  nonempty_dt <- class_dt[
    n > 0L
  ]

  min_class_pct <- if (
    nrow(nonempty_dt) > 0L
  ) {

    min(
      nonempty_dt$pct,
      na.rm = TRUE
    )

  } else {

    NA_real_

  }

  min_appa <- if (
    nrow(nonempty_dt) > 0L &&
    any(
      is.finite(
        nonempty_dt$APPA
      )
    )
  ) {

    min(
      nonempty_dt$APPA[
        is.finite(
          nonempty_dt$APPA
        )
      ]
    )

  } else {

    NA_real_

  }

  max_post <- apply(
    post,
    1L,
    max,
    na.rm = TRUE
  )

  structure_exact <- (
    stored_ng ==
      requested_g
  ) &&
    (
      posterior_ncol ==
        requested_g
    )

  all_occupied <- (
    occupied_groups ==
      requested_g
  )

  full_g <- structure_exact &&
    all_occupied

  summary_dt <- data.table(
    indicator = indicator,
    requested_g = as.integer(
      requested_g
    ),
    stored_ng = as.integer(
      stored_ng
    ),
    posterior_ncol = as.integer(
      posterior_ncol
    ),
    occupied_groups = as.integer(
      occupied_groups
    ),
    structure_exact = as.integer(
      structure_exact
    ),
    all_classes_occupied = as.integer(
      all_occupied
    ),
    full_g = as.integer(
      full_g
    ),
    empty_classes = if (
      length(
        empty_classes
      ) > 0L
    ) {

      paste(
        empty_classes,
        collapse = ","
      )

    } else {

      ""

    },
    assignment_posterior_mismatch_n = as.integer(
      mismatch_n
    ),
    min_class_pct_nonempty = as.numeric(
      min_class_pct
    ),
    min_APPA_nonempty = as.numeric(
      min_appa
    ),
    mean_max_posterior = as.numeric(
      mean(
        max_post,
        na.rm = TRUE
      )
    ),
    median_max_posterior = as.numeric(
      median(
        max_post,
        na.rm = TRUE
      )
    ),
    pct_max_posterior_ge_0_70 = as.numeric(
      100 *
        mean(
          max_post >= 0.70,
          na.rm = TRUE
        )
    ),
    pct_max_posterior_ge_0_90 = as.numeric(
      100 *
        mean(
          max_post >= 0.90,
          na.rm = TRUE
        )
    )
  )

  list(
    class_table = class_dt,
    summary = summary_dt
  )

}


extract_trajectory_table <- function(
  fit,
  indicator,
  requested_g
) {

  if (!inherits(fit, "gbmt")) {
    stop(
      paste0(
        indicator,
        " G",
        requested_g,
        ": fitted object is not class 'gbmt'."
      )
    )
  }

  fit_x_names <- as.character(
    fit$call$x.names
  )

  if (length(fit_x_names) != 1L) {
    stop(
      paste0(
        indicator,
        " G",
        requested_g,
        ": Stage-1 model must contain exactly one indicator."
      )
    )
  }

  variable_name <- fit_x_names[1L]

  time_col <- as.character(
    fit$call$time
  )

  if (
    length(time_col) != 1L ||
    !time_col %in%
      names(fit$data.orig)
  ) {
    stop(
      paste0(
        indicator,
        " G",
        requested_g,
        ": fitted-object time column cannot be resolved."
      )
    )
  }

  original_times <- sort(
    unique(
      suppressWarnings(
        as.numeric(
          fit$data.orig[[time_col]]
        )
      )
    )
  )

  if (
    length(original_times) == 0L ||
    anyNA(original_times)
  ) {
    stop(
      paste0(
        indicator,
        " G",
        requested_g,
        ": invalid original time values."
      )
    )
  }

  # gbmt 0.1.4 stores final estimated group trajectories directly:
  #
  #   fit$fitted[[group]]
  #
  # These are constructed internally as Xmat %*% beta with original
  # time-level names assigned as row names. Using this stored object
  # avoids predict.gbmt() and preserves the exact fitted trajectories.
  fitted_list <- fit$fitted

  if (
    !is.list(fitted_list) ||
    length(fitted_list) == 0L
  ) {
    stop(
      paste0(
        indicator,
        " G",
        requested_g,
        ": fit$fitted is missing or empty."
      )
    )
  }

  stored_g <- length(
    fitted_list
  )

  if (stored_g != requested_g) {
    stop(
      paste0(
        indicator,
        " G",
        requested_g,
        ": fit$fitted contains ",
        stored_g,
        " group(s), not requested G=",
        requested_g,
        "."
      )
    )
  }

  out_list <- vector(
    "list",
    stored_g
  )

  for (
    k in seq_len(
      stored_g
    )
  ) {

    group_fit <- as.data.frame(
      fitted_list[[k]]
    )

    if (nrow(group_fit) == 0L) {
      stop(
        paste0(
          indicator,
          " G",
          requested_g,
          " class ",
          k,
          ": fit$fitted group dataframe is empty."
        )
      )
    }

    if (
      !variable_name %in%
        names(group_fit)
    ) {
      stop(
        paste0(
          indicator,
          " G",
          requested_g,
          " class ",
          k,
          ": indicator '",
          variable_name,
          "' is absent from fit$fitted."
        )
      )
    }

    predicted_value <- suppressWarnings(
      as.numeric(
        group_fit[[variable_name]]
      )
    )

    fitted_times <- suppressWarnings(
      as.numeric(
        rownames(
          group_fit
        )
      )
    )

    if (
      length(fitted_times) != nrow(group_fit) ||
      anyNA(fitted_times)
    ) {

      xmat_times <- suppressWarnings(
        as.numeric(
          rownames(
            fit$Xmat
          )
        )
      )

      if (
        length(xmat_times) ==
          nrow(group_fit) &&
        !anyNA(xmat_times)
      ) {

        fitted_times <- xmat_times

      } else {

        stop(
          paste0(
            indicator,
            " G",
            requested_g,
            " class ",
            k,
            ": time labels cannot be recovered from fit$fitted or fit$Xmat."
          )
        )

      }

    }

    if (
      length(fitted_times) !=
        length(original_times)
    ) {
      stop(
        paste0(
          indicator,
          " G",
          requested_g,
          " class ",
          k,
          ": fitted time count differs from original unique time count."
        )
      )
    }

    if (
      !isTRUE(
        all.equal(
          sort(fitted_times),
          original_times,
          tolerance = 1e-8,
          check.attributes = FALSE
        )
      )
    ) {
      stop(
        paste0(
          indicator,
          " G",
          requested_g,
          " class ",
          k,
          ": fit$fitted times do not match original fitted-model time values."
        )
      )
    }

    if (
      length(predicted_value) !=
        length(fitted_times) ||
      anyNA(predicted_value) ||
      any(
        !is.finite(
          predicted_value
        )
      )
    ) {
      stop(
        paste0(
          indicator,
          " G",
          requested_g,
          " class ",
          k,
          ": invalid fitted trajectory values."
        )
      )
    }

    time_hour_center <- 12 +
      12 *
      fitted_times

    out_list[[k]] <- data.table(
      indicator = indicator,
      variable = variable_name,
      requested_g = as.integer(
        requested_g
      ),
      stored_g = as.integer(
        stored_g
      ),
      class = as.integer(
        k
      ),
      time_model = as.numeric(
        fitted_times
      ),
      time_hour_center = as.numeric(
        time_hour_center
      ),
      predicted_z = as.numeric(
        predicted_value
      )
    )

  }

  out <- rbindlist(
    out_list,
    fill = TRUE
  )

  setorder(
    out,
    class,
    time_model
  )

  expected_rows <- requested_g *
    length(
      original_times
    )

  if (
    nrow(out) !=
      expected_rows
  ) {
    stop(
      paste0(
        indicator,
        " G",
        requested_g,
        ": trajectory output row count = ",
        nrow(out),
        "; expected ",
        expected_rows,
        "."
      )
    )
  }

  out

}


task_is_complete <- function(
  indicator,
  g,
  output_dir
) {

  summary_path <- file.path(
    output_dir,
    paste0(
      "V267_task_summary_",
      indicator,
      "_G",
      g,
      ".csv"
    )
  )

  fit_path <- file.path(
    output_dir,
    paste0(
      "V267_screen_fit_",
      indicator,
      "_G",
      g,
      ".rds"
    )
  )

  if (
    !file.exists(summary_path) ||
    !file.exists(fit_path)
  ) {
    return(FALSE)
  }

  d <- tryCatch(
    fread(
      summary_path
    ),
    error = function(e) NULL
  )

  if (
    is.null(d) ||
    nrow(d) < 1L ||
    !"status" %in%
      names(d)
  ) {
    return(FALSE)
  }

  completed_statuses <- c(
    "CONVERGED_VALID_FULL_G",
    "CONVERGED_EMPTY_CLASS",
    "CONVERGED_STRUCTURE_MISMATCH",
    "NONCONVERGED_FULL_G",
    "NONCONVERGED_EMPTY_CLASS",
    "NONCONVERGED_STRUCTURE_MISMATCH"
  )

  as.character(
    d$status[1L]
  ) %in%
    completed_statuses

}


# ============================================================
# 10. One-model worker
# ============================================================

run_one_stage1_model <- function(
  task_row,
  output_dir,
  version_tag,
  model_degree,
  screen_tol,
  nstart_g1,
  nstart_g2plus,
  stagger_seconds,
  n_workers,
  nonbinding_maxit,
  master_seed,
  expected_patients,
  preferred_min_class_pct,
  preferred_min_appa
) {

  task_id <- as.integer(
    task_row$task_id[1L]
  )

  indicator <- as.character(
    task_row$indicator[1L]
  )

  variable <- as.character(
    task_row$variable[1L]
  )

  g <- as.integer(
    task_row$requested_g[1L]
  )

  indicator_input_file <- as.character(
    task_row$input_file[1L]
  )

  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )

  suppressPackageStartupMessages({
    library(data.table)
    library(gbmt)
  })

  data.table::setDTthreads(
    1L
  )

  fit_path <- file.path(
    output_dir,
    paste0(
      "V267_screen_fit_",
      indicator,
      "_G",
      g,
      ".rds"
    )
  )

  partial_fit_path <- file.path(
    output_dir,
    paste0(
      "V267_screen_fit_partial_",
      indicator,
      "_G",
      g,
      ".rds"
    )
  )

  summary_path <- file.path(
    output_dir,
    paste0(
      "V267_task_summary_",
      indicator,
      "_G",
      g,
      ".csv"
    )
  )

  class_path <- file.path(
    output_dir,
    paste0(
      "V267_class_quality_",
      indicator,
      "_G",
      g,
      ".csv"
    )
  )

  em_runs_path <- file.path(
    output_dir,
    paste0(
      "V267_EM_runs_",
      indicator,
      "_G",
      g,
      ".csv"
    )
  )

  trajectory_path <- file.path(
    output_dir,
    paste0(
      "V267_trajectory_",
      indicator,
      "_G",
      g,
      ".csv"
    )
  )

  em_progress_path <- file.path(
    output_dir,
    paste0(
      "V267_EM_progress_",
      indicator,
      "_G",
      g,
      ".log"
    )
  )

  stage_path <- file.path(
    output_dir,
    paste0(
      "V267_stage_",
      indicator,
      "_G",
      g,
      ".log"
    )
  )

  running_path <- file.path(
    output_dir,
    paste0(
      "V267_RUNNING_",
      indicator,
      "_G",
      g,
      ".txt"
    )
  )

  done_path <- file.path(
    output_dir,
    paste0(
      "V267_DONE_",
      indicator,
      "_G",
      g,
      ".txt"
    )
  )

  error_path <- file.path(
    output_dir,
    paste0(
      "V267_ERROR_",
      indicator,
      "_G",
      g,
      ".txt"
    )
  )

  completed_statuses <- c(
    "CONVERGED_VALID_FULL_G",
    "CONVERGED_EMPTY_CLASS",
    "CONVERGED_STRUCTURE_MISMATCH",
    "NONCONVERGED_FULL_G",
    "NONCONVERGED_EMPTY_CLASS",
    "NONCONVERGED_STRUCTURE_MISMATCH"
  )

  if (
    file.exists(summary_path) &&
    file.exists(fit_path)
  ) {

    old <- tryCatch(
      fread(
        summary_path
      ),
      error = function(e) NULL
    )

    if (
      !is.null(old) &&
      nrow(old) >= 1L &&
      "status" %in%
        names(old) &&
      as.character(
        old$status[1L]
      ) %in%
        completed_statuses
    ) {

      old[
        ,
        resumed_skip := 1L
      ]

      return(
        old
      )

    }

  }

  cleanup_paths <- c(
    partial_fit_path,
    running_path,
    done_path,
    error_path
  )

  for (f in cleanup_paths) {

    if (
      file.exists(f)
    ) {
      file.remove(f)
    }

  }

  if (file.exists(stage_path)) {
    file.remove(
      stage_path
    )
  }

  if (file.exists(em_progress_path)) {
    file.remove(
      em_progress_path
    )
  }

  nstart_used <- if (
    g == 1L
  ) {
    as.integer(
      nstart_g1
    )
  } else {
    as.integer(
      nstart_g2plus
    )
  }

  seed_used <- as.integer(
    master_seed +
      task_id * 100003L +
      g * 1009L
  )

  # Each load-balanced batch is staggered only within the worker count,
  # avoiding an ever-growing sleep for later tasks.
  stagger_slot <- (
    task_id - 1L
  ) %%
    n_workers

  sleep_seconds <- as.integer(
    stagger_slot *
      stagger_seconds
  )

  start_time <- Sys.time()

  writeLines(
    c(
      paste0(
        "Version: ",
        version_tag
      ),
      "Stage: 1 single-indicator GBTM",
      paste0(
        "Indicator: ",
        indicator
      ),
      paste0(
        "Variable: ",
        variable
      ),
      paste0(
        "Requested G: ",
        g
      ),
      paste0(
        "nstart: ",
        nstart_used
      ),
      paste0(
        "tol: ",
        screen_tol
      ),
      "delete.empty: FALSE",
      "Automatic rerun: NONE",
      paste0(
        "Technical maxit: ",
        nonbinding_maxit
      ),
      paste0(
        "Seed: ",
        seed_used
      ),
      paste0(
        "PID: ",
        Sys.getpid()
      ),
      paste0(
        "Start: ",
        format(
          start_time,
          "%Y-%m-%d %H:%M:%S"
        )
      )
    ),
    running_path
  )

  append_stage(
    stage_path,
    "worker entered"
  )

  if (
    sleep_seconds > 0L
  ) {

    append_stage(
      stage_path,
      paste0(
        "stagger sleep ",
        sleep_seconds,
        " sec"
      )
    )

    Sys.sleep(
      sleep_seconds
    )

  }

  append_stage(
    stage_path,
    "before readRDS indicator input"
  )

  model_data <- readRDS(
    indicator_input_file
  )

  append_stage(
    stage_path,
    "after readRDS indicator input"
  )

  if (
    !is.data.frame(
      model_data
    )
  ) {

    model_data <- as.data.frame(
      model_data
    )

  }

  model_data$unit_id <- as.character(
    model_data$unit_id
  )

  model_data$time_model <- as.numeric(
    model_data$time_model
  )

  model_data[[variable]] <- as.numeric(
    model_data[[variable]]
  )

  if (
    anyNA(
      model_data[[variable]]
    )
  ) {
    stop(
      paste0(
        indicator,
        " worker input unexpectedly contains NA values."
      )
    )
  }

  n_units <- length(
    unique(
      model_data$unit_id
    )
  )

  if (
    n_units !=
      expected_patients
  ) {
    stop(
      paste0(
        indicator,
        " G",
        g,
        " worker patient count = ",
        n_units,
        "; expected ",
        expected_patients,
        "."
      )
    )
  }

  append_stage(
    stage_path,
    "worker input validation passed"
  )

  gc(
    full = TRUE
  )

  set.seed(
    seed_used
  )

  cat(
    paste0(
      version_tag,
      " ",
      indicator,
      " G",
      g,
      " fitting started at ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      "\n"
    ),
    file = em_progress_path,
    append = TRUE
  )

  append_stage(
    stage_path,
    "before gbmt call"
  )

  fit <- NULL
  fit_error <- NA_character_

  fit_start <- Sys.time()

  sink(
    em_progress_path,
    append = TRUE,
    type = "output"
  )

  fit <- tryCatch(
    {

      gbmt::gbmt(
        x.names = variable,
        unit = "unit_id",
        time = "time_model",
        ng = as.integer(g),
        d = as.integer(
          model_degree
        ),
        data = model_data,
        scaling = 0,
        pruning = TRUE,
        nstart = as.integer(
          nstart_used
        ),
        tol = as.numeric(
          screen_tol
        ),
        maxit = as.integer(
          nonbinding_maxit
        ),
        delete.empty = FALSE,
        quiet = FALSE
      )

    },
    error = function(e) {

      fit_error <<- conditionMessage(
        e
      )

      NULL

    },
    finally = {

      while (
        sink.number(
          type = "output"
        ) > 0L
      ) {

        sink(
          type = "output"
        )

      }

    }
  )

  fit_end <- Sys.time()

  fit_seconds <- as.numeric(
    difftime(
      fit_end,
      fit_start,
      units = "secs"
    )
  )

  if (is.null(fit)) {

    append_stage(
      stage_path,
      paste0(
        "gbmt ERROR: ",
        fit_error
      )
    )

    error_dt <- data.table(
      version = version_tag,
      stage = 1L,
      indicator = indicator,
      variable = variable,
      requested_g = as.integer(g),
      status = "ERROR",
      error_message = as.character(
        fit_error
      ),
      nstart_used = as.integer(
        nstart_used
      ),
      screening_tol = as.numeric(
        screen_tol
      ),
      seed = as.integer(
        seed_used
      ),
      fit_seconds = as.numeric(
        fit_seconds
      ),
      start_time = format(
        start_time,
        "%Y-%m-%d %H:%M:%S"
      ),
      end_time = format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    )

    fwrite(
      error_dt,
      summary_path
    )

    if (
      file.exists(
        running_path
      )
    ) {
      file.remove(
        running_path
      )
    }

    writeLines(
      c(
        paste0(
          "Version: ",
          version_tag
        ),
        paste0(
          "Indicator: ",
          indicator
        ),
        paste0(
          "G: ",
          g
        ),
        "Status: ERROR",
        paste0(
          "Error: ",
          fit_error
        )
      ),
      error_path
    )

    return(
      error_dt
    )

  }

  append_stage(
    stage_path,
    "gbmt call returned"
  )

  # ----------------------------------------------------------
  # Preserve exact fit object before diagnostics
  # ----------------------------------------------------------

  saveRDS(
    fit,
    partial_fit_path
  )

  if (
    !file.exists(
      partial_fit_path
    ) ||
    is.na(
      file.info(
        partial_fit_path
      )$size
    ) ||
    file.info(
      partial_fit_path
    )$size <=
      0
  ) {
    stop(
      paste0(
        indicator,
        " G",
        g,
        " partial model RDS was not written correctly."
      )
    )
  }

  if (
    file.exists(
      fit_path
    )
  ) {
    file.remove(
      fit_path
    )
  }

  renamed <- file.rename(
    partial_fit_path,
    fit_path
  )

  if (!renamed) {

    copied <- file.copy(
      partial_fit_path,
      fit_path,
      overwrite = TRUE
    )

    if (!copied) {
      stop(
        paste0(
          "Could not finalize model RDS for ",
          indicator,
          " G",
          g,
          "."
        )
      )
    }

    file.remove(
      partial_fit_path
    )

  }

  # ----------------------------------------------------------
  # Diagnostics
  # ----------------------------------------------------------

  em_dt <- extract_em_table(
    fit
  )

  fwrite(
    em_dt,
    em_runs_path
  )

  em_status <- extract_em_status(
    fit
  )

  class_obj <- extract_class_quality(
    fit = fit,
    requested_g = g,
    indicator = indicator
  )

  fwrite(
    class_obj$class_table,
    class_path
  )

  trajectory_dt <- extract_trajectory_table(
    fit = fit,
    indicator = indicator,
    requested_g = g
  )

  fwrite(
    trajectory_dt,
    trajectory_path
  )

  q <- class_obj$summary

  package_converged <- isTRUE(
    em_status$any_converged
  )

  structure_exact <- (
    q$structure_exact[1L] ==
      1L
  )

  all_occupied <- (
    q$all_classes_occupied[1L] ==
      1L
  )

  full_g <- (
    q$full_g[1L] ==
      1L
  )

  status <- if (
    package_converged &&
    full_g
  ) {

    "CONVERGED_VALID_FULL_G"

  } else if (
    package_converged &&
    structure_exact &&
    !all_occupied
  ) {

    "CONVERGED_EMPTY_CLASS"

  } else if (
    package_converged &&
    !structure_exact
  ) {

    "CONVERGED_STRUCTURE_MISMATCH"

  } else if (
    !package_converged &&
    full_g
  ) {

    "NONCONVERGED_FULL_G"

  } else if (
    !package_converged &&
    structure_exact &&
    !all_occupied
  ) {

    "NONCONVERGED_EMPTY_CLASS"

  } else {

    "NONCONVERGED_STRUCTURE_MISMATCH"

  }

  final_logLik <- suppressWarnings(
    as.numeric(
      logLik(
        fit
      )
    )
  )

  AIC_value <- safe_ic_value(
    fit,
    "AIC"
  )

  BIC_value <- safe_ic_value(
    fit,
    "BIC"
  )

  CAIC_value <- safe_ic_value(
    fit,
    "CAIC"
  )

  ssBIC_value <- safe_ic_value(
    fit,
    "ssBIC"
  )

  HQIC_value <- safe_ic_value(
    fit,
    "HQIC"
  )

  if (
    !is.finite(
      AIC_value
    )
  ) {

    AIC_value <- tryCatch(
      as.numeric(
        AIC(
          fit
        )
      ),
      error = function(e) NA_real_
    )

  }

  if (
    !is.finite(
      BIC_value
    )
  ) {

    BIC_value <- tryCatch(
      as.numeric(
        BIC(
          fit
        )
      ),
      error = function(e) NA_real_
    )

  }

  min_class_pct <- q$min_class_pct_nonempty[1L]
  min_appa <- q$min_APPA_nonempty[1L]

  preferred_class_size_ok <- is.finite(
    min_class_pct
  ) &&
    min_class_pct >=
      preferred_min_class_pct

  preferred_appa_ok <- is.finite(
    min_appa
  ) &&
    min_appa >=
      preferred_min_appa

  ic_comparison_eligible <- (
    status ==
      "CONVERGED_VALID_FULL_G"
  )

  end_time <- Sys.time()

  summary_dt <- data.table(
    version = version_tag,
    stage = 1L,
    indicator = indicator,
    variable = variable,
    requested_g = as.integer(
      g
    ),
    status = as.character(
      status
    ),
    package_any_restart_converged = as.integer(
      em_status$any_converged
    ),
    package_all_restarts_converged = as.integer(
      em_status$all_converged
    ),
    n_restarts = as.integer(
      em_status$n_restarts
    ),
    n_converged_restarts = as.integer(
      em_status$n_converged
    ),
    best_restart = as.integer(
      em_status$best_restart
    ),
    best_restart_logLik = as.numeric(
      em_status$best_logLik
    ),
    best_restart_iterations = as.numeric(
      em_status$best_iterations
    ),
    stored_ng = as.integer(
      q$stored_ng[1L]
    ),
    posterior_ncol = as.integer(
      q$posterior_ncol[1L]
    ),
    occupied_groups = as.integer(
      q$occupied_groups[1L]
    ),
    structure_exact = as.integer(
      structure_exact
    ),
    all_classes_occupied = as.integer(
      all_occupied
    ),
    full_g = as.integer(
      full_g
    ),
    empty_classes = as.character(
      q$empty_classes[1L]
    ),
    screening_tol = as.numeric(
      screen_tol
    ),
    nstart_used = as.integer(
      nstart_used
    ),
    delete_empty = 0L,
    technical_maxit = as.integer(
      nonbinding_maxit
    ),
    automatic_rerun = 0L,
    final_logLik = as.numeric(
      final_logLik
    ),
    AIC = as.numeric(
      AIC_value
    ),
    BIC = as.numeric(
      BIC_value
    ),
    CAIC = as.numeric(
      CAIC_value
    ),
    ssBIC = as.numeric(
      ssBIC_value
    ),
    HQIC = as.numeric(
      HQIC_value
    ),
    min_class_pct_nonempty = as.numeric(
      min_class_pct
    ),
    min_APPA_nonempty = as.numeric(
      min_appa
    ),
    mean_max_posterior = as.numeric(
      q$mean_max_posterior[1L]
    ),
    median_max_posterior = as.numeric(
      q$median_max_posterior[1L]
    ),
    pct_max_posterior_ge_0_70 = as.numeric(
      q$pct_max_posterior_ge_0_70[1L]
    ),
    pct_max_posterior_ge_0_90 = as.numeric(
      q$pct_max_posterior_ge_0_90[1L]
    ),
    assignment_posterior_mismatch_n = as.integer(
      q$assignment_posterior_mismatch_n[1L]
    ),
    preferred_class_size_ok = as.integer(
      preferred_class_size_ok
    ),
    preferred_APPA_ok = as.integer(
      preferred_appa_ok
    ),
    IC_comparison_eligible = as.integer(
      ic_comparison_eligible
    ),
    seed = as.integer(
      seed_used
    ),
    fit_seconds = as.numeric(
      fit_seconds
    ),
    total_task_seconds = as.numeric(
      difftime(
        end_time,
        start_time,
        units = "secs"
      )
    ),
    n_units = as.integer(
      n_units
    ),
    n_observed_rows = as.integer(
      nrow(
        model_data
      )
    ),
    pid = as.integer(
      Sys.getpid()
    ),
    start_time = format(
      start_time,
      "%Y-%m-%d %H:%M:%S"
    ),
    end_time = format(
      end_time,
      "%Y-%m-%d %H:%M:%S"
    ),
    model_file = normalizePath(
      fit_path,
      winslash = "/",
      mustWork = TRUE
    ),
    resumed_skip = 0L,
    error_message = NA_character_
  )

  fwrite(
    summary_dt,
    summary_path
  )

  if (
    file.exists(
      running_path
    )
  ) {
    file.remove(
      running_path
    )
  }

  writeLines(
    c(
      paste0(
        "Version: ",
        version_tag
      ),
      "Stage: 1",
      paste0(
        "Indicator: ",
        indicator
      ),
      paste0(
        "Requested G: ",
        g
      ),
      paste0(
        "Status: ",
        status
      ),
      paste0(
        "nstart used: ",
        nstart_used
      ),
      paste0(
        "Converged restarts: ",
        em_status$n_converged,
        "/",
        em_status$n_restarts
      ),
      paste0(
        "Best restart: ",
        em_status$best_restart
      ),
      paste0(
        "Requested G: ",
        g
      ),
      paste0(
        "Stored ng: ",
        q$stored_ng[1L]
      ),
      paste0(
        "Occupied groups: ",
        q$occupied_groups[1L]
      ),
      paste0(
        "Empty classes: ",
        ifelse(
          nzchar(
            q$empty_classes[1L]
          ),
          q$empty_classes[1L],
          "none"
        )
      ),
      paste0(
        "Final LL: ",
        final_logLik
      ),
      paste0(
        "BIC: ",
        BIC_value
      ),
      paste0(
        "Min class %: ",
        min_class_pct
      ),
      paste0(
        "Min APPA: ",
        min_appa
      ),
      paste0(
        "IC comparison eligible: ",
        ic_comparison_eligible
      ),
      paste0(
        "End: ",
        format(
          end_time,
          "%Y-%m-%d %H:%M:%S"
        )
      )
    ),
    done_path
  )

  append_stage(
    stage_path,
    paste0(
      "task finished: ",
      status
    )
  )

  summary_dt

}



# ============================================================
# 11. Build task table
# ============================================================

task_dt <- CJ(
  indicator = INDICATOR_MAP$indicator,
  requested_g = GROUPS_TO_SCREEN,
  unique = TRUE
)

task_dt <- merge(
  task_dt,
  as.data.table(
    INDICATOR_MAP
  )[
    ,
    .(
      indicator,
      variable
    )
  ],
  by = "indicator",
  all.x = TRUE,
  sort = FALSE
)

task_dt <- merge(
  task_dt,
  indicator_input_dt[
    ,
    .(
      indicator,
      input_file
    )
  ],
  by = "indicator",
  all.x = TRUE,
  sort = FALSE
)

# Longest-first scheduling:
# higher G first improves total makespan under load balancing.
setorder(
  task_dt,
  -requested_g,
  indicator
)

task_dt[
  ,
  task_id := seq_len(.N)
]

setcolorder(
  task_dt,
  c(
    "task_id",
    "indicator",
    "variable",
    "requested_g",
    "input_file"
  )
)

fwrite(
  task_dt,
  file.path(
    output_dir,
    "V267_task_plan.csv"
  )
)


# ============================================================
# 12. Determine pending tasks
# ============================================================

task_dt[
  ,
  complete := mapply(
    FUN = task_is_complete,
    indicator = indicator,
    g = requested_g,
    MoreArgs = list(
      output_dir = output_dir
    )
  )
]

pending_dt <- task_dt[
  complete == FALSE
]

cat(
  "Total Stage-1 models: ",
  nrow(task_dt),
  "\n",
  sep = ""
)

cat(
  "Pending models: ",
  nrow(pending_dt),
  "\n\n",
  sep = ""
)


# ============================================================
# 13. Master markers
# ============================================================

master_start <- Sys.time()

master_running_path <- file.path(
  output_dir,
  "V267_MASTER_RUNNING.txt"
)

master_done_path <- file.path(
  output_dir,
  "V267_MASTER_DONE.txt"
)

if (
  file.exists(
    master_done_path
  )
) {
  file.remove(
    master_done_path
  )
}

writeLines(
  c(
    "V267 Stage-1 single-indicator GBTM master run is active.",
    paste0(
      "Start: ",
      format(
        master_start,
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    "Indicators: MSI, Na, K, Cl, HCO3",
    "Each indicator modeled separately.",
    "G range: G1-G6",
    "Workers: 8",
    "tol: 1",
    "G1 nstart: 1",
    "G2-G6 nstart: 3",
    "delete.empty: FALSE",
    "Automatic rerun: NONE",
    "Outcome information: NONE"
  ),
  master_running_path
)


# ============================================================
# 14. Run pending Stage-1 models
# ============================================================

if (
  nrow(
    pending_dt
  ) > 0L
) {

  workers_now <- min(
    N_WORKERS,
    nrow(
      pending_dt
    )
  )

  cat(
    "Starting ",
    workers_now,
    " PSOCK workers...\n\n",
    sep = ""
  )

  worker_console_path <- file.path(
    output_dir,
    "V267_PSOCK_worker_console.log"
  )

  cl <- parallel::makePSOCKcluster(
    workers_now,
    outfile = worker_console_path
  )

  cluster_stopped <- FALSE

  on.exit(
    {

      if (
        exists("cl") &&
        !is.null(cl) &&
        !cluster_stopped
      ) {

        try(
          parallel::stopCluster(
            cl
          ),
          silent = TRUE
        )

      }

    },
    add = TRUE
  )

  # ----------------------------------------------------------
  # PSOCK preflight
  # ----------------------------------------------------------
  # Before launching any GBMTM fit, verify that EVERY worker:
  #   1) receives the explicit output_dir argument,
  #   2) can see the output directory,
  #   3) can load data.table and gbmt,
  #   4) is running with nested numerical threading disabled.
  #
  # This prevents the old "object 'output_dir' not found" failure
  # from reaching the expensive model-fitting stage.
  # ----------------------------------------------------------

  preflight <- parallel::clusterCall(
    cl,
    function(test_output_dir) {

      Sys.setenv(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        BLAS_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        NUMEXPR_NUM_THREADS = "1"
      )

      dt_ok <- requireNamespace(
        "data.table",
        quietly = TRUE
      )

      gbmt_ok <- requireNamespace(
        "gbmt",
        quietly = TRUE
      )

      list(
        pid = Sys.getpid(),
        output_dir_received = is.character(test_output_dir) &&
          length(test_output_dir) == 1L &&
          nzchar(test_output_dir),
        output_dir_exists = dir.exists(test_output_dir),
        data_table_ok = dt_ok,
        gbmt_ok = gbmt_ok,
        gbmt_version = if (gbmt_ok) {
          as.character(
            packageVersion("gbmt")
          )
        } else {
          NA_character_
        }
      )

    },
    output_dir
  )

  preflight_dt <- data.table::rbindlist(
    lapply(
      seq_along(preflight),
      function(i) {

        x <- preflight[[i]]

        data.table::data.table(
          worker = as.integer(i),
          pid = as.integer(x$pid),
          output_dir_received = as.logical(
            x$output_dir_received
          ),
          output_dir_exists = as.logical(
            x$output_dir_exists
          ),
          data_table_ok = as.logical(
            x$data_table_ok
          ),
          gbmt_ok = as.logical(
            x$gbmt_ok
          ),
          gbmt_version = as.character(
            x$gbmt_version
          )
        )

      }
    ),
    fill = TRUE
  )

  data.table::fwrite(
    preflight_dt,
    file.path(
      output_dir,
      "V267_PSOCK_preflight.csv"
    )
  )

  preflight_ok <- all(
    preflight_dt$output_dir_received &
      preflight_dt$output_dir_exists &
      preflight_dt$data_table_ok &
      preflight_dt$gbmt_ok
  )

  if (!preflight_ok) {

    parallel::stopCluster(
      cl
    )

    cluster_stopped <- TRUE
    cl <- NULL

    stop(
      paste0(
        "V267 PSOCK preflight failed. ",
        "No GBMTM model was started.\n",
        "Check:\n",
        file.path(
          output_dir,
          "V267_PSOCK_preflight.csv"
        )
      )
    )

  }

  cat(
    "PSOCK preflight passed on all ",
    workers_now,
    " workers.\n\n",
    sep = ""
  )

  # ----------------------------------------------------------
  # Export ONLY functions required by the real worker.
  # No configuration object is required from worker globals.
  # ----------------------------------------------------------

  parallel::clusterExport(
    cl,
    varlist = c(
      "run_one_stage1_model",
      "append_stage",
      "safe_ic_value",
      "extract_em_table",
      "extract_em_status",
      "extract_class_quality",
      "extract_trajectory_table"
    ),
    envir = environment()
  )

  # ----------------------------------------------------------
  # Build self-contained task records.
  # Each list element contains task_id, indicator, variable,
  # requested_g, and the exact indicator input-file path.
  # ----------------------------------------------------------

  task_records <- lapply(
    seq_len(
      nrow(
        pending_dt
      )
    ),
    function(i) {

      as.list(
        pending_dt[
          i,
          .(
            task_id,
            indicator,
            variable,
            requested_g,
            input_file
          )
        ]
      )

    }
  )

  # ----------------------------------------------------------
  # Direct V251-style parallel call.
  #
  # Crucially, output_dir and every other configuration value
  # are explicit parLapplyLB arguments. There is no anonymous
  # closure and no wrapper function.
  # ----------------------------------------------------------

  results_new <- parallel::parLapplyLB(
    cl,
    task_records,
    fun = run_one_stage1_model,
    output_dir = output_dir,
    version_tag = VERSION_TAG,
    model_degree = MODEL_DEGREE,
    screen_tol = SCREEN_TOL,
    nstart_g1 = NSTART_G1,
    nstart_g2plus = NSTART_G2PLUS,
    stagger_seconds = START_STAGGER_SECONDS,
    n_workers = N_WORKERS,
    nonbinding_maxit = NONBINDING_MAXIT,
    master_seed = MASTER_SEED,
    expected_patients = EXPECTED_PATIENTS,
    preferred_min_class_pct = PREFERRED_MIN_CLASS_PCT,
    preferred_min_appa = PREFERRED_MIN_APPA
  )

  parallel::stopCluster(
    cl
  )

  cluster_stopped <- TRUE
  cl <- NULL

} else {

  results_new <- list()

}


# ============================================================
# 15. Aggregate all task summaries
# ============================================================

summary_list <- lapply(
  seq_len(
    nrow(
      task_dt
    )
  ),
  function(i) {

    indicator <- task_dt$indicator[i]
    g <- task_dt$requested_g[i]

    f <- file.path(
      output_dir,
      paste0(
        "V267_task_summary_",
        indicator,
        "_G",
        g,
        ".csv"
      )
    )

    if (!file.exists(f)) {

      return(
        data.table(
          version = VERSION_TAG,
          stage = 1L,
          indicator = indicator,
          requested_g = as.integer(g),
          status = "NOT_FINISHED_OR_MISSING"
        )
      )

    }

    tryCatch(
      fread(
        f
      ),
      error = function(e) {

        data.table(
          version = VERSION_TAG,
          stage = 1L,
          indicator = indicator,
          requested_g = as.integer(g),
          status = "SUMMARY_READ_ERROR",
          error_message = conditionMessage(
            e
          )
        )

      }
    )

  }
)

comparison_all_dt <- rbindlist(
  summary_list,
  fill = TRUE
)

setorder(
  comparison_all_dt,
  indicator,
  requested_g
)

fwrite(
  comparison_all_dt,
  file.path(
    output_dir,
    "V267_model_comparison_all.csv"
  )
)


# ============================================================
# 16. Per-indicator comparison tables
# ============================================================

for (
  indicator_name in
    INDICATOR_MAP$indicator
) {

  indicator_comparison_dt <- comparison_all_dt[
    indicator ==
      indicator_name
  ]

  setorder(
    indicator_comparison_dt,
    requested_g
  )

  fwrite(
    indicator_comparison_dt,
    file.path(
      output_dir,
      paste0(
        "V267_model_comparison_",
        indicator_name,
        ".csv"
      )
    )
  )

}


# ============================================================
# 17. Adjacent-G IC changes within each indicator
# ============================================================

adjacent_all_list <- list()

for (
  indicator_name in
    INDICATOR_MAP$indicator
) {

  d <- comparison_all_dt[
    indicator ==
      indicator_name &
      status ==
      "CONVERGED_VALID_FULL_G"
  ][
    order(
      requested_g
    )
  ]

  needed <- c(
    "requested_g",
    "BIC",
    "AIC",
    "CAIC",
    "ssBIC",
    "HQIC"
  )

  if (
    nrow(d) > 0L &&
    all(
      needed %in%
        names(d)
    )
  ) {

    adj <- d[
      ,
      .(
        indicator,
        requested_g,
        BIC,
        delta_BIC_vs_previous_G =
          BIC -
          shift(BIC),
        AIC,
        delta_AIC_vs_previous_G =
          AIC -
          shift(AIC),
        CAIC,
        delta_CAIC_vs_previous_G =
          CAIC -
          shift(CAIC),
        ssBIC,
        delta_ssBIC_vs_previous_G =
          ssBIC -
          shift(ssBIC),
        HQIC,
        delta_HQIC_vs_previous_G =
          HQIC -
          shift(HQIC),
        min_class_pct_nonempty,
        min_APPA_nonempty,
        mean_max_posterior
      )
    ]

    adjacent_all_list[[indicator_name]] <- adj

    fwrite(
      adj,
      file.path(
        output_dir,
        paste0(
          "V267_adjacent_IC_",
          indicator_name,
          ".csv"
        )
      )
    )

  } else {

    empty_adj <- data.table()

    adjacent_all_list[[indicator_name]] <- empty_adj

    fwrite(
      empty_adj,
      file.path(
        output_dir,
        paste0(
          "V267_adjacent_IC_",
          indicator_name,
          ".csv"
        )
      )
    )

  }

}

adjacent_all_dt <- rbindlist(
  adjacent_all_list,
  fill = TRUE
)

fwrite(
  adjacent_all_dt,
  file.path(
    output_dir,
    "V267_adjacent_IC_all.csv"
  )
)


# ============================================================
# 18. Preliminary per-indicator ranking
# ============================================================

ranking_all_list <- list()

for (
  indicator_name in
    INDICATOR_MAP$indicator
) {

  d <- comparison_all_dt[
    indicator ==
      indicator_name &
      status ==
      "CONVERGED_VALID_FULL_G"
  ]

  if (
    nrow(d) > 0L &&
    "BIC" %in%
      names(d)
  ) {

    setorder(
      d,
      BIC,
      requested_g,
      na.last = TRUE
    )

  }

  ranking_all_list[[indicator_name]] <- d

  fwrite(
    d,
    file.path(
      output_dir,
      paste0(
        "V267_preliminary_ranking_",
        indicator_name,
        ".csv"
      )
    )
  )

}

ranking_all_dt <- rbindlist(
  ranking_all_list,
  fill = TRUE
)

fwrite(
  ranking_all_dt,
  file.path(
    output_dir,
    "V267_preliminary_ranking_all.csv"
  )
)


# ============================================================
# 19. Aggregate class-quality files
# ============================================================

class_quality_list <- list()

idx <- 0L

for (
  indicator_name in
    INDICATOR_MAP$indicator
) {

  for (
    g in GROUPS_TO_SCREEN
  ) {

    f <- file.path(
      output_dir,
      paste0(
        "V267_class_quality_",
        indicator_name,
        "_G",
        g,
        ".csv"
      )
    )

    if (
      file.exists(f)
    ) {

      idx <- idx + 1L

      class_quality_list[[idx]] <- tryCatch(
        fread(
          f
        ),
        error = function(e) NULL
      )

    }

  }

}

class_quality_all_dt <- rbindlist(
  class_quality_list,
  fill = TRUE
)

if (
  nrow(
    class_quality_all_dt
  ) > 0L
) {

  setorder(
    class_quality_all_dt,
    indicator,
    requested_g,
    class
  )

}

fwrite(
  class_quality_all_dt,
  file.path(
    output_dir,
    "V267_class_quality_all.csv"
  )
)


# ============================================================
# 20. Aggregate predicted trajectories
# ============================================================

trajectory_list <- list()

idx <- 0L

for (
  indicator_name in
    INDICATOR_MAP$indicator
) {

  for (
    g in GROUPS_TO_SCREEN
  ) {

    f <- file.path(
      output_dir,
      paste0(
        "V267_trajectory_",
        indicator_name,
        "_G",
        g,
        ".csv"
      )
    )

    if (
      file.exists(f)
    ) {

      idx <- idx + 1L

      trajectory_list[[idx]] <- tryCatch(
        fread(
          f
        ),
        error = function(e) NULL
      )

    }

  }

}

trajectory_all_dt <- rbindlist(
  trajectory_list,
  fill = TRUE
)

if (
  nrow(
    trajectory_all_dt
  ) > 0L
) {

  setorder(
    trajectory_all_dt,
    indicator,
    requested_g,
    class,
    time_model
  )

}

fwrite(
  trajectory_all_dt,
  file.path(
    output_dir,
    "V267_trajectory_all.csv"
  )
)


# ============================================================
# 21. Human-readable screen report
# ============================================================

master_end <- Sys.time()

elapsed_hours <- as.numeric(
  difftime(
    master_end,
    master_start,
    units = "hours"
  )
)

n_total <- nrow(
  comparison_all_dt
)

n_valid <- comparison_all_dt[
  status ==
    "CONVERGED_VALID_FULL_G",
  .N
]

n_empty <- comparison_all_dt[
  grepl(
    "EMPTY_CLASS",
    status,
    fixed = TRUE
  ),
  .N
]

n_structure_mismatch <- comparison_all_dt[
  grepl(
    "STRUCTURE_MISMATCH",
    status,
    fixed = TRUE
  ),
  .N
]

n_error <- comparison_all_dt[
  status ==
    "ERROR",
  .N
]

report_lines <- c(
  "V267 STAGE 1 SINGLE-INDICATOR GBTM SCREEN",
  "==========================================",
  "",
  paste0(
    "Start: ",
    format(
      master_start,
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste0(
    "End: ",
    format(
      master_end,
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste0(
    "Elapsed hours: ",
    round(
      elapsed_hours,
      3
    )
  ),
  "",
  "ROLE",
  "Stage 1 exploratory single-indicator GBTM.",
  "These classes are not final joint phenotypes.",
  "",
  "INDICATORS",
  "MSI",
  "Na",
  "K",
  "Cl",
  "HCO3",
  "",
  "MODEL RANGE",
  "Each indicator: G1-G6",
  "",
  "SETTINGS",
  "Polynomial degree: 2",
  "scaling: 0",
  "pruning: TRUE",
  "tol: 1",
  "G1 nstart: 1",
  "G2-G6 nstart: 3",
  "delete.empty: FALSE",
  "Automatic rerun: NONE",
  "Wall-clock timeout: NONE",
  "Practical maxit budget: NONE",
  paste0(
    "Technical maxit: ",
    NONBINDING_MAXIT
  ),
  "",
  "MISSINGNESS HANDLING",
  "For each single-indicator model, rows with NA for that indicator",
  "were removed before gbmt().",
  "No interpolation, LOCF, NOCB, or external imputation was added.",
  "",
  "RESULT COUNTS",
  paste0(
    "Total models: ",
    n_total
  ),
  paste0(
    "Converged valid full-G: ",
    n_valid
  ),
  paste0(
    "Empty-class fits: ",
    n_empty
  ),
  paste0(
    "Structure mismatch: ",
    n_structure_mismatch
  ),
  paste0(
    "Errors: ",
    n_error
  ),
  "",
  "INTERPRETATION",
  "Use BIC/other ICs, class proportions, APPA, restart behavior,",
  "and trajectory shapes together.",
  "Stage-1 results are a trajectory vocabulary to interpret Stage-2",
  "five-indicator joint models.",
  "Do not combine the five single-indicator class labels into a",
  "Cartesian-product phenotype.",
  "",
  "KEY FILES",
  "V267_model_comparison_all.csv",
  "V267_adjacent_IC_all.csv",
  "V267_preliminary_ranking_all.csv",
  "V267_class_quality_all.csv",
  "V267_trajectory_all.csv",
  "V267_indicator_input_summary.csv",
  "V267_screen_report.txt"
)

writeLines(
  report_lines,
  file.path(
    output_dir,
    "V267_screen_report.txt"
  )
)


# ============================================================
# 22. Master completion marker
# ============================================================

if (
  file.exists(
    master_running_path
  )
) {
  file.remove(
    master_running_path
  )
}

master_status <- if (
  n_valid ==
    length(
      GROUPS_TO_SCREEN
    ) *
      nrow(
        INDICATOR_MAP
      )
) {

  "ALL_STAGE1_MODELS_CONVERGED_VALID_FULL_G"

} else {

  "STAGE1_FINISHED_WITH_ITEMS_TO_REVIEW"

}

writeLines(
  c(
    paste0(
      "Version: ",
      VERSION_TAG
    ),
    paste0(
      "Status: ",
      master_status
    ),
    paste0(
      "Start: ",
      format(
        master_start,
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "End: ",
      format(
        master_end,
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "Elapsed hours: ",
      round(
        elapsed_hours,
        3
      )
    ),
    paste0(
      "Valid full-G models: ",
      n_valid,
      "/",
      n_total
    ),
    paste0(
      "Empty-class fits: ",
      n_empty
    ),
    paste0(
      "Structure mismatches: ",
      n_structure_mismatch
    ),
    paste0(
      "Errors: ",
      n_error
    )
  ),
  master_done_path
)


# ============================================================
# 23. Final console output
# ============================================================

cat("\n")
cat("============================================================\n")
cat("V267 STAGE 1 SINGLE-INDICATOR SCREEN FINISHED\n")
cat("============================================================\n\n")

cat(
  "Elapsed hours: ",
  round(
    elapsed_hours,
    3
  ),
  "\n",
  sep = ""
)

cat(
  "Valid full-G models: ",
  n_valid,
  "/",
  n_total,
  "\n\n",
  sep = ""
)

cat("Per-indicator BIC rankings among valid models:\n\n")

for (
  indicator_name in
    INDICATOR_MAP$indicator
) {

  cat(
    "---- ",
    indicator_name,
    " ----\n",
    sep = ""
  )

  d <- ranking_all_dt[
    indicator ==
      indicator_name
  ]

  cols_to_show <- intersect(
    c(
      "requested_g",
      "BIC",
      "AIC",
      "min_class_pct_nonempty",
      "min_APPA_nonempty",
      "n_converged_restarts",
      "n_restarts"
    ),
    names(d)
  )

  if (
    nrow(d) > 0L &&
    length(
      cols_to_show
    ) > 0L
  ) {

    print(
      d[
        ,
        ..cols_to_show
      ]
    )

  } else {

    cat(
      "No valid models available.\n"
    )

  }

  cat("\n")

}

cat("Output folder:\n")
cat(
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)

cat("Most useful files to send back:\n")
cat("V267_model_comparison_all.csv\n")
cat("V267_adjacent_IC_all.csv\n")
cat("V267_class_quality_all.csv\n")
cat("V267_trajectory_all.csv\n")
cat("V267_indicator_input_summary.csv\n")
cat("V267_screen_report.txt\n\n")

cat("For one specific model:\n")
cat("V267_task_summary_Indicator_GX.csv\n")
cat("V267_class_quality_Indicator_GX.csv\n")
cat("V267_EM_runs_Indicator_GX.csv\n")
cat("V267_trajectory_Indicator_GX.csv\n")
cat("V267_EM_progress_Indicator_GX.log\n")
cat("V267_screen_fit_Indicator_GX.rds\n\n")

cat("Reminder:\n")
cat("Stage-1 classes are interpretive single-indicator trajectory vocabularies.\n")
cat("They are NOT the final joint phenotypes.\n")
cat("============================================================\n")
