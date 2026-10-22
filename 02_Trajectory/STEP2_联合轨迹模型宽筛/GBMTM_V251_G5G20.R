# ============================================================
# V251 - GBMTM G5-G20 JOINT TRAJECTORY FAST SCREEN
# ============================================================
#
# Study model:
#   MSI + Na + K + Cl + HCO3 jointly enter ONE GBMTM model.
#
# Purpose:
#   Broad first-stage screening of requested class numbers G5-G20.
#
# Key settings fixed for V251:
#   - Requested G values: 5:20 (16 models)
#   - 16 PSOCK workers
#   - 10-second stagger between G values
#   - d = 2
#   - scaling = 0
#   - pruning = TRUE
#   - nstart = 1 (one randomized initialization)
#   - tol = 10 (deliberately relaxed FIRST-STAGE screening tolerance)
#   - delete.empty = FALSE
#   - NO automatic random-seed retry
#   - NO wall-clock timeout
#   - NO small/G-specific maxit budget
#   - maxit uses .Machine$integer.max only because gbmt() requires
#     a finite integer argument. It is intended to be non-binding.
#
# IMPORTANT:
#   1) A requested GX will NOT be silently converted into X-1 by this
#      script. delete.empty=FALSE prevents the gbmt empty-class deletion
#      step from intentionally removing an empty component.
#   2) If GX finishes with X stored components but one or more modal
#      groups have zero patients, the fit is SAVED and labeled
#      CONVERGED_EMPTY_CLASS. It is not accepted as a valid GX for IC
#      ranking, but it is also NOT relabeled as X-1 and NOT rerun.
#   3) If the stored model structure itself does not equal requested X
#      despite delete.empty=FALSE, the fit is saved and labeled
#      CONVERGED_STRUCTURE_MISMATCH for diagnosis.
#   4) tol=10 is for broad screening only. Final shortlisted models must
#      later be refitted with strict tolerance (e.g. 1e-4) and multiple
#      independent starts.
#
# Windows path rule:
#   Never hard-code a Windows username.
#   USERPROFILE is used as the current Windows profile root.
#   Required relative hierarchy is preserved exactly:
#   Desktop/Research/MSI/提取轨迹文件/GBMTM_stage1_output
#
# Output naming rule:
#   Version_Content_GX
#   e.g.
#     V251_task_summary_G5.csv
#     V251_class_quality_G5.csv
#     V251_EM_runs_G5.csv
#     V251_EM_progress_G5.log
#
# ============================================================


options(stringsAsFactors = FALSE)


# ============================================================
# 0. Fixed configuration
# ============================================================

VERSION_TAG <- "V251"

GROUPS_TO_SCREEN <- 5:20
N_WORKERS <- 16L

MODEL_DEGREE <- 2L

SCREEN_TOL <- 10

NSTART <- 1L

START_STAGGER_SECONDS <- 10L

# gbmt() requires finite integer maxit.
# This is a technical sentinel rather than a practical stopping budget.
NONBINDING_MAXIT <- .Machine$integer.max

MASTER_SEED <- 20260904L

EXPECTED_PATIENTS <- 9804L
EXPECTED_LONG_ROWS <- 117648L

PREFERRED_MIN_CLASS_PCT <- 5
PREFERRED_MIN_APPA <- 0.70

Z_MARKERS <- c(
  "msi_z",
  "sodium_z",
  "potassium_z",
  "chloride_z",
  "bicarbonate_z"
)

REQUIRED_COLS <- c(
  "unit_id",
  "time_model",
  Z_MARKERS
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
      "\nInstall them before running V251."
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(gbmt)
  library(parallel)
})


# ============================================================
# 2. Prevent nested threading inside each PSOCK worker
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
# 3. Verify gbmt() interface before running
# ============================================================

gbmt_formals <- names(
  formals(
    gbmt::gbmt
  )
)

required_gbmt_args <- c(
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
      as.character(packageVersion("gbmt")),
      "\nStop here. Do not run with a different interface."
    )
  )
}


# ============================================================
# 4. Resolve exact portable path using USERPROFILE
# ============================================================

if (.Platform$OS.type == "windows") {

  user_home <- Sys.getenv(
    "USERPROFILE",
    unset = ""
  )

  if (!nzchar(user_home)) {
    stop(
      paste0(
        "Windows USERPROFILE is empty.\n",
        "The script will not guess another user directory."
      )
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
  "GBMTM_V251_G5_G20_screen"
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
cat("V251 GBMTM G5-G20 JOINT TRAJECTORY FAST SCREEN\n")
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
cat("Requested G            : G5-G20\n")
cat("Workers                : 16\n")
cat("Stagger                : 10 sec/G\n")
cat("Screening tolerance    : 10\n")
cat("nstart                  : 1\n")
cat("delete.empty            : FALSE\n")
cat("Wall-clock timeout      : NONE\n")
cat("Practical maxit budget : NONE\n")
cat(
  "Technical maxit value  : ",
  NONBINDING_MAXIT,
  "\n",
  sep = ""
)
cat("Automatic retry        : NONE\n")
cat("============================================================\n\n")


# ============================================================
# 6. Validate full input
# ============================================================

cat("Reading and validating MIMIC long-format input...\n")

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

for (v in Z_MARKERS) {
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

model_input <- as.data.frame(
  full_dt[
    ,
    ..REQUIRED_COLS
  ]
)

minimal_input_rds <- file.path(
  output_dir,
  "V251_model_input.rds"
)

saveRDS(
  model_input,
  minimal_input_rds,
  compress = FALSE
)

if (!file.exists(minimal_input_rds)) {
  stop(
    "Failed to create V251_model_input.rds."
  )
}

cat(
  "Input validation passed: ",
  n_patients,
  " patients; ",
  n_rows,
  " rows.\n",
  sep = ""
)

cat(
  "Minimal worker input size: ",
  round(
    file.info(minimal_input_rds)$size / 1024^2,
    2
  ),
  " MB\n\n",
  sep = ""
)

rm(
  full_dt,
  model_input,
  duplicate_time_rows
)

gc(full = TRUE)


# ============================================================
# 7. Save configuration
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
    "N_WORKERS",
    "MODEL_DEGREE",
    "SCREEN_TOL",
    "NSTART",
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
    paste(GROUPS_TO_SCREEN, collapse = ","),
    as.character(N_WORKERS),
    as.character(MODEL_DEGREE),
    as.character(SCREEN_TOL),
    as.character(NSTART),
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
    as.character(packageVersion("gbmt")),
    as.character(packageVersion("data.table"))
  )
)

fwrite(
  config_dt,
  file.path(
    output_dir,
    "V251_config.csv"
  )
)


# ============================================================
# 8. Helper functions
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
    tolower(names(x)) ==
      tolower(target)
  )

  if (length(idx) == 0L) {
    return(NA_real_)
  }

  as.numeric(
    x[idx[1L]]
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

  as.data.table(
    as.data.frame(
      x,
      stringsAsFactors = FALSE
    )
  )

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
        converged = FALSE,
        iterations = NA_real_
      )
    )
  }

  nm <- names(em_dt)

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

  iterations <- NA_real_

  if (length(iter_col) > 0L) {

    vals <- suppressWarnings(
      as.numeric(
        em_dt[[iter_col[1L]]]
      )
    )

    vals <- vals[
      is.finite(vals)
    ]

    if (length(vals) > 0L) {
      iterations <- max(vals)
    }

  }

  converged <- FALSE

  if (length(conv_col) > 0L) {

    vals <- suppressWarnings(
      as.numeric(
        em_dt[[conv_col[1L]]]
      )
    )

    vals <- vals[
      is.finite(vals)
    ]

    if (length(vals) > 0L) {
      converged <- all(vals == 1)
    }

  }

  list(
    converged = isTRUE(converged),
    iterations = as.numeric(iterations)
  )

}


extract_class_quality <- function(
  fit,
  requested_g
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

  posterior_ncol <- ncol(post)

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
      seq_len(posterior_ncol),
      function(k) {

        idx <- which(
          assign_vec == k
        )

        n_k <- length(idx)

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
          requested_g = as.integer(requested_g),
          stored_ng = as.integer(stored_ng),
          class = as.integer(k),
          n = as.integer(n_k),
          pct = as.numeric(
            100 * n_k / length(assign_vec)
          ),
          APPA = as.numeric(appa),
          posterior_soft_mass = as.numeric(soft_mass),
          posterior_soft_mass_pct = as.numeric(
            100 * soft_mass / nrow(post)
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

  min_class_pct <- if (nrow(nonempty_dt) > 0L) {
    min(
      nonempty_dt$pct,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  min_appa <- if (
    nrow(nonempty_dt) > 0L &&
    any(is.finite(nonempty_dt$APPA))
  ) {
    min(
      nonempty_dt$APPA[
        is.finite(nonempty_dt$APPA)
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
    stored_ng == requested_g
  ) &&
    (
      posterior_ncol == requested_g
    )

  all_occupied <- (
    occupied_groups == requested_g
  )

  full_g <- structure_exact &&
    all_occupied

  summary_dt <- data.table(
    requested_g = as.integer(requested_g),
    stored_ng = as.integer(stored_ng),
    posterior_ncol = as.integer(posterior_ncol),
    occupied_groups = as.integer(occupied_groups),
    structure_exact = as.integer(structure_exact),
    all_classes_occupied = as.integer(all_occupied),
    full_g = as.integer(full_g),
    empty_classes = if (
      length(empty_classes) > 0L
    ) {
      paste(
        empty_classes,
        collapse = ","
      )
    } else {
      ""
    },
    assignment_posterior_mismatch_n = as.integer(mismatch_n),
    min_class_pct_nonempty = as.numeric(min_class_pct),
    min_APPA_nonempty = as.numeric(min_appa),
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


task_is_complete <- function(
  g,
  output_dir
) {

  summary_path <- file.path(
    output_dir,
    paste0(
      "V251_task_summary_G",
      g,
      ".csv"
    )
  )

  fit_path <- file.path(
    output_dir,
    paste0(
      "V251_screen_fit_G",
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
    fread(summary_path),
    error = function(e) NULL
  )

  if (
    is.null(d) ||
    nrow(d) < 1L ||
    !"status" %in% names(d)
  ) {
    return(FALSE)
  }

  # Any completed fit is intentionally treated as finished.
  # V251 does not silently retry empty-class solutions.
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
# 9. One-G worker
# ============================================================

run_one_g_v251 <- function(
  g,
  minimal_input_rds,
  output_dir,
  z_markers,
  model_degree,
  screen_tol,
  nstart,
  stagger_seconds_per_g,
  nonbinding_maxit,
  master_seed,
  expected_patients,
  expected_long_rows,
  preferred_min_class_pct,
  preferred_min_appa,
  minimum_g
) {

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

  data.table::setDTthreads(1L)

  fit_path <- file.path(
    output_dir,
    paste0(
      "V251_screen_fit_G",
      g,
      ".rds"
    )
  )

  partial_fit_path <- file.path(
    output_dir,
    paste0(
      "V251_screen_fit_partial_G",
      g,
      ".rds"
    )
  )

  summary_path <- file.path(
    output_dir,
    paste0(
      "V251_task_summary_G",
      g,
      ".csv"
    )
  )

  class_path <- file.path(
    output_dir,
    paste0(
      "V251_class_quality_G",
      g,
      ".csv"
    )
  )

  em_runs_path <- file.path(
    output_dir,
    paste0(
      "V251_EM_runs_G",
      g,
      ".csv"
    )
  )

  em_progress_path <- file.path(
    output_dir,
    paste0(
      "V251_EM_progress_G",
      g,
      ".log"
    )
  )

  stage_path <- file.path(
    output_dir,
    paste0(
      "V251_stage_G",
      g,
      ".log"
    )
  )

  running_path <- file.path(
    output_dir,
    paste0(
      "V251_RUNNING_G",
      g,
      ".txt"
    )
  )

  done_path <- file.path(
    output_dir,
    paste0(
      "V251_DONE_G",
      g,
      ".txt"
    )
  )

  error_path <- file.path(
    output_dir,
    paste0(
      "V251_ERROR_G",
      g,
      ".txt"
    )
  )

  # ----------------------------------------------------------
  # Resume protection
  # ----------------------------------------------------------

  if (
    file.exists(summary_path) &&
    file.exists(fit_path)
  ) {

    old <- tryCatch(
      fread(summary_path),
      error = function(e) NULL
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
      !is.null(old) &&
      nrow(old) >= 1L &&
      "status" %in% names(old) &&
      as.character(old$status[1L]) %in% completed_statuses
    ) {

      old[
        ,
        resumed_skip := 1L
      ]

      return(old)

    }

  }

  cleanup_paths <- c(
    partial_fit_path,
    running_path,
    done_path,
    error_path
  )

  for (f in cleanup_paths) {
    if (file.exists(f)) {
      file.remove(f)
    }
  }

  if (file.exists(stage_path)) {
    file.remove(stage_path)
  }

  if (file.exists(em_progress_path)) {
    file.remove(em_progress_path)
  }

  start_time <- Sys.time()

  stagger_seconds <- as.integer(
    max(
      0L,
      g - minimum_g
    ) *
      stagger_seconds_per_g
  )

  seed_used <- as.integer(
    master_seed +
      g * 100003L
  )

  writeLines(
    c(
      paste0("Version: V251"),
      paste0("Task: G", g),
      "Status: RUNNING",
      paste0(
        "Start: ",
        format(
          start_time,
          "%Y-%m-%d %H:%M:%S"
        )
      ),
      paste0(
        "PID: ",
        Sys.getpid()
      ),
      paste0(
        "Random seed: ",
        seed_used
      ),
      paste0(
        "tol: ",
        screen_tol
      ),
      paste0(
        "nstart: ",
        nstart
      ),
      "delete.empty: FALSE",
      paste0(
        "technical maxit: ",
        nonbinding_maxit
      ),
      paste0(
        "stagger: ",
        stagger_seconds,
        " sec"
      ),
      "automatic retry: NONE"
    ),
    running_path
  )

  append_stage(
    stage_path,
    "worker entered"
  )

  if (stagger_seconds > 0L) {

    append_stage(
      stage_path,
      paste0(
        "stagger sleep ",
        stagger_seconds,
        " sec"
      )
    )

    Sys.sleep(
      stagger_seconds
    )

  }

  append_stage(
    stage_path,
    "before readRDS minimal input"
  )

  model_data <- readRDS(
    minimal_input_rds
  )

  append_stage(
    stage_path,
    "after readRDS minimal input"
  )

  if (!is.data.frame(model_data)) {
    model_data <- as.data.frame(model_data)
  }

  model_data$unit_id <- as.character(
    model_data$unit_id
  )

  model_data$time_model <- as.numeric(
    model_data$time_model
  )

  for (v in z_markers) {
    model_data[[v]] <- as.numeric(
      model_data[[v]]
    )
  }

  n_units <- length(
    unique(
      model_data$unit_id
    )
  )

  n_model_rows <- nrow(
    model_data
  )

  if (n_units != expected_patients) {
    stop(
      paste0(
        "G",
        g,
        " worker patient count = ",
        n_units,
        "; expected ",
        expected_patients,
        "."
      )
    )
  }

  if (n_model_rows != expected_long_rows) {
    stop(
      paste0(
        "G",
        g,
        " worker row count = ",
        n_model_rows,
        "; expected ",
        expected_long_rows,
        "."
      )
    )
  }

  append_stage(
    stage_path,
    "worker input validation passed"
  )

  gc(full = TRUE)

  set.seed(
    seed_used
  )

  cat(
    paste0(
      "V251 G",
      g,
      " gbmt fitting started at ",
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

  fit_error <- NA_character_
  fit <- NULL

  fit_start <- Sys.time()

  sink(
    em_progress_path,
    append = TRUE,
    type = "output"
  )

  fit <- tryCatch(
    {

      gbmt::gbmt(
        x.names = z_markers,
        unit = "unit_id",
        time = "time_model",
        ng = as.integer(g),
        d = as.integer(model_degree),
        data = model_data,
        scaling = 0,
        pruning = TRUE,
        nstart = as.integer(nstart),
        tol = as.numeric(screen_tol),
        maxit = as.integer(nonbinding_maxit),
        delete.empty = FALSE,
        quiet = FALSE
      )

    },
    error = function(e) {

      fit_error <<- conditionMessage(e)

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
      requested_g = as.integer(g),
      status = "ERROR",
      error_message = as.character(fit_error),
      start_time = format(
        start_time,
        "%Y-%m-%d %H:%M:%S"
      ),
      end_time = format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      fit_seconds = as.numeric(fit_seconds),
      seed = as.integer(seed_used)
    )

    fwrite(
      error_dt,
      summary_path
    )

    if (file.exists(running_path)) {
      file.remove(running_path)
    }

    writeLines(
      c(
        paste0("Version: V251"),
        paste0("Task: G", g),
        "Status: ERROR",
        paste0(
          "Error: ",
          fit_error
        ),
        paste0(
          "End: ",
          format(
            Sys.time(),
            "%Y-%m-%d %H:%M:%S"
          )
        )
      ),
      error_path
    )

    return(error_dt)

  }

  append_stage(
    stage_path,
    "gbmt call returned"
  )

  # ----------------------------------------------------------
  # Save raw fit BEFORE classification checks.
  # This preserves the exact requested-G fit even if an empty
  # modal class is present.
  # ----------------------------------------------------------

  saveRDS(
    fit,
    partial_fit_path
  )

  if (
    !file.exists(partial_fit_path) ||
    is.na(
      file.info(partial_fit_path)$size
    ) ||
    file.info(partial_fit_path)$size <= 0
  ) {
    stop(
      paste0(
        "G",
        g,
        " model partial RDS was not written correctly."
      )
    )
  }

  if (file.exists(fit_path)) {
    file.remove(fit_path)
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
          "Could not finalize model RDS for G",
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
  # Extract EM and class diagnostics
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
    fit,
    requested_g = g
  )

  fwrite(
    class_obj$class_table,
    class_path
  )

  q <- class_obj$summary

  package_converged <- isTRUE(
    em_status$converged
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

  # ----------------------------------------------------------
  # Status logic
  # ----------------------------------------------------------

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

  # ----------------------------------------------------------
  # Information criteria
  # ----------------------------------------------------------

  final_logLik <- suppressWarnings(
    as.numeric(
      logLik(fit)
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

  if (!is.finite(AIC_value)) {
    AIC_value <- tryCatch(
      as.numeric(AIC(fit)),
      error = function(e) NA_real_
    )
  }

  if (!is.finite(BIC_value)) {
    BIC_value <- tryCatch(
      as.numeric(BIC(fit)),
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
    version = VERSION_TAG,
    requested_g = as.integer(g),
    status = as.character(status),
    package_converged = as.integer(package_converged),
    stored_ng = as.integer(q$stored_ng[1L]),
    posterior_ncol = as.integer(q$posterior_ncol[1L]),
    occupied_groups = as.integer(q$occupied_groups[1L]),
    structure_exact = as.integer(structure_exact),
    all_classes_occupied = as.integer(all_occupied),
    full_g = as.integer(full_g),
    empty_classes = as.character(q$empty_classes[1L]),
    screening_tol = as.numeric(screen_tol),
    nstart = as.integer(nstart),
    delete_empty = 0L,
    technical_maxit = as.integer(nonbinding_maxit),
    automatic_retry = 0L,
    em_iterations = as.numeric(em_status$iterations),
    final_logLik = as.numeric(final_logLik),
    AIC = as.numeric(AIC_value),
    BIC = as.numeric(BIC_value),
    CAIC = as.numeric(CAIC_value),
    ssBIC = as.numeric(ssBIC_value),
    HQIC = as.numeric(HQIC_value),
    min_class_pct_nonempty = as.numeric(min_class_pct),
    min_APPA_nonempty = as.numeric(min_appa),
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
    seed = as.integer(seed_used),
    fit_seconds = as.numeric(fit_seconds),
    total_task_seconds = as.numeric(
      difftime(
        end_time,
        start_time,
        units = "secs"
      )
    ),
    n_units = as.integer(n_units),
    n_rows = as.integer(n_model_rows),
    pid = as.integer(Sys.getpid()),
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

  if (file.exists(running_path)) {
    file.remove(running_path)
  }

  done_lines <- c(
    paste0("Version: ", VERSION_TAG),
    paste0("Task: G", g),
    paste0("Status: ", status),
    paste0(
      "Package converged: ",
      package_converged
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
      "Posterior columns: ",
      q$posterior_ncol[1L]
    ),
    paste0(
      "Occupied groups: ",
      q$occupied_groups[1L]
    ),
    paste0(
      "Empty classes: ",
      ifelse(
        nzchar(q$empty_classes[1L]),
        q$empty_classes[1L],
        "none"
      )
    ),
    paste0(
      "EM iterations: ",
      em_status$iterations
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
      "Screening tol: ",
      screen_tol
    ),
    paste0(
      "nstart: ",
      nstart
    ),
    "delete.empty: FALSE",
    "Automatic retry: NONE",
    paste0(
      "End: ",
      format(
        end_time,
        "%Y-%m-%d %H:%M:%S"
      )
    )
  )

  writeLines(
    done_lines,
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
# 10. Determine pending G values
# ============================================================

pending_groups <- GROUPS_TO_SCREEN[
  !vapply(
    GROUPS_TO_SCREEN,
    task_is_complete,
    logical(1),
    output_dir = output_dir
  )
]

cat("Pending G values: ")

if (length(pending_groups) == 0L) {

  cat("none\n\n")

} else {

  cat(
    paste0(
      "G",
      pending_groups,
      collapse = ", "
    ),
    "\n\n",
    sep = ""
  )

}


# ============================================================
# 11. Master-run markers
# ============================================================

master_start <- Sys.time()

master_running_path <- file.path(
  output_dir,
  "V251_MASTER_RUNNING.txt"
)

master_done_path <- file.path(
  output_dir,
  "V251_MASTER_DONE.txt"
)

if (file.exists(master_done_path)) {
  file.remove(master_done_path)
}

writeLines(
  c(
    "V251 GBMTM G5-G20 master run is active.",
    paste0(
      "Start: ",
      format(
        master_start,
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    "Joint model: MSI + Na + K + Cl + HCO3",
    "Requested G: 5-20",
    "Workers: 16",
    "tol: 10",
    "nstart: 1",
    "delete.empty: FALSE",
    "Automatic retry: NONE",
    "Wall-clock timeout: NONE",
    "Practical maxit budget: NONE"
  ),
  master_running_path
)


# ============================================================
# 12. Run pending models in parallel
# ============================================================

if (length(pending_groups) > 0L) {

  workers_now <- min(
    N_WORKERS,
    length(pending_groups)
  )

  cat(
    "Starting ",
    workers_now,
    " PSOCK workers...\n\n",
    sep = ""
  )

  worker_console_path <- file.path(
    output_dir,
    "V251_PS0CK_worker_console.log"
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
          parallel::stopCluster(cl),
          silent = TRUE
        )

      }

    },
    add = TRUE
  )

  parallel::clusterExport(
    cl,
    varlist = c(
      "run_one_g_v251",
      "append_stage",
      "safe_ic_value",
      "extract_em_table",
      "extract_em_status",
      "extract_class_quality",
      "VERSION_TAG"
    ),
    envir = environment()
  )

  results_new <- parallel::parLapplyLB(
    cl,
    pending_groups,
    fun = run_one_g_v251,
    minimal_input_rds = minimal_input_rds,
    output_dir = output_dir,
    z_markers = Z_MARKERS,
    model_degree = MODEL_DEGREE,
    screen_tol = SCREEN_TOL,
    nstart = NSTART,
    stagger_seconds_per_g = START_STAGGER_SECONDS,
    nonbinding_maxit = NONBINDING_MAXIT,
    master_seed = MASTER_SEED,
    expected_patients = EXPECTED_PATIENTS,
    expected_long_rows = EXPECTED_LONG_ROWS,
    preferred_min_class_pct = PREFERRED_MIN_CLASS_PCT,
    preferred_min_appa = PREFERRED_MIN_APPA,
    minimum_g = min(GROUPS_TO_SCREEN)
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
# 13. Aggregate per-G summaries
# ============================================================

summary_files <- file.path(
  output_dir,
  paste0(
    "V251_task_summary_G",
    GROUPS_TO_SCREEN,
    ".csv"
  )
)

summary_list <- lapply(
  seq_along(summary_files),
  function(i) {

    f <- summary_files[i]
    g <- GROUPS_TO_SCREEN[i]

    if (!file.exists(f)) {

      return(
        data.table(
          version = VERSION_TAG,
          requested_g = as.integer(g),
          status = "NOT_FINISHED_OR_MISSING"
        )
      )

    }

    tryCatch(
      {

        fread(f)

      },
      error = function(e) {

        data.table(
          version = VERSION_TAG,
          requested_g = as.integer(g),
          status = "SUMMARY_READ_ERROR",
          error_message = conditionMessage(e)
        )

      }
    )

  }
)

comparison_dt <- rbindlist(
  summary_list,
  fill = TRUE
)

setorder(
  comparison_dt,
  requested_g
)

fwrite(
  comparison_dt,
  file.path(
    output_dir,
    "V251_model_comparison.csv"
  )
)


# ============================================================
# 14. Compact run-status table
# ============================================================

desired_status_cols <- c(
  "version",
  "requested_g",
  "status",
  "package_converged",
  "stored_ng",
  "posterior_ncol",
  "occupied_groups",
  "empty_classes",
  "em_iterations",
  "final_logLik",
  "BIC",
  "min_class_pct_nonempty",
  "min_APPA_nonempty",
  "IC_comparison_eligible",
  "fit_seconds",
  "start_time",
  "end_time"
)

existing_status_cols <- intersect(
  desired_status_cols,
  names(comparison_dt)
)

run_status_dt <- comparison_dt[
  ,
  ..existing_status_cols
]

fwrite(
  run_status_dt,
  file.path(
    output_dir,
    "V251_run_status.csv"
  )
)


# ============================================================
# 15. Preliminary ranking
# ============================================================

if (
  all(
    c(
      "status",
      "BIC"
    ) %in%
      names(comparison_dt)
  )
) {

  preliminary_ranking_dt <- comparison_dt[
    status ==
      "CONVERGED_VALID_FULL_G"
  ]

  if (nrow(preliminary_ranking_dt) > 0L) {

    setorder(
      preliminary_ranking_dt,
      BIC,
      requested_g,
      na.last = TRUE
    )

  }

} else {

  preliminary_ranking_dt <- data.table()

}

fwrite(
  preliminary_ranking_dt,
  file.path(
    output_dir,
    "V251_preliminary_ranking.csv"
  )
)


# ============================================================
# 16. Adjacent-G IC changes
# ============================================================

adjacent_dt <- data.table()

if (
  nrow(comparison_dt) > 0L &&
  all(
    c(
      "requested_g",
      "status",
      "BIC",
      "AIC",
      "CAIC",
      "ssBIC",
      "HQIC"
    ) %in%
      names(comparison_dt)
  )
) {

  valid_for_adjacent <- comparison_dt[
    status ==
      "CONVERGED_VALID_FULL_G"
  ][
    order(requested_g)
  ]

  if (nrow(valid_for_adjacent) > 0L) {

    adjacent_dt <- valid_for_adjacent[
      ,
      .(
        requested_g,
        BIC,
        delta_BIC_vs_previous_valid_G = BIC - shift(BIC),
        AIC,
        delta_AIC_vs_previous_valid_G = AIC - shift(AIC),
        CAIC,
        delta_CAIC_vs_previous_valid_G = CAIC - shift(CAIC),
        ssBIC,
        delta_ssBIC_vs_previous_valid_G = ssBIC - shift(ssBIC),
        HQIC,
        delta_HQIC_vs_previous_valid_G = HQIC - shift(HQIC)
      )
    ]

  }

}

fwrite(
  adjacent_dt,
  file.path(
    output_dir,
    "V251_adjacent_IC_change.csv"
  )
)


# ============================================================
# 17. Human-readable report
# ============================================================

master_end <- Sys.time()

elapsed_hours <- as.numeric(
  difftime(
    master_end,
    master_start,
    units = "hours"
  )
)

n_valid <- if (
  "status" %in%
    names(comparison_dt)
) {
  comparison_dt[
    status ==
      "CONVERGED_VALID_FULL_G",
    .N
  ]
} else {
  0L
}

n_empty <- if (
  "status" %in%
    names(comparison_dt)
) {
  comparison_dt[
    grepl(
      "EMPTY_CLASS",
      status,
      fixed = TRUE
    ),
    .N
  ]
} else {
  0L
}

n_structure_mismatch <- if (
  "status" %in%
    names(comparison_dt)
) {
  comparison_dt[
    grepl(
      "STRUCTURE_MISMATCH",
      status,
      fixed = TRUE
    ),
    .N
  ]
} else {
  0L
}

n_error <- if (
  "status" %in%
    names(comparison_dt)
) {
  comparison_dt[
    status ==
      "ERROR",
    .N
  ]
} else {
  0L
}

report_lines <- c(
  "V251 GBMTM G5-G20 JOINT TRAJECTORY FAST SCREEN",
  "================================================",
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
  "MODEL",
  "MSI + Na + K + Cl + HCO3 jointly modeled in one GBMTM.",
  "",
  "SCREEN SETTINGS",
  "Requested G: G5-G20",
  "Workers: 16",
  "Stagger: 10 sec/G",
  "Polynomial degree: 2",
  "scaling: 0",
  "pruning: TRUE",
  "tol: 10",
  "nstart: 1",
  "delete.empty: FALSE",
  "Automatic seed retry: NONE",
  "Wall-clock timeout: NONE",
  "Small/G-specific maxit budget: NONE",
  paste0(
    "Technical maxit: ",
    NONBINDING_MAXIT
  ),
  "",
  "RESULT COUNTS",
  paste0(
    "Converged valid full-G: ",
    n_valid,
    "/",
    length(GROUPS_TO_SCREEN)
  ),
  paste0(
    "Empty-class fits: ",
    n_empty
  ),
  paste0(
    "Structure-mismatch fits: ",
    n_structure_mismatch
  ),
  paste0(
    "Errors: ",
    n_error
  ),
  "",
  "INTERPRETATION RULE",
  "Only CONVERGED_VALID_FULL_G fits are eligible for preliminary IC ranking.",
  "CONVERGED_EMPTY_CLASS means the requested number of components was retained",
  "but one or more modal groups had zero assigned patients.",
  "Such a fit is saved exactly as estimated and is NOT relabeled as G-1.",
  "No automatic random-seed retry is performed.",
  "",
  "IMPORTANT",
  "tol=10 is deliberately relaxed for broad first-stage screening.",
  "Do not use this run alone to lock the final class number.",
  "Shortlisted G values require strict multi-start refinement.",
  "",
  "KEY FILES",
  "V251_model_comparison.csv",
  "V251_run_status.csv",
  "V251_preliminary_ranking.csv",
  "V251_adjacent_IC_change.csv",
  "V251_screen_report.txt"
)

writeLines(
  report_lines,
  file.path(
    output_dir,
    "V251_screen_report.txt"
  )
)


# ============================================================
# 18. Master completion marker
# ============================================================

if (file.exists(master_running_path)) {
  file.remove(master_running_path)
}

master_status <- if (
  n_valid ==
    length(GROUPS_TO_SCREEN)
) {

  "ALL_G5_G20_CONVERGED_VALID_FULL_G"

} else {

  "G5_G20_FINISHED_WITH_ITEMS_TO_REVIEW"

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
      "Converged valid full-G: ",
      n_valid,
      "/",
      length(GROUPS_TO_SCREEN)
    ),
    paste0(
      "Empty-class fits: ",
      n_empty
    ),
    paste0(
      "Structure-mismatch fits: ",
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
# 19. Final console output
# ============================================================

cat("\n")
cat("============================================================\n")
cat("V251 GBMTM G5-G20 SCREEN FINISHED\n")
cat("============================================================\n\n")

print(
  run_status_dt
)

cat("\nPreliminary ranking among CONVERGED_VALID_FULL_G only:\n\n")

print(
  preliminary_ranking_dt
)

cat("\nAdjacent-G IC changes among valid full-G models:\n\n")

print(
  adjacent_dt
)

cat("\nOutput folder:\n")
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
cat("V251_model_comparison.csv\n")
cat("V251_run_status.csv\n")
cat("V251_preliminary_ranking.csv\n")
cat("V251_adjacent_IC_change.csv\n")
cat("V251_screen_report.txt\n\n")

cat("For individual G diagnostics:\n")
cat("V251_task_summary_GX.csv\n")
cat("V251_class_quality_GX.csv\n")
cat("V251_EM_runs_GX.csv\n")
cat("V251_EM_progress_GX.log\n")
cat("V251_stage_GX.log\n\n")

cat("Reminder: V251 is a broad screen with tol=10.\n")
cat("Final shortlisted models must later use strict multi-start refinement.\n")
cat("============================================================\n")
