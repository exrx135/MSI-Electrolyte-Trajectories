#!/usr/bin/env Rscript

# =============================================================================
# V292_step2_visual_G2G20.R
# FOURTH STUDY / MSI-Electrolyte GBMTM
#
# PURPOSE
#   Read the archived completed G2-G20 Stage-2 GBMTM objects from:
#
#   <current Windows profile>/Desktop/Research/MSI/提取轨迹文件/
#   GBMTM_stage2_output/STEP2 SAVE
#
#   Extract compact trajectory / assignment / posterior-quality information,
#   reproduce the native-scale and baseline-centered visualization workflow,
#   calculate level-versus-shape decomposition and adjacent-G transitions, and
#   create one uploadable ZIP.
#
# IMPORTANT
#   - NO GBMTM refitting.
#   - NO MAKE30, mortality, KRT, or other outcome data.
#   - NO source file is modified, moved, renamed, or deleted.
#   - Large RDS files are processed one G at a time.
#   - Per-G extraction supports validated checkpoint/resume.
#   - Exact archived filenames from V291 inventory are used; no fuzzy matching.
#   - The frozen MIMIC scaling values are the exact project-frozen values.
#
# REQUIRED ARCHIVE STRUCTURE
#   STEP2 SAVE/
#     01_MODEL_OBJECTS/
#       V286_STEP2_model_fit_source_V251_G2.rds
#       ...
#       V286_STEP2_model_fit_source_V251_G20.rds
#     02_PER_G_RESULTS/
#       G2/
#         V286_STEP2_class_quality_source_V251_G2.csv
#         V286_STEP2_task_summary_source_V251_G2.csv
#       ...
#       G20/
#
# OUTPUT
#   .../GBMTM_stage2_output/V292_step2_visual_G2G20/
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(zip)
})

VERSION <- "V292"
G_VALUES <- 2:20
LOW_G_FOCUS <- 2:6

INDICATORS <- c("MSI", "Na", "K", "Cl", "HCO3")
ELECTROLYTES <- c("Na", "K", "Cl", "HCO3")
REAL_ELECTROLYTE_ANCHORS <- c(3, 9, 15, 21)

EXPECTED_HOUR_GRID <- seq(1, 23, by = 2)
EXPECTED_MODEL_GRID <- (EXPECTED_HOUR_GRID - 12) / 12

Z_MARKER_TO_INDICATOR <- c(
  msi_z = "MSI",
  sodium_z = "Na",
  potassium_z = "K",
  chloride_z = "Cl",
  bicarbonate_z = "HCO3"
)

# Exact frozen project scaling values.
SCALING <- data.table(
  marker = c("msi", "sodium", "potassium", "chloride", "bicarbonate"),
  indicator = INDICATORS,
  mimic_mean = c(
    1.229422,
    138.931869,
    4.221718,
    105.125243,
    22.905997
  ),
  mimic_sd = c(
    0.3219380,
    5.3785897,
    0.6426219,
    6.5980074,
    4.7412984
  )
)

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

user_home <- normalizePath(user_home, winslash = "/", mustWork = TRUE)

stage2_base <- file.path(
  user_home,
  "Desktop", "Research", "MSI", "提取轨迹文件", "GBMTM_stage2_output"
)

archive_root <- file.path(stage2_base, "STEP2 SAVE")
model_dir <- file.path(archive_root, "01_MODEL_OBJECTS")
per_g_dir <- file.path(archive_root, "02_PER_G_RESULTS")

output_dir <- file.path(stage2_base, "V292_step2_visual_G2G20")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
resume_dir <- file.path(output_dir, "_resume")

required_dirs <- c(stage2_base, archive_root, model_dir, per_g_dir)
missing_dirs <- required_dirs[!dir.exists(required_dirs)]

if (length(missing_dirs) > 0L) {
  stop(
    paste0(
      "Required directory does not exist. Checked exact path(s):\n",
      paste(
        normalizePath(missing_dirs, winslash = "/", mustWork = FALSE),
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}

for (d in c(output_dir, table_dir, figure_dir, resume_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d)) {
    stop(
      paste0(
        "Could not create output directory:\n",
        normalizePath(d, winslash = "/", mustWork = FALSE)
      ),
      call. = FALSE
    )
  }
}

log_path <- file.path(output_dir, "V292_visual_extract_log_G2G20.txt")
zip_path <- file.path(output_dir, "V292_step2_visual_G2G20.zip")

# =============================================================================
# 2. Helpers
# =============================================================================

log_message <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_path, append = TRUE, sep = "")
  flush.console()
}

required_columns <- function(x, cols, label) {
  missing_cols <- setdiff(cols, names(x))
  if (length(missing_cols) > 0L) {
    stop(
      paste0(
        label,
        " is missing exact required field(s): ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  xx <- x[ok]
  yy <- y[ok]
  sx <- safe_sd(xx)
  sy <- safe_sd(yy)
  if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) {
    return(NA_real_)
  }
  stats::cor(xx, yy)
}

write_csv_exact <- function(x, path) {
  fwrite(x, path, na = "")
  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size <= 0
  ) {
    stop("Failed to write output file: ", path, call. = FALSE)
  }
}

md5_one <- function(path) {
  value <- tolower(unname(tools::md5sum(path)))
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop("Failed to calculate MD5: ", path, call. = FALSE)
  }
  value
}

adjusted_rand_index <- function(x, y) {
  if (length(x) != length(y)) {
    stop("ARI input vectors have different lengths.", call. = FALSE)
  }
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  if (n < 2L) return(NA_real_)

  sum_ij <- sum(choose2(tab))
  sum_i <- sum(choose2(rowSums(tab)))
  sum_j <- sum(choose2(colSums(tab)))
  total <- choose2(n)

  if (!is.finite(total) || total == 0) return(NA_real_)

  expected <- sum_i * sum_j / total
  maximum <- 0.5 * (sum_i + sum_j)
  denominator <- maximum - expected

  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  (sum_ij - expected) / denominator
}

normalised_mutual_information <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  if (n <= 0L) return(NA_real_)

  pxy <- tab / n
  px <- rowSums(pxy)
  py <- colSums(pxy)

  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(
    pxy[nz] *
      log(
        pxy[nz] /
          (px[nz[, 1L]] * py[nz[, 2L]])
      )
  )

  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  denominator <- sqrt(hx * hy)

  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  mi / denominator
}

resume_marker_is_valid <- function(marker_path, fingerprint, output_paths) {
  if (!file.exists(marker_path)) return(FALSE)
  if (any(!file.exists(output_paths))) return(FALSE)

  marker <- tryCatch(readRDS(marker_path), error = function(e) NULL)
  if (is.null(marker)) return(FALSE)
  if (!identical(marker$version, VERSION)) return(FALSE)
  if (!identical(marker$fingerprint, fingerprint)) return(FALSE)

  current_paths <- normalizePath(
    output_paths,
    winslash = "/",
    mustWork = TRUE
  )

  if (!identical(marker$output_paths, current_paths)) return(FALSE)

  current_md5 <- vapply(output_paths, md5_one, character(1))
  identical(unname(marker$output_md5), unname(current_md5))
}

# =============================================================================
# 3. Exact archived input files
# =============================================================================

model_paths <- file.path(
  model_dir,
  paste0(
    "V286_STEP2_model_fit_source_V251_G",
    G_VALUES,
    ".rds"
  )
)

quality_paths <- file.path(
  per_g_dir,
  paste0("G", G_VALUES),
  paste0(
    "V286_STEP2_class_quality_source_V251_G",
    G_VALUES,
    ".csv"
  )
)

task_paths <- file.path(
  per_g_dir,
  paste0("G", G_VALUES),
  paste0(
    "V286_STEP2_task_summary_source_V251_G",
    G_VALUES,
    ".csv"
  )
)

all_required_files <- c(model_paths, quality_paths, task_paths)
missing_files <- all_required_files[!file.exists(all_required_files)]

if (length(missing_files) > 0L) {
  stop(
    paste0(
      "Required archived Stage-2 file(s) are missing.\n",
      "No alternate filename was guessed.\n",
      "Checked exact path(s):\n",
      paste(
        normalizePath(missing_files, winslash = "/", mustWork = FALSE),
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}

source_inventory <- data.table(
  requested_g = G_VALUES,
  model_file = normalizePath(model_paths, winslash = "/", mustWork = TRUE),
  model_size_bytes = as.numeric(file.info(model_paths)$size),
  model_size_mb = as.numeric(file.info(model_paths)$size) / 1024^2,
  model_modified = format(
    file.info(model_paths)$mtime,
    "%Y-%m-%d %H:%M:%S"
  ),
  class_quality_file = normalizePath(
    quality_paths,
    winslash = "/",
    mustWork = TRUE
  ),
  task_summary_file = normalizePath(
    task_paths,
    winslash = "/",
    mustWork = TRUE
  )
)

write_csv_exact(
  source_inventory,
  file.path(table_dir, "V292_source_inventory_G2G20.csv")
)

scaling_bridge <- copy(SCALING)
scaling_bridge[
  ,
  conversion_formula :=
    "fitted_native = mimic_mean + mimic_sd * fitted_z"
]
scaling_bridge[
  ,
  scaling_source :=
    "FOURTH STUDY frozen MIMIC development-cohort parameters"
]

write_csv_exact(
  scaling_bridge,
  file.path(table_dir, "V292_scaling_bridge_G2G20.csv")
)

# =============================================================================
# 4. Sequential extraction with per-G resume
# =============================================================================

trajectory_list <- list()
quality_list <- list()
summary_list <- list()
assignment_list <- list()

log_message(VERSION, " started")
log_message(
  "Archive root: ",
  normalizePath(archive_root, winslash = "/", mustWork = TRUE)
)
log_message(
  "Output directory: ",
  normalizePath(output_dir, winslash = "/", mustWork = TRUE)
)

for (g in G_VALUES) {
  key <- as.character(g)

  model_path <- file.path(
    model_dir,
    paste0(
      "V286_STEP2_model_fit_source_V251_G",
      g,
      ".rds"
    )
  )

  official_quality_path <- file.path(
    per_g_dir,
    paste0("G", g),
    paste0(
      "V286_STEP2_class_quality_source_V251_G",
      g,
      ".csv"
    )
  )

  task_summary_path <- file.path(
    per_g_dir,
    paste0("G", g),
    paste0(
      "V286_STEP2_task_summary_source_V251_G",
      g,
      ".csv"
    )
  )

  trajectory_path <- file.path(
    table_dir,
    paste0("V292_trajectory_G", g, ".csv")
  )

  quality_path <- file.path(
    table_dir,
    paste0("V292_class_quality_G", g, ".csv")
  )

  summary_path <- file.path(
    table_dir,
    paste0("V292_model_summary_G", g, ".csv")
  )

  assignment_path <- file.path(
    table_dir,
    paste0("V292_assignment_G", g, ".csv")
  )

  marker_path <- file.path(
    resume_dir,
    paste0("V292_extract_done_G", g, ".rds")
  )

  output_paths <- c(
    trajectory_path,
    quality_path,
    summary_path,
    assignment_path
  )

  fingerprint <- list(
    requested_g = as.integer(g),
    model_path = normalizePath(
      model_path,
      winslash = "/",
      mustWork = TRUE
    ),
    model_size_bytes = as.numeric(file.info(model_path)$size),
    model_mtime_num = as.numeric(file.info(model_path)$mtime),
    quality_path = normalizePath(
      official_quality_path,
      winslash = "/",
      mustWork = TRUE
    ),
    quality_size_bytes = as.numeric(file.info(official_quality_path)$size),
    quality_mtime_num = as.numeric(file.info(official_quality_path)$mtime),
    task_path = normalizePath(
      task_summary_path,
      winslash = "/",
      mustWork = TRUE
    ),
    task_size_bytes = as.numeric(file.info(task_summary_path)$size),
    task_mtime_num = as.numeric(file.info(task_summary_path)$mtime),
    frozen_scaling = paste(
      paste0(
        scaling_bridge$indicator,
        "=",
        scaling_bridge$mimic_mean,
        "/",
        scaling_bridge$mimic_sd
      ),
      collapse = "|"
    )
  )

  if (
    resume_marker_is_valid(
      marker_path,
      fingerprint,
      output_paths
    )
  ) {
    log_message(
      "G",
      g,
      " valid compact extraction found; large RDS read skipped"
    )

    trajectory_list[[key]] <- fread(trajectory_path)
    quality_list[[key]] <- fread(quality_path)
    summary_list[[key]] <- fread(summary_path)
    assignment_list[[key]] <- fread(assignment_path)
    next
  }

  log_message(
    "Reading G",
    g,
    " RDS (",
    sprintf("%.1f", file.info(model_path)$size / 1024^2),
    " MB)"
  )

  extraction_start <- Sys.time()
  fit <- readRDS(model_path)

  if (!inherits(fit, "gbmt")) {
    stop(
      "G",
      g,
      " archived object does not inherit from class 'gbmt'.",
      call. = FALSE
    )
  }

  required_fit_fields <- c(
    "call",
    "posterior",
    "fitted",
    "assign",
    "logLik",
    "npar",
    "em",
    "data.orig"
  )

  missing_fit_fields <- setdiff(required_fit_fields, names(fit))
  if (length(missing_fit_fields) > 0L) {
    stop(
      paste0(
        "G",
        g,
        " GBMT object is missing exact required field(s): ",
        paste(missing_fit_fields, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  unit_name <- as.character(fit$call$unit)[1L]
  time_name <- as.character(fit$call$time)[1L]
  x_names <- as.character(fit$call$x.names)
  stored_ng <- as.integer(fit$call$ng)[1L]

  if (!identical(x_names, names(Z_MARKER_TO_INDICATOR))) {
    stop(
      paste0(
        "G",
        g,
        " call$x.names differs from the frozen marker order.\n",
        "Expected: ",
        paste(names(Z_MARKER_TO_INDICATOR), collapse = ", "),
        "\nObserved: ",
        paste(x_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (
    !identical(unit_name, "unit_id") ||
    !identical(time_name, "time_model")
  ) {
    stop(
      paste0(
        "G",
        g,
        " call metadata differs from frozen model metadata.\n",
        "unit = ",
        unit_name,
        "; time = ",
        time_name
      ),
      call. = FALSE
    )
  }

  data_orig <- as.data.table(fit$data.orig)
  required_columns(
    data_orig,
    c(unit_name, time_name, x_names),
    paste0("G", g, " fit$data.orig")
  )

  posterior <- as.matrix(fit$posterior)
  storage.mode(posterior) <- "double"
  assigned_class <- as.integer(fit$assign)

  if (
    stored_ng != g ||
    ncol(posterior) != g ||
    nrow(posterior) != length(assigned_class) ||
    length(fit$fitted) != g
  ) {
    stop(
      paste0(
        "G",
        g,
        " model structure validation failed: stored_ng=",
        stored_ng,
        ", posterior=",
        nrow(posterior),
        "x",
        ncol(posterior),
        ", assignment_length=",
        length(assigned_class),
        ", fitted_classes=",
        length(fit$fitted)
      ),
      call. = FALSE
    )
  }

  unit_factor <- factor(data_orig[[unit_name]])
  unit_ids <- levels(unit_factor)

  if (length(unit_ids) != nrow(posterior)) {
    stop(
      paste0(
        "G",
        g,
        " unit ordering validation failed: factor levels=",
        length(unit_ids),
        ", posterior rows=",
        nrow(posterior)
      ),
      call. = FALSE
    )
  }

  modal_class <- max.col(posterior, ties.method = "first")
  max_posterior <- apply(
    posterior,
    1L,
    max,
    na.rm = TRUE
  )

  posterior_for_entropy <- pmax(
    posterior,
    .Machine$double.xmin
  )

  entropy_raw <- -rowSums(
    posterior_for_entropy *
      log(posterior_for_entropy)
  )

  entropy_normalised <- entropy_raw / log(g)

  assignment_dt <- data.table(
    requested_g = g,
    unit_id = as.character(unit_ids),
    assigned_class = assigned_class,
    modal_class = as.integer(modal_class),
    assignment_modal_mismatch = as.integer(
      assigned_class != modal_class
    ),
    max_posterior = as.numeric(max_posterior),
    posterior_entropy_normalised = as.numeric(entropy_normalised)
  )

  quality_dt <- rbindlist(
    lapply(
      seq_len(g),
      function(k) {
        idx <- which(assigned_class == k)
        n_k <- length(idx)

        data.table(
          requested_g = g,
          class = k,
          n = n_k,
          class_pct = 100 * n_k / length(assigned_class),
          APPA_model = if (n_k > 0L) {
            mean(posterior[idx, k], na.rm = TRUE)
          } else {
            NA_real_
          },
          posterior_soft_mass = sum(
            posterior[, k],
            na.rm = TRUE
          ),
          posterior_soft_mass_pct =
            100 *
              sum(posterior[, k], na.rm = TRUE) /
              nrow(posterior),
          mean_max_posterior_in_assigned_class =
            if (n_k > 0L) {
              mean(max_posterior[idx])
            } else {
              NA_real_
            },
          median_max_posterior_in_assigned_class =
            if (n_k > 0L) {
              stats::median(max_posterior[idx])
            } else {
              NA_real_
            }
        )
      }
    )
  )

  quality_dt[
    ,
    class_label := paste0(
      "C",
      class,
      " (",
      sprintf("%.2f", class_pct),
      "%)"
    )
  ]

  # Exact archived V251-compatible class-quality structure.
  official_quality <- fread(official_quality_path)

  required_columns(
    official_quality,
    c(
      "requested_g",
      "class",
      "n",
      "pct",
      "APPA"
    ),
    paste0(
      "Archived class quality G",
      g
    )
  )

  if (
    nrow(official_quality) != g ||
    any(as.integer(official_quality$requested_g) != g)
  ) {
    stop(
      "G",
      g,
      " archived class-quality row structure is invalid.",
      call. = FALSE
    )
  }

  official_quality[
    ,
    `:=`(
      class = as.integer(class),
      n = as.integer(n),
      pct = as.numeric(pct),
      APPA = as.numeric(APPA)
    )
  ]

  quality_check <- merge(
    quality_dt[
      ,
      .(
        class,
        n_extracted = n,
        pct_extracted = class_pct,
        APPA_extracted = APPA_model
      )
    ],
    official_quality[
      ,
      .(
        class,
        n_official = n,
        pct_official = pct,
        APPA_official = APPA
      )
    ],
    by = "class",
    all = TRUE,
    sort = TRUE
  )

  if (
    nrow(quality_check) != g ||
    any(quality_check$n_extracted != quality_check$n_official) ||
    max(
      abs(
        quality_check$pct_extracted -
          quality_check$pct_official
      ),
      na.rm = TRUE
    ) > 1e-8 ||
    max(
      abs(
        quality_check$APPA_extracted -
          quality_check$APPA_official
      ),
      na.rm = TRUE
    ) > 1e-8
  ) {
    stop(
      paste0(
        "G",
        g,
        " extracted posterior/class quality does not reproduce ",
        "the archived official class-quality output."
      ),
      call. = FALSE
    )
  }

  fitted_parts <- vector(
    "list",
    g * length(x_names)
  )

  part_index <- 0L

  for (k in seq_len(g)) {
    fitted_k <- as.data.frame(
      fit$fitted[[k]],
      stringsAsFactors = FALSE
    )

    required_columns(
      fitted_k,
      x_names,
      paste0("G", g, " fitted class ", k)
    )

    model_time <- safe_numeric(rownames(fitted_k))

    if (
      length(model_time) != nrow(fitted_k) ||
      any(!is.finite(model_time)) ||
      anyDuplicated(model_time) > 0L ||
      !isTRUE(
        all.equal(
          model_time,
          EXPECTED_MODEL_GRID,
          tolerance = 1e-12,
          check.attributes = FALSE
        )
      )
    ) {
      stop(
        paste0(
          "G",
          g,
          " fitted class ",
          k,
          " does not use the exact frozen time grid."
        ),
        call. = FALSE
      )
    }

    for (marker in x_names) {
      part_index <- part_index + 1L

      fitted_parts[[part_index]] <- data.table(
        requested_g = g,
        class = k,
        marker = marker,
        indicator = unname(
          Z_MARKER_TO_INDICATOR[marker]
        ),
        model_time = model_time,
        time_hour_center = EXPECTED_HOUR_GRID,
        fitted_z = safe_numeric(fitted_k[[marker]])
      )
    }
  }

  trajectory_dt <- rbindlist(
    fitted_parts,
    use.names = TRUE
  )

  if (
    any(!is.finite(trajectory_dt$fitted_z)) ||
    !setequal(
      unique(trajectory_dt$time_hour_center),
      EXPECTED_HOUR_GRID
    )
  ) {
    stop(
      "G",
      g,
      " extracted fitted trajectories failed numeric/time validation.",
      call. = FALSE
    )
  }

  trajectory_dt <- merge(
    trajectory_dt,
    scaling_bridge[
      ,
      .(
        indicator,
        mimic_mean,
        mimic_sd
      )
    ],
    by = "indicator",
    all.x = TRUE,
    sort = FALSE
  )

  if (
    any(!is.finite(trajectory_dt$mimic_mean)) ||
    any(!is.finite(trajectory_dt$mimic_sd))
  ) {
    stop(
      "G",
      g,
      " frozen scaling join failed.",
      call. = FALSE
    )
  }

  trajectory_dt[
    ,
    fitted_native :=
      mimic_mean +
      mimic_sd * fitted_z
  ]

  trajectory_dt[
    ,
    c("mimic_mean", "mimic_sd") := NULL
  ]

  trajectory_dt <- merge(
    trajectory_dt,
    quality_dt[
      ,
      .(
        class,
        class_label,
        n,
        class_pct,
        APPA_model
      )
    ],
    by = "class",
    all.x = TRUE,
    sort = FALSE
  )

  trajectory_dt[
    ,
    indicator_order := match(
      indicator,
      INDICATORS
    )
  ]

  setorder(
    trajectory_dt,
    class,
    indicator_order,
    time_hour_center
  )

  trajectory_dt[
    ,
    indicator_order := NULL
  ]

  setcolorder(
    trajectory_dt,
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
      "fitted_native"
    )
  )

  task_summary <- fread(task_summary_path)

  if (nrow(task_summary) != 1L) {
    stop(
      "G",
      g,
      " task summary must contain exactly one data row.",
      call. = FALSE
    )
  }

  required_columns(
    task_summary,
    c(
      "requested_g",
      "status",
      "n_units",
      "BIC",
      "min_class_pct_nonempty",
      "min_APPA_nonempty"
    ),
    paste0("Archived task summary G", g)
  )

  if (
    as.integer(task_summary$requested_g[1L]) != g
  ) {
    stop(
      "G",
      g,
      " task-summary requested_g does not match its directory.",
      call. = FALSE
    )
  }

  summary_dt <- data.table(
    requested_g = g,
    n_patients = as.integer(
      task_summary$n_units[1L]
    ),
    status = as.character(
      task_summary$status[1L]
    ),
    final_logLik = as.numeric(fit$logLik),
    npar = as.numeric(fit$npar),
    BIC = as.numeric(task_summary$BIC[1L]),
    min_class_pct = as.numeric(
      task_summary$min_class_pct_nonempty[1L]
    ),
    min_APPA = as.numeric(
      task_summary$min_APPA_nonempty[1L]
    ),
    posterior_rows = nrow(posterior),
    posterior_columns = ncol(posterior),
    assignment_modal_mismatch_n = sum(
      assigned_class != modal_class
    ),
    mean_max_posterior = mean(max_posterior),
    median_max_posterior =
      stats::median(max_posterior),
    mean_posterior_entropy_normalised =
      mean(entropy_normalised),
    em_iterations = if (
      !is.null(fit$em[["n.iter"]])
    ) {
      as.integer(fit$em[["n.iter"]])
    } else {
      NA_integer_
    },
    package_converged = if (
      !is.null(fit$em[["converged"]])
    ) {
      as.integer(
        as.logical(
          fit$em[["converged"]]
        )[1L]
      )
    } else {
      NA_integer_
    },
    source_model_file = normalizePath(
      model_path,
      winslash = "/",
      mustWork = TRUE
    ),
    extraction_seconds = as.numeric(
      difftime(
        Sys.time(),
        extraction_start,
        units = "secs"
      )
    )
  )

  if (
    summary_dt$n_patients != nrow(posterior)
  ) {
    stop(
      "G",
      g,
      " task-summary n_units differs from posterior row count.",
      call. = FALSE
    )
  }

  write_csv_exact(
    trajectory_dt,
    trajectory_path
  )
  write_csv_exact(
    quality_dt,
    quality_path
  )
  write_csv_exact(
    summary_dt,
    summary_path
  )
  write_csv_exact(
    assignment_dt,
    assignment_path
  )

  marker <- list(
    version = VERSION,
    completed_at = format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    fingerprint = fingerprint,
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

  saveRDS(marker, marker_path)

  trajectory_list[[key]] <- trajectory_dt
  quality_list[[key]] <- quality_dt
  summary_list[[key]] <- summary_dt
  assignment_list[[key]] <- assignment_dt

  rm(
    fit,
    data_orig,
    posterior,
    posterior_for_entropy,
    assigned_class,
    modal_class,
    max_posterior,
    entropy_raw,
    entropy_normalised,
    assignment_dt,
    quality_dt,
    official_quality,
    quality_check,
    fitted_parts,
    fitted_k,
    trajectory_dt,
    task_summary,
    summary_dt,
    unit_factor,
    unit_ids
  )

  invisible(gc(full = TRUE))

  log_message(
    "G",
    g,
    " extraction completed and resume marker saved"
  )
}

# =============================================================================
# 5. Combined compact tables
# =============================================================================

all_trajectories <- rbindlist(
  trajectory_list,
  use.names = TRUE,
  fill = TRUE
)

all_quality <- rbindlist(
  quality_list,
  use.names = TRUE,
  fill = TRUE
)

all_summaries <- rbindlist(
  summary_list,
  use.names = TRUE,
  fill = TRUE
)

all_assignments <- rbindlist(
  assignment_list,
  use.names = TRUE,
  fill = TRUE
)

all_trajectories[
  ,
  indicator_order := match(
    indicator,
    INDICATORS
  )
]

setorder(
  all_trajectories,
  requested_g,
  class,
  indicator_order,
  time_hour_center
)

all_trajectories[
  ,
  indicator_order := NULL
]

all_trajectories[
  ,
  centered_z :=
    fitted_z - first(fitted_z),
  by = .(
    requested_g,
    class,
    indicator
  )
]

setorder(
  all_quality,
  requested_g,
  class
)

setorder(
  all_summaries,
  requested_g
)

setorder(
  all_assignments,
  requested_g,
  unit_id
)

write_csv_exact(
  all_trajectories,
  file.path(
    table_dir,
    "V292_trajectory_all_G2G20.csv"
  )
)

write_csv_exact(
  all_quality,
  file.path(
    table_dir,
    "V292_class_quality_all_G2G20.csv"
  )
)

write_csv_exact(
  all_summaries,
  file.path(
    table_dir,
    "V292_model_summary_all_G2G20.csv"
  )
)

write_csv_exact(
  all_assignments,
  file.path(
    table_dir,
    "V292_assignment_all_G2G20.csv"
  )
)

trajectory_summary <- all_trajectories[
  ,
  .(
    n = first(n),
    class_pct = first(class_pct),
    APPA_model = first(APPA_model),
    start_native = first(fitted_native),
    end_native = last(fitted_native),
    delta_native =
      last(fitted_native) -
      first(fitted_native),
    mean_native = mean(fitted_native),
    min_native = min(fitted_native),
    max_native = max(fitted_native),
    start_z = first(fitted_z),
    end_z = last(fitted_z),
    delta_z =
      last(fitted_z) -
      first(fitted_z),
    centered_shape_range_z =
      max(centered_z) -
      min(centered_z)
  ),
  by = .(
    requested_g,
    class,
    indicator
  )
]

write_csv_exact(
  trajectory_summary,
  file.path(
    table_dir,
    "V292_trajectory_summary_G2G20.csv"
  )
)

# =============================================================================
# 6. Pairwise level-versus-shape decomposition
# =============================================================================

pairwise_indicator_parts <- list()
pair_idx <- 0L

for (g in G_VALUES) {
  classes <- sort(
    unique(
      all_trajectories[
        requested_g == g,
        class
      ]
    )
  )

  if (length(classes) != g) {
    stop(
      "Trajectory class count mismatch for G",
      g,
      ".",
      call. = FALSE
    )
  }

  class_pairs <- utils::combn(
    classes,
    2L
  )

  for (j in seq_len(ncol(class_pairs))) {
    class_a <- class_pairs[1L, j]
    class_b <- class_pairs[2L, j]

    for (indicator_now in INDICATORS) {
      a <- all_trajectories[
        requested_g == g &
          class == class_a &
          indicator == indicator_now
      ][order(time_hour_center)]

      b <- all_trajectories[
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
            "Trajectory alignment failed for G",
            g,
            " C",
            class_a,
            " versus C",
            class_b,
            " indicator ",
            indicator_now,
            "."
          ),
          call. = FALSE
        )
      }

      diff_z <- a$fitted_z - b$fitted_z
      diff_native <-
        a$fitted_native -
        b$fitted_native

      centered_a <-
        a$fitted_z -
        a$fitted_z[1L]

      centered_b <-
        b$fitted_z -
        b$fitted_z[1L]

      centered_diff <-
        centered_a -
        centered_b

      pair_idx <- pair_idx + 1L

      pairwise_indicator_parts[[pair_idx]] <- data.table(
        requested_g = g,
        class_a = class_a,
        class_b = class_b,
        indicator = indicator_now,
        mean_level_offset_z =
          mean(diff_z),
        absolute_level_offset_z =
          abs(mean(diff_z)),
        mean_level_offset_native =
          mean(diff_native),
        absolute_level_offset_native =
          abs(mean(diff_native)),
        total_RMS_difference_z =
          sqrt(mean(diff_z^2)),
        baseline_centered_shape_RMS_z =
          sqrt(mean(centered_diff^2)),
        baseline_centered_shape_max_abs_z =
          max(abs(centered_diff)),
        baseline_centered_shape_correlation =
          safe_cor(
            centered_a,
            centered_b
          )
      )
    }
  }
}

level_shape_indicator <- rbindlist(
  pairwise_indicator_parts,
  use.names = TRUE,
  fill = TRUE
)

level_shape_joint <- level_shape_indicator[
  ,
  .(
    indicator_count = .N,
    joint_level_RMS_z =
      sqrt(
        mean(
          mean_level_offset_z^2
        )
      ),
    joint_total_RMS_z =
      sqrt(
        mean(
          total_RMS_difference_z^2
        )
      ),
    joint_centered_shape_RMS_z =
      sqrt(
        mean(
          baseline_centered_shape_RMS_z^2
        )
      ),
    mean_centered_shape_correlation =
      safe_mean(
        baseline_centered_shape_correlation
      )
  ),
  by = .(
    requested_g,
    class_a,
    class_b
  )
]

level_shape_joint[
  ,
  shape_share_of_level_plus_shape :=
    fifelse(
      (
        joint_level_RMS_z +
          joint_centered_shape_RMS_z
      ) > 0,
      joint_centered_shape_RMS_z /
        (
          joint_level_RMS_z +
            joint_centered_shape_RMS_z
        ),
      NA_real_
    )
]

level_shape_indicator[
  ,
  indicator_order := match(
    indicator,
    INDICATORS
  )
]

setorder(
  level_shape_indicator,
  requested_g,
  class_a,
  class_b,
  indicator_order
)

level_shape_indicator[
  ,
  indicator_order := NULL
]

setorder(
  level_shape_joint,
  requested_g,
  class_a,
  class_b
)

nearest_shape_pair <- level_shape_joint[
  order(
    requested_g,
    joint_centered_shape_RMS_z,
    -mean_centered_shape_correlation,
    class_a,
    class_b
  ),
  .SD[1L],
  by = requested_g
]

write_csv_exact(
  level_shape_indicator,
  file.path(
    table_dir,
    "V292_level_shape_by_indicator_G2G20.csv"
  )
)

write_csv_exact(
  level_shape_joint,
  file.path(
    table_dir,
    "V292_level_shape_joint_G2G20.csv"
  )
)

write_csv_exact(
  nearest_shape_pair,
  file.path(
    table_dir,
    "V292_nearest_shape_pair_G2G20.csv"
  )
)

# =============================================================================
# 7. Adjacent-G transitions + ARI/NMI
# =============================================================================

transition_parts <- list()
metric_parts <- list()
transition_idx <- 0L
metric_idx <- 0L

for (g in G_VALUES[-length(G_VALUES)]) {
  g2 <- g + 1L

  a <- all_assignments[
    requested_g == g,
    .(
      unit_id,
      class_from = assigned_class
    )
  ]

  b <- all_assignments[
    requested_g == g2,
    .(
      unit_id,
      class_to = assigned_class
    )
  ]

  if (
    nrow(a) == 0L ||
    nrow(b) == 0L ||
    anyDuplicated(a$unit_id) > 0L ||
    anyDuplicated(b$unit_id) > 0L
  ) {
    stop(
      "Invalid assignment table for G",
      g,
      " -> G",
      g2,
      ".",
      call. = FALSE
    )
  }

  m <- merge(
    a,
    b,
    by = "unit_id",
    all = FALSE,
    sort = FALSE
  )

  if (
    nrow(m) != nrow(a) ||
    nrow(m) != nrow(b)
  ) {
    stop(
      "Patient set mismatch for G",
      g,
      " -> G",
      g2,
      ".",
      call. = FALSE
    )
  }

  trans <- m[
    ,
    .(N = .N),
    by = .(
      class_from,
      class_to
    )
  ]

  trans[
    ,
    row_total := sum(N),
    by = class_from
  ]

  trans[
    ,
    column_total := sum(N),
    by = class_to
  ]

  trans[
    ,
    `:=`(
      from_g = g,
      to_g = g2,
      row_pct =
        100 * N / row_total,
      column_pct =
        100 * N / column_total
    )
  ]

  setcolorder(
    trans,
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
    )
  )

  setorder(
    trans,
    class_from,
    -N,
    class_to
  )

  transition_idx <- transition_idx + 1L
  transition_parts[[transition_idx]] <- trans

  metric_idx <- metric_idx + 1L
  metric_parts[[metric_idx]] <- data.table(
    from_g = g,
    to_g = g2,
    n_units = nrow(m),
    adjusted_rand_index =
      adjusted_rand_index(
        m$class_from,
        m$class_to
      ),
    normalised_mutual_information =
      normalised_mutual_information(
        m$class_from,
        m$class_to
      )
  )
}

transition_all <- rbindlist(
  transition_parts,
  use.names = TRUE,
  fill = TRUE
)

transition_metrics <- rbindlist(
  metric_parts,
  use.names = TRUE,
  fill = TRUE
)

write_csv_exact(
  transition_all,
  file.path(
    table_dir,
    "V292_transition_G2G20.csv"
  )
)

write_csv_exact(
  transition_metrics,
  file.path(
    table_dir,
    "V292_transition_metrics_G2G20.csv"
  )
)

# =============================================================================
# 8. Visualization
# =============================================================================

plot_dt <- copy(all_trajectories)

plot_dt[
  ,
  indicator := factor(
    indicator,
    levels = INDICATORS
  )
]

plot_dt[
  ,
  g_label := factor(
    paste0("G", requested_g),
    levels = paste0("G", G_VALUES)
  )
]

# Global per-indicator bounds across G2-G20, preserving the V274 visual logic.
native_bounds <- plot_dt[
  ,
  .(
    lo = min(fitted_native),
    hi = max(fitted_native)
  ),
  by = indicator
]

native_bounds <- melt(
  native_bounds,
  id.vars = "indicator",
  measure.vars = c("lo", "hi"),
  value.name = "fitted_native"
)

native_bounds[
  ,
  time_hour_center := 1
]

centered_bounds <- plot_dt[
  ,
  .(
    lo = min(centered_z),
    hi = max(centered_z)
  ),
  by = indicator
]

centered_bounds <- melt(
  centered_bounds,
  id.vars = "indicator",
  measure.vars = c("lo", "hi"),
  value.name = "centered_z"
)

centered_bounds[
  ,
  time_hour_center := 1
]

p_native <- function(g) {
  d <- plot_dt[requested_g == g]
  q <- all_quality[requested_g == g]
  sm <- all_summaries[requested_g == g]

  class_levels <- q[
    order(class),
    class_label
  ]

  d[
    ,
    class_label := factor(
      class_label,
      levels = class_levels
    )
  ]

  pts <- d[
    indicator == "MSI" |
      (
        as.character(indicator) %in% ELECTROLYTES &
          time_hour_center %in% REAL_ELECTROLYTE_ANCHORS
      )
  ]

  ggplot(
    d,
    aes(
      x = time_hour_center,
      y = fitted_native,
      group = class_label,
      colour = class_label
    )
  ) +
    geom_blank(
      data = native_bounds,
      aes(
        x = time_hour_center,
        y = fitted_native
      ),
      inherit.aes = FALSE
    ) +
    geom_line(
      linewidth = 0.78,
      alpha = 0.90
    ) +
    geom_point(
      data = pts,
      aes(
        x = time_hour_center,
        y = fitted_native,
        group = class_label,
        colour = class_label
      ),
      inherit.aes = FALSE,
      size = 1.6
    ) +
    facet_wrap(
      ~indicator,
      scales = "free_y",
      ncol = 2
    ) +
    scale_x_continuous(
      breaks = EXPECTED_HOUR_GRID
    ) +
    labs(
      title = paste0(
        "Stage 2 Joint Trajectories: G",
        g,
        " — Native Scale"
      ),
      subtitle = paste0(
        "n=",
        format(sm$n_patients[1L], big.mark = ","),
        " | BIC=",
        sprintf("%.1f", sm$BIC[1L]),
        " | min class=",
        sprintf("%.2f", sm$min_class_pct[1L]),
        "% | min APPA=",
        sprintf("%.3f", sm$min_APPA[1L])
      ),
      x = "Time window center (hours)",
      y = "Fitted value (native scale)",
      colour = "Class",
      caption = paste0(
        "Electrolyte larger points mark real 3/9/15/21 h anchors; ",
        "lines show fitted trajectories."
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 7.8),
      strip.background = element_rect(fill = "grey93"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(
        hjust = 0,
        size = 8.2
      )
    )
}

p_centered <- function(g) {
  d <- plot_dt[requested_g == g]
  q <- all_quality[requested_g == g]

  class_levels <- q[
    order(class),
    class_label
  ]

  d[
    ,
    class_label := factor(
      class_label,
      levels = class_levels
    )
  ]

  pts <- d[
    indicator == "MSI" |
      (
        as.character(indicator) %in% ELECTROLYTES &
          time_hour_center %in% REAL_ELECTROLYTE_ANCHORS
      )
  ]

  ggplot(
    d,
    aes(
      x = time_hour_center,
      y = centered_z,
      group = class_label,
      colour = class_label
    )
  ) +
    geom_blank(
      data = centered_bounds,
      aes(
        x = time_hour_center,
        y = centered_z
      ),
      inherit.aes = FALSE
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_line(
      linewidth = 0.78,
      alpha = 0.90
    ) +
    geom_point(
      data = pts,
      aes(
        x = time_hour_center,
        y = centered_z,
        group = class_label,
        colour = class_label
      ),
      inherit.aes = FALSE,
      size = 1.6
    ) +
    facet_wrap(
      ~indicator,
      scales = "free_y",
      ncol = 2
    ) +
    scale_x_continuous(
      breaks = EXPECTED_HOUR_GRID
    ) +
    labs(
      title = paste0(
        "Stage 2 Joint Trajectories: G",
        g,
        " — Baseline-Centered"
      ),
      subtitle = paste0(
        "Each class/indicator subtracts its first fitted z-value. ",
        "Persistent vertical level offsets are removed to emphasize temporal shape."
      ),
      x = "Time window center (hours)",
      y = "Change from first fitted point (SD units)",
      colour = "Class",
      caption = paste0(
        "Diagnostic view only; centering does not refit or alter the GBMTM. ",
        "Electrolyte clinical interpretation remains anchored to 3/9/15/21 h."
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 7.8),
      strip.background = element_rect(fill = "grey93"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(
        hjust = 0,
        size = 8.2
      )
    )
}

p_cross <- function(indicator_now, centered = FALSE, gs = G_VALUES) {
  d <- plot_dt[
    as.character(indicator) == indicator_now &
      requested_g %in% gs
  ]

  d[
    ,
    g_label := factor(
      paste0("G", requested_g),
      levels = paste0("G", gs)
    )
  ]

  if (!centered) {
    p <- ggplot(
      d,
      aes(
        x = time_hour_center,
        y = fitted_native,
        group = interaction(
          requested_g,
          class
        ),
        colour = factor(class)
      )
    ) +
      geom_line(
        linewidth = 0.58,
        alpha = 0.82,
        show.legend = FALSE
      ) +
      labs(
        title = paste0(
          indicator_now,
          ": G",
          min(gs),
          "-G",
          max(gs),
          " Native-Scale Comparison"
        ),
        y = paste0(
          indicator_now,
          " fitted value"
        )
      )
  } else {
    p <- ggplot(
      d,
      aes(
        x = time_hour_center,
        y = centered_z,
        group = interaction(
          requested_g,
          class
        ),
        colour = factor(class)
      )
    ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.35
      ) +
      geom_line(
        linewidth = 0.58,
        alpha = 0.82,
        show.legend = FALSE
      ) +
      labs(
        title = paste0(
          indicator_now,
          ": G",
          min(gs),
          "-G",
          max(gs),
          " Baseline-Centered Shape Comparison"
        ),
        y = "Change from first point (SD)"
      )
  }

  p +
    facet_wrap(
      ~g_label,
      ncol = if (length(gs) <= 5L) 3 else 4
    ) +
    scale_x_continuous(
      breaks = c(3, 9, 15, 21)
    ) +
    labs(
      subtitle = paste0(
        "Each line is one latent class; all G panels for this indicator ",
        "share the same y-scale."
      ),
      x = "Time window center (hours)"
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      strip.background = element_rect(fill = "grey93"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

log_message("Generating native and centered figures")

pdf(
  file.path(
    figure_dir,
    "V292_native_all_G2G20.pdf"
  ),
  width = 13,
  height = 10,
  onefile = TRUE
)

for (g in G_VALUES) {
  p <- p_native(g)
  print(p)

  ggsave(
    file.path(
      figure_dir,
      paste0(
        "V292_traj_native_G",
        g,
        ".png"
      )
    ),
    p,
    width = 13,
    height = 10,
    units = "in",
    dpi = 320
  )
}

dev.off()

pdf(
  file.path(
    figure_dir,
    "V292_centered_all_G2G20.pdf"
  ),
  width = 13,
  height = 10,
  onefile = TRUE
)

for (g in G_VALUES) {
  p <- p_centered(g)
  print(p)

  ggsave(
    file.path(
      figure_dir,
      paste0(
        "V292_traj_centered_G",
        g,
        ".png"
      )
    ),
    p,
    width = 13,
    height = 10,
    units = "in",
    dpi = 320
  )
}

dev.off()

for (indicator_now in INDICATORS) {
  p_all_native <- p_cross(
    indicator_now,
    centered = FALSE,
    gs = G_VALUES
  )

  p_all_centered <- p_cross(
    indicator_now,
    centered = TRUE,
    gs = G_VALUES
  )

  p_low_native <- p_cross(
    indicator_now,
    centered = FALSE,
    gs = LOW_G_FOCUS
  )

  p_low_centered <- p_cross(
    indicator_now,
    centered = TRUE,
    gs = LOW_G_FOCUS
  )

  ggsave(
    file.path(
      figure_dir,
      paste0(
        "V292_compare_native_",
        indicator_now,
        "_G2G20.png"
      )
    ),
    p_all_native,
    width = 15,
    height = 12,
    units = "in",
    dpi = 320
  )

  ggsave(
    file.path(
      figure_dir,
      paste0(
        "V292_compare_centered_",
        indicator_now,
        "_G2G20.png"
      )
    ),
    p_all_centered,
    width = 15,
    height = 12,
    units = "in",
    dpi = 320
  )

  ggsave(
    file.path(
      figure_dir,
      paste0(
        "V292_compare_native_",
        indicator_now,
        "_G2G6.png"
      )
    ),
    p_low_native,
    width = 15,
    height = 9,
    units = "in",
    dpi = 320
  )

  ggsave(
    file.path(
      figure_dir,
      paste0(
        "V292_compare_centered_",
        indicator_now,
        "_G2G6.png"
      )
    ),
    p_low_centered,
    width = 15,
    height = 9,
    units = "in",
    dpi = 320
  )
}

# =============================================================================
# 9. README, session info, manifest, ZIP
# =============================================================================

readme_lines <- c(
  "FOURTH STUDY / MSI-Electrolyte GBMTM",
  "V292 archived G2-G20 visualization + compact structural extraction",
  "",
  "SOURCE ROOT",
  normalizePath(
    archive_root,
    winslash = "/",
    mustWork = TRUE
  ),
  "",
  "SCOPE",
  "- Reads exact V286 archived copies of completed V251-compatible G2-G20 models.",
  "- No GBMTM refitting.",
  "- No MAKE30 or other outcome.",
  "- Processes one large RDS at a time.",
  "- Per-G compact extraction supports validated resume.",
  "- Native conversion uses exact frozen MIMIC means/SDs embedded in this script.",
  "",
  "KEY FIGURES",
  "- V292_native_all_G2G20.pdf",
  "- V292_centered_all_G2G20.pdf",
  "- V292_compare_*_G2G20.png",
  "- V292_compare_*_G2G6.png",
  "",
  "KEY STEP-3 INPUT TABLES",
  "- V292_trajectory_all_G2G20.csv",
  "- V292_assignment_all_G2G20.csv",
  "- V292_class_quality_all_G2G20.csv",
  "- V292_model_summary_all_G2G20.csv",
  "- V292_level_shape_by_indicator_G2G20.csv",
  "- V292_level_shape_joint_G2G20.csv",
  "- V292_nearest_shape_pair_G2G20.csv",
  "- V292_transition_G2G20.csv",
  "- V292_transition_metrics_G2G20.csv",
  "",
  "NOTE",
  "Baseline-centered plots subtract the first fitted point for structural comparison.",
  "Electrolyte clinical interpretation remains restricted to the real 3/9/15/21 h anchors."
)

writeLines(
  readme_lines,
  file.path(
    output_dir,
    "V292_README_G2G20.txt"
  ),
  useBytes = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "V292_sessionInfo_G2G20.txt"
  ),
  useBytes = TRUE
)

manifest_path <- file.path(
  table_dir,
  "V292_output_manifest_G2G20.csv"
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
    basename(manifest_files) != basename(zip_path) &
    normalizePath(
      manifest_files,
      winslash = "/",
      mustWork = TRUE
    ) != normalizePath(
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
    file.info(manifest_files)$size
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
  unlink(zip_path, force = TRUE)
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
    "Final ZIP creation failed: ",
    zip_path,
    call. = FALSE
  )
}

log_message(
  "V292 completed successfully."
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
cat("V292 completed successfully.\n")
cat("Output directory:\n")
cat(
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)
cat("Uploadable ZIP:\n")
cat(
  normalizePath(
    zip_path,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)
cat("============================================================\n")
