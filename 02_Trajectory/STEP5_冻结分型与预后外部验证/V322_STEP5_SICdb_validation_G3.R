# V322: SICdb external validation of the frozen V295 G3 S01 phenotype.
# 在 RStudio 打开本文件，Ctrl+A 全选，然后 Ctrl+Enter；不弹出选文件窗口。
# 输入和输出均使用当前用户 Desktop/Research/MSI/提取轨迹文件/GBMTM_stage5_output。
# 完整独立脚本。固定使用已核验的 V318 分类，不重新拟合或改变类别标签。
# eICU 的两条 Unknown/Other 性别沿用 V317：MICE 20 次、5 轮、种子 20260927。
# 普通协变量插补不改变任何轨迹观测值、冻结后验概率或原始 CSV/RDS。

main_external_validation <- function() {
  SELF <- sys.function()
  old_options <- options(stringsAsFactors = FALSE, scipen = 999, error = NULL,
                         contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(old_options), add = TRUE)
  VERSION <- "V322"
  DB <- "SICdb"
  ID <- "case_id"
  PATIENT_ID <- "patient_id"
  RAW_NAME <- "SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv"
  RAW_MD5 <- "d52a3ecdbb636bba2859eb9f102fdf45"
  CLASS_NAME <- "V318_SICdb_fixed_scoring_G3.csv"
  CLASS_MD5 <- "c3ea4af8b771bbe91e7db7fca41e70e4"
  EXPECTED_N <- 3831L
  EXPECTED_CLASS_COUNTS <- c(1125L, 951L, 1755L)
  CLASSES <- c("C1", "C2", "C3")
  CORE <- c("age", "sex", "sepsis_sofa_score", "baseline_egfr_mdrd",
            "ckd_or_renal_disease", "diabetes",
            "mechanical_ventilation_0_24h", "vasopressor_0_24h")
  Z <- c("msi_z", "sodium_z", "potassium_z", "chloride_z", "bicarbonate_z")
  INDICATORS <- c("MSI", "Na", "K", "Cl", "HCO3")
  RAW_GROUPS <- list(
    MSI = paste0("msi_m", 1:12), Na = paste0("sodium_t", 1:4),
    K = paste0("potassium_t", 1:4), Cl = paste0("chloride_t", 1:4),
    HCO3 = paste0("bicarbonate_t", 1:4))
  HOURS <- list(MSI = seq(1, 23, 2), Na = c(3, 9, 15, 21),
                K = c(3, 9, 15, 21), Cl = c(3, 9, 15, 21),
                HCO3 = c(3, 9, 15, 21))
  MI_M <- 20L
  MI_MAXIT <- 5L
  MI_SEED <- 20260927L
  CLASS_COLORS <- c(C1 = "#287C8E", C2 = "#C5543E", C3 = "#6B68A0")
  MARKER_COLORS <- c(MSI = "#287C8E", Na = "#CD653B", K = "#7D68A2",
                     Cl = "#659140", HCO3 = "#BF4E84")

  profile <- Sys.getenv("USERPROFILE", unset = "")
  if (!nzchar(profile) || !dir.exists(profile)) {
    profile <- path.expand("~")
    if (identical(basename(profile), "Documents"))
      stop("USERPROFILE unavailable and ~ resolves to Documents; correct USERPROFILE first.")
  }
  root <- file.path(profile, "Desktop", "Research", "MSI",
                    "提取轨迹文件", "GBMTM_stage5_output")
  if (!dir.exists(root)) stop("Exact STEP5 root does not exist: ", root)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  kit <- file.path(root, "V307_STEP5_laptop_kit_G3")
  scoring <- file.path(root, "V318_STEP5_scoring_G3_89f684d006_f1eb2d7f")
  out_dir <- file.path(root, paste0(VERSION, "_STEP5_", DB, "_validation_G3"))
  dir.create(out_dir, recursive = FALSE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Cannot create output directory: ", out_dir)
  filename <- function(content, extension = "csv")
    file.path(out_dir, paste0(VERSION, "_", content, "_G3.", extension))
  progress_log <- function(...) {
    line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                   paste0(..., collapse = ""))
    cat(line, "\n", sep = "")
    cat(line, "\n", sep = "", file = filename("validation_log", "txt"), append = TRUE)
  }
  replace_generated <- function(tmp, dest) {
    if (file.exists(dest) && unlink(dest) != 0L)
      stop("Cannot replace generated output: ", dest)
    if (!file.rename(tmp, dest)) stop("Cannot save generated output: ", dest)
    invisible(dest)
  }
  write_csv <- function(d, content) {
    dest <- filename(content)
    tmp <- paste0(dest, ".tmp_", Sys.getpid())
    oo <- options(scipen = 0)
    on.exit(options(oo), add = TRUE)
    utils::write.csv(d, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
    replace_generated(tmp, dest)
  }
  write_text <- function(lines, content, extension = "txt") {
    dest <- filename(content, extension)
    tmp <- paste0(dest, ".tmp_", Sys.getpid())
    writeLines(enc2utf8(lines), tmp, useBytes = TRUE)
    replace_generated(tmp, dest)
  }
  save_checkpoint <- function(x, content) {
    dest <- filename(content, "rds")
    tmp <- paste0(dest, ".tmp_", Sys.getpid())
    saveRDS(x, tmp, compress = TRUE)
    replace_generated(tmp, dest)
  }
  read_csv <- function(path, id = NULL) {
    utils::read.csv(path, check.names = FALSE,
      colClasses = if (is.null(id)) NA else stats::setNames(rep("character", length(id)), id),
      na.strings = c("", "NA", "NULL"))
  }
  assert_columns <- function(d, required, label) {
    if (anyDuplicated(names(d)) || !all(required %in% names(d)))
      stop(label, " exact schema differs; missing: ", paste(setdiff(required, names(d)), collapse = ", "))
  }
  check_binary <- function(x, label, allow_missing = TRUE) {
    if (all(is.na(x)) && allow_missing) return(invisible(TRUE))
    if (!is.numeric(x) || any(!is.na(x) & !(x %in% c(0, 1))) ||
        (!allow_missing && anyNA(x))) stop("Invalid binary coding: ", label)
  }
  check_numeric <- function(x, label, allow_missing = TRUE) {
    if ((!is.numeric(x) && !all(is.na(x))) ||
        any(!is.finite(x[!is.na(x)])) || (!allow_missing && anyNA(x)))
      stop("Invalid numeric values: ", label)
  }
  wilson <- function(events, n) {
    if (!n) return(c(NA_real_, NA_real_))
    z <- stats::qnorm(0.975); p <- events / n
    center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    c(center - half, center + half)
  }
  read_classes <- function(path, raw, id, n, expected_counts = NULL, rounded = FALSE) {
    cl <- read_csv(path, "unit_id")
    fields <- c("unit_id", "assigned_class", "p_class1", "p_class2", "p_class3")
    assert_columns(cl, fields, basename(path))
    if (nrow(raw) != n || nrow(cl) != n || anyNA(raw[[id]]) ||
        anyNA(cl$unit_id) || anyDuplicated(raw[[id]]) || anyDuplicated(cl$unit_id) ||
        !setequal(as.character(raw[[id]]), cl$unit_id) ||
        anyNA(cl$assigned_class) || !all(cl$assigned_class %in% 1:3))
      stop("Exact cohort IDs/classifications differ: ", basename(path))
    p <- as.matrix(cl[c("p_class1", "p_class2", "p_class3")])
    tolerance <- if (rounded) 0.000151 else 1e-10
    if (!is.numeric(p) || anyNA(p) || any(!is.finite(p) | p < 0 | p > 1) ||
        max(abs(rowSums(p) - 1)) > tolerance)
      stop("Invalid stored posterior: ", basename(path))
    if (!rounded && any(cl$assigned_class != max.col(p, ties.method = "first")))
      stop("Stored class is inconsistent with the maximum posterior.")
    ix <- match(as.character(raw[[id]]), cl$unit_id)
    raw$assigned_class <- cl$assigned_class[ix]
    raw$phenotype <- factor(paste0("C", raw$assigned_class), levels = CLASSES)
    raw$max_posterior <- apply(p, 1L, max)[ix]
    for (v in c("p_class1", "p_class2", "p_class3")) raw[[v]] <- cl[[v]][ix]
    if (!is.null(expected_counts) &&
        !identical(as.integer(table(raw$phenotype)), as.integer(expected_counts)))
      stop("Frozen class counts differ from the inspected archive.")
    list(data = raw, classes = cl)
  }
  observed_summary <- function(d, db) {
    rows <- list(); j <- 0L
    for (k in CLASSES) for (indicator in INDICATORS) {
      sub <- d[d$phenotype == k, , drop = FALSE]
      columns <- RAW_GROUPS[[indicator]]
      variable <- Z[match(indicator, INDICATORS)]
      for (t in seq_along(columns)) {
        x <- sub[[columns[t]]]; known <- x[!is.na(x)]; n <- length(known)
        q <- if (n) stats::quantile(known, c(0.25, 0.5, 0.75), names = FALSE) else rep(NA_real_, 3L)
        j <- j + 1L
        rows[[j]] <- data.frame(dataset = db, phenotype = k, indicator = indicator,
          variable = variable, source_column = columns[t], time_hour = HOURS[[indicator]][t],
          class_n = nrow(sub), observed_n = n, missing_n = sum(is.na(x)),
          missing_pct = 100 * sum(is.na(x)) / nrow(sub),
          mean_native = if (n) mean(known) else NA_real_,
          sd_native = if (n > 1L) stats::sd(known) else NA_real_,
          Q25_native = q[1L], median_native = q[2L], Q75_native = q[3L])
      }
    }
    ans <- do.call(rbind, rows)
    ans$mean_z <- (ans$mean_native - state$raw_center[ans$variable]) / state$raw_scale[ans$variable]
    ans$baseline_centered_z <- NA_real_
    for (k in CLASSES) for (indicator in INDICATORS) {
      ii <- which(ans$phenotype == k & ans$indicator == indicator)
      ans$baseline_centered_z[ii] <- ans$mean_z[ii] - ans$mean_z[ii[1L]]
    }
    ans
  }
  endpoint_rates <- function(d, db) {
    endpoints <- c(MAKE30 = "make30", Death = "make30_death", New_KRT = "make30_new_krt",
      Persistent_kidney_dysfunction = "make30_persistent_kidney_dysfunction")
    do.call(rbind, lapply(names(endpoints), function(endpoint) {
      v <- endpoints[[endpoint]]
      do.call(rbind, lapply(CLASSES, function(k) {
        sub <- d[d$phenotype == k & d$make30_evaluable_window == 1L, , drop = FALSE]
        x <- sub[[v]][!is.na(sub[[v]])]; n <- length(x); events <- sum(x); ci <- wilson(events, n)
        data.frame(dataset = db, endpoint = endpoint, phenotype = k,
          eligible_n = nrow(sub), known_outcome_n = n, missing_outcome_n = sum(is.na(sub[[v]])),
          events = events, event_pct = if (n) 100 * events / n else NA_real_,
          CI95_low_pct = 100 * ci[1L], CI95_high_pct = 100 * ci[2L])
      }))
    }))
  }
  baseline_table <- function(pool) {
    continuous <- c("age", "weight_kg", "height_cm", "bmi", "sepsis_sofa_score",
                    "baseline_creatinine", "baseline_egfr_mdrd")
    if (DB == "SICdb") continuous <- c(continuous, "saps3")
    binary <- c("ckd_or_renal_disease", "diabetes", "congestive_heart_failure",
                "liver_disease", "malignancy", "mechanical_ventilation_0_24h", "vasopressor_0_24h")
    if (DB == "eICU") binary <- c(binary, "immunosuppression")
    categorical <- c("sex", "race_ethnicity", "baseline_creatinine_source")
    assert_columns(pool, c(continuous, binary, categorical), paste0(DB, " Table1"))
    cell <- function(x, kind, level) {
      known <- x[!is.na(x)]; n <- length(known)
      if (!n) return("not available")
      if (kind == "continuous") {
        q <- stats::quantile(known, c(0.25, 0.5, 0.75), names = FALSE)
        return(sprintf("%.2f [%.2f, %.2f]", q[2L], q[1L], q[3L]))
      }
      e <- if (kind == "binary") sum(known == 1) else sum(as.character(known) == level)
      sprintf("%d/%d (%.1f%%)", e, n, 100 * e / n)
    }
    rows <- list(); j <- 0L
    add_row <- function(v, kind, level = "") {
      row <- data.frame(dataset = DB, variable = v,
        statistic = if (kind == "continuous") "median [IQR]" else "n/known (%)",
        level = level, missing_n = sum(is.na(pool[[v]])),
        status = if (all(is.na(pool[[v]]))) "unavailable" else "observed")
      row$Overall <- cell(pool[[v]], kind, level)
      for (k in CLASSES) {
        sub <- pool[[v]][pool$phenotype == k]
        row[[k]] <- cell(sub, kind, level)
        row[[paste0(k, "_missing_n")]] <- sum(is.na(sub))
      }
      row
    }
    for (v in continuous) { j <- j + 1L; rows[[j]] <- add_row(v, "continuous") }
    for (v in binary) {
      check_binary(pool[[v]], v)
      j <- j + 1L; rows[[j]] <- add_row(v, "binary", "1")
    }
    for (v in categorical) {
      levels <- sort(unique(as.character(pool[[v]][!is.na(pool[[v]])])))
      if (!length(levels)) levels <- ""
      for (level in levels) { j <- j + 1L; rows[[j]] <- add_row(v, "categorical", level) }
    }
    do.call(rbind, rows)
  }

    fit_one <- function(data, estimator, adjusted, model_name, imputation = 0L) {
    rhs <- c("phenotype", if (adjusted) CORE else character(0L))
    formula <- stats::as.formula(paste("make30 ~", paste(rhs, collapse = " + ")), env = baseenv())
    formula_text <- paste(deparse(formula), collapse = "")
    model_signature <- list(base = signature, estimator = estimator, model_name = model_name,
      formula = formula_text, imputation = imputation, n = nrow(data), events = sum(data$make30),
      binary_sex_at_missing_records = as.character(data$sex[match(a[[ID]][is.na(a$sex)], data[[ID]])]),
      mice_version = if (imputation > 0L) as.character(utils::packageVersion("mice")) else NULL)
    content <- paste0("model_", estimator, "_", model_name,
                      if (imputation > 0L) sprintf("_MI%02d", imputation) else "")
    dest <- filename(content, "rds")
    if (file.exists(dest)) {
      saved <- readRDS(dest)
      if (!identical(saved$signature, model_signature) || !isTRUE(saved$converged) ||
          !identical(saved$formula, formula_text))
        stop("Existing model checkpoint is incompatible: ", dest,
             "\nKeep original inputs; move only this script's output folder aside before a new run.")
      progress_log("Resumed ", estimator, " / ", model_name,
                   if (imputation > 0L) paste0(" / MI", imputation) else "", ".")
      return(saved)
    }
    warnings <- character(0L)
    fit <- withCallingHandlers(stats::glm(formula, data = data,
      family = if (estimator == "Poisson_RR") stats::poisson(link = "log") else stats::binomial(link = "logit"),
      na.action = stats::na.fail, x = TRUE, y = TRUE), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    b <- stats::coef(fit)
    if (!isTRUE(fit$converged) || anyNA(b) || any(!is.finite(b)) ||
        stats::nobs(fit) != nrow(data) || fit$rank != ncol(fit$x))
      stop("Regression convergence/rank failed: ", estimator, " / ", model_name)
    covariance <- stats::vcov(fit)
    if (estimator == "Poisson_RR") {
      # Identical to the reviewed V319/V320 MIMIC sandwich covariance.
      scores <- fit$x * as.numeric(fit$y - stats::fitted(fit))
      covariance <- covariance %*% crossprod(scores) %*% covariance *
        nrow(fit$x) / (nrow(fit$x) - ncol(fit$x))
    }
    if (any(!is.finite(covariance)) || any(diag(covariance) <= 0))
      stop("Invalid regression covariance: ", estimator, " / ", model_name)
    result <- list(signature = model_signature, estimator = estimator, model_name = model_name,
      formula = formula_text, coefficients = b, covariance = covariance,
      df = stats::setNames(rep(Inf, length(b)), names(b)), n = nrow(data),
      events = sum(data$make30), npar = ncol(fit$x), converged = fit$converged,
      imputation = imputation, imputations = 0L, pooling = "none",
      fitted_min = min(stats::fitted(fit)), fitted_max = max(stats::fitted(fit)),
      n_fitted_above_one = sum(stats::fitted(fit) > 1), warnings = warnings)
    save_checkpoint(result, content)
    progress_log("Finished ", estimator, " / ", model_name,
                 if (imputation > 0L) paste0(" / MI", imputation) else "", ".")
    result
  }
  make_imputations <- function() {
    if (!requireNamespace("mice", quietly = TRUE)) {
      progress_log("Installing mice for the two eICU binary-sex covariates; trajectories remain untouched.")
      oo <- options(timeout = max(300, getOption("timeout")))
      tryCatch(utils::install.packages("mice", repos = "https://cloud.r-project.org"),
               finally = options(oo))
    }
    if (!requireNamespace("mice", quietly = TRUE))
      stop("mice could not be installed. Run install.packages('mice'), then rerun this entire script.")
    variables <- unique(c("sex", setdiff(CORE, "sex"), "phenotype", "make30"))
    input <- a[variables]
    methods <- stats::setNames(rep("", length(variables)), variables)
    methods["sex"] <- "logreg"
    predictors <- matrix(0L, length(variables), length(variables), dimnames = list(variables, variables))
    predictors["sex", setdiff(variables, "sex")] <- 1L
    mi_signature <- list(base = signature, variables = variables, method = methods,
      predictorMatrix = predictors, m = MI_M, maxit = MI_MAXIT, seed = MI_SEED,
      RNGkind = c("Mersenne-Twister", "Inversion", "Rejection"),
      mice_version = as.character(utils::packageVersion("mice")))
    dest <- filename("sex_MICE_checkpoint", "rds")
    if (file.exists(dest)) {
      saved <- readRDS(dest)
      if (!identical(saved$signature, mi_signature) || !inherits(saved$mids, "mids") ||
          saved$mids$m != MI_M)
        stop("Existing MICE checkpoint is incompatible: ", dest)
      imp <- saved$mids
      progress_log("Resumed the completed 20-imputation binary-sex checkpoint.")
    } else {
      previous_rng <- RNGkind()
      had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      previous_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
      on.exit({
        do.call(RNGkind, as.list(previous_rng))
        if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv)
      }, add = TRUE)
      RNGkind("Mersenne-Twister", "Inversion", "Rejection")
      # Retain V317's exact method, predictors, m/maxit/seed. MAKE30 only predicts
      # an ordinary covariate after all phenotype assignments have been fixed.
      imp <- mice::mice(input, m = MI_M, maxit = MI_MAXIT, method = methods,
                        predictorMatrix = predictors, seed = MI_SEED, printFlag = FALSE)
      save_checkpoint(list(signature = mi_signature, mids = imp), "sex_MICE_checkpoint")
      progress_log("Completed MICE: 20 imputations, 5 iterations; seed ", MI_SEED, ".")
    }
    write_csv(data.frame(variable = variables, method = unname(methods),
      predictor_names = vapply(seq_along(variables), function(i)
        paste(variables[predictors[i, ] == 1L], collapse = " | "), character(1L))), "MICE_methods")
    events <- imp$loggedEvents
    if (is.null(events)) events <- data.frame(note = "No logged MICE events.")
    write_csv(events, "MICE_logged_events")
    missing_ix <- which(is.na(a$sex))
    datasets <- vector("list", MI_M); record <- vector("list", MI_M)
    for (j in seq_len(MI_M)) {
      completed <- mice::complete(imp, action = j)
      z <- a
      z$sex <- factor(as.character(completed$sex), levels = c("Female", "Male"))
      if (anyNA(z$sex)) stop("Incomplete binary-sex imputation: ", j)
      datasets[[j]] <- z
      record[[j]] <- data.frame(imputation = j, unit_id = as.character(a[[ID]][missing_ix]),
        reported_sex = a$sex_reported[missing_ix],
        imputed_binary_sex = as.character(z$sex[missing_ix]))
    }
    write_csv(do.call(rbind, record), "sex_MICE_20x")
    datasets
  }
  pool_models <- function(fits, model_name) {
    terms <- names(fits[[1L]]$coefficients)
    if (!all(vapply(fits, function(x) identical(names(x$coefficients), terms), logical(1L))))
      stop("Coefficient terms differ among imputations.")
    Q <- do.call(rbind, lapply(fits, function(x) x$coefficients))
    ubar <- Reduce("+", lapply(fits, function(x) x$covariance)) / length(fits)
    between <- stats::cov(Q)
    total <- ubar + (1 + 1 / length(fits)) * between
    scalar <- lapply(seq_along(terms), function(j)
      mice::pool.scalar(Q = Q[, j], U = vapply(fits, function(x) x$covariance[j, j], numeric(1L)),
        n = fits[[1L]]$n, k = fits[[1L]]$npar, rule = "rubin1987"))
    b <- stats::setNames(vapply(scalar, function(x) x$qbar, numeric(1L)), terms)
    dfs <- stats::setNames(vapply(scalar, function(x) x$df, numeric(1L)), terms)
    scalar_variance <- vapply(scalar, function(x) x$t, numeric(1L))
    if (max(abs(diag(total) - scalar_variance)) > 1e-10 ||
        any(!is.finite(total)) || any(diag(total) <= 0) || anyNA(dfs) || any(dfs <= 0))
      stop("Rubin-pooled covariance or degrees of freedom are invalid.")
    result <- list(signature = signature, estimator = fits[[1L]]$estimator,
      model_name = model_name, formula = fits[[1L]]$formula, coefficients = b,
      covariance = total, df = dfs, n = fits[[1L]]$n, events = fits[[1L]]$events,
      npar = fits[[1L]]$npar, converged = all(vapply(fits, function(x) x$converged, logical(1L))),
      imputation = 0L, imputations = length(fits), pooling = "Rubin1987; Barnard-Rubin df",
      fitted_min = min(vapply(fits, function(x) x$fitted_min, numeric(1L))),
      fitted_max = max(vapply(fits, function(x) x$fitted_max, numeric(1L))),
      n_fitted_above_one = NA_integer_, warnings = unique(unlist(lapply(fits, function(x) x$warnings))))
    save_checkpoint(result, paste0("pooled_", result$estimator, "_", model_name))
    result
  }
  contrast_rows <- function(z) {
    terms <- c("phenotypeC2", "phenotypeC3")
    if (!all(terms %in% names(z$coefficients))) stop("Expected frozen C2/C3 contrasts unavailable.")
    b <- z$coefficients[terms]; se <- sqrt(diag(z$covariance))[terms]; df <- z$df[terms]
    critical <- stats::qt(0.975, df = df)
    data.frame(dataset = DB, estimator = z$estimator, model = z$model_name,
      comparison = c("C2 vs C1", "C3 vs C1"), estimate = exp(b),
      CI95_low = exp(b - critical * se), CI95_high = exp(b + critical * se),
      SE_log_effect = se, p_value = 2 * stats::pt(-abs(b / se), df = df), df = df,
      n = z$n, events = z$events, imputations = z$imputations, pooling = z$pooling, row.names = NULL)
  }
    plot_save <- function(content, draw, width = 10, height = 7) {
    grDevices::pdf(filename(content, "pdf"), width = width, height = height, useDingbats = FALSE)
    tryCatch(draw(), finally = grDevices::dev.off())
    if (isTRUE(capabilities("png"))) {
      args <- list(filename = filename(content, "png"), width = round(width * 300),
                   height = round(height * 300), res = 300)
      if (isTRUE(capabilities("cairo"))) args$type <- "cairo"
      do.call(grDevices::png, args)
      tryCatch(draw(), finally = grDevices::dev.off())
    }
  }
  padded_range <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (any(!is.finite(r))) return(c(-1, 1))
    margin <- if (diff(r) > 0) 0.12 * diff(r) else max(0.1, abs(r[1L]) * 0.02)
    r + c(-margin, margin)
  }
  draw_native <- function() {
    graphics::par(mfrow = c(2, 3), mar = c(4.0, 4.5, 3.0, 0.7), oma = c(0, 0, 2, 0))
    for (indicator in INDICATORS) {
      s <- observed[observed$dataset == DB & observed$indicator == indicator, ]
      graphics::plot(NA, NA, xlim = c(0, 24), ylim = padded_range(s$mean_native),
        xlab = "Hours after sepsis onset", ylab = if (indicator == "MSI") "MSI" else "mmol/L",
        main = indicator, xaxt = "n")
      graphics::axis(1, at = c(0, 6, 12, 18, 24))
      graphics::grid(col = "#E5E5E5")
      for (k in CLASSES) {
        p <- s[s$phenotype == k, ]
        graphics::lines(p$time_hour, p$mean_native, type = "b", pch = 16,
                        col = CLASS_COLORS[k], lwd = 2, cex = 0.7)
      }
    }
    graphics::plot.new()
    graphics::legend("center", legend = paste0(CLASSES, ": n=", EXPECTED_CLASS_COUNTS),
                      col = CLASS_COLORS, pch = 16, lwd = 2, bty = "n", cex = 1.1)
    graphics::mtext(paste0(DB, ": observed native-scale trajectories; frozen G3 labels"),
                     outer = TRUE, line = 0.5, font = 2)
  }
  draw_joint <- function(centered = FALSE) {
    field <- if (centered) "baseline_centered_z" else "mean_z"
    s <- observed[observed$dataset == DB, ]
    yl <- padded_range(s[[field]])
    graphics::par(mfrow = c(1, 3), mar = c(4.2, 4.5, 3.0, 0.8), oma = c(3.0, 0, 2.0, 0))
    for (k in CLASSES) {
      p <- s[s$phenotype == k, ]
      graphics::plot(NA, NA, xlim = c(0, 24), ylim = yl, xaxt = "n",
        xlab = "Hours after sepsis onset",
        ylab = if (centered) "Change from first observed mean, z" else "Observed mean, MIMIC-standardized z",
        main = paste0(k, " (n=", EXPECTED_CLASS_COUNTS[match(k, CLASSES)], ")"))
      graphics::axis(1, at = c(0, 6, 12, 18, 24)); graphics::grid(col = "#E5E5E5")
      graphics::abline(h = 0, lty = 3, col = "#999999")
      for (j in seq_along(INDICATORS)) {
        indicator <- INDICATORS[j]; z <- p[p$indicator == indicator, ]
        graphics::lines(z$time_hour, z[[field]], type = "b", pch = c(16, 17, 15, 18, 1)[j],
                        col = MARKER_COLORS[indicator], lwd = 2, cex = 0.65)
      }
      if (k == "C1") graphics::legend("bottomleft", legend = INDICATORS,
        col = MARKER_COLORS, pch = c(16, 17, 15, 18, 1), lwd = 2, bty = "n", cex = 0.75)
    }
    graphics::mtext(paste0(DB, ": five indicators together within each frozen phenotype"),
                     outer = TRUE, line = 0.5, font = 2)
    graphics::mtext(if (centered) "Each class/indicator starts at zero; this graph emphasizes temporal shape."
      else "All z values use the frozen MIMIC center and SD; no external restandardization.",
      side = 1, outer = TRUE, line = 1, cex = 0.8)
  }
  draw_comparison <- function(centered = FALSE) {
    field <- if (centered) "baseline_centered_z" else "mean_native"
    graphics::par(mfrow = c(5, 3), mar = c(3.3, 4.1, 2.4, 0.7), oma = c(2.6, 0, 2.2, 0))
    for (indicator in INDICATORS) {
      s <- observed[observed$indicator == indicator, ]
      yl <- padded_range(s[[field]])
      for (k in CLASSES) {
        p <- s[s$phenotype == k, ]
        graphics::plot(NA, NA, xlim = c(0, 24), ylim = yl, xaxt = "n",
          xlab = "Hours after sepsis onset",
          ylab = paste0(indicator, if (centered) " (change, z)" else if (indicator == "MSI") "" else " (mmol/L)"),
          main = paste0(k, " / ", indicator))
        graphics::axis(1, at = c(0, 6, 12, 18, 24)); graphics::grid(col = "#E5E5E5")
        if (centered) graphics::abline(h = 0, lty = 3, col = "#BBBBBB")
        ref <- p[p$dataset == "MIMIC", ]; ext <- p[p$dataset == DB, ]
        graphics::lines(ref$time_hour, ref[[field]], type = "b", lty = 2, pch = 1,
                        col = "#777777", lwd = 1.6, cex = 0.6)
        graphics::lines(ext$time_hour, ext[[field]], type = "b", lty = 1, pch = 16,
                        col = CLASS_COLORS[k], lwd = 2, cex = 0.6)
        if (indicator == "MSI" && k == "C1") graphics::legend("topright",
          legend = c("MIMIC", DB), col = c("#777777", CLASS_COLORS[k]),
          lty = c(2, 1), pch = c(1, 16), lwd = c(1.6, 2), bty = "n", cex = 0.7)
      }
    }
    graphics::mtext(paste0("MIMIC versus ", DB, ": ",
      if (centered) "baseline-centered observed shapes" else "observed native-scale trajectories"),
      outer = TRUE, line = 0.6, font = 2)
    graphics::mtext("Dashed gray: MIMIC; solid class color: external cohort. Shared scale within each indicator.",
      side = 1, outer = TRUE, line = 1, cex = 0.75)
  }
  draw_rates <- function(components = FALSE) {
    endpoints <- if (components) c("MAKE30", "Death", "New_KRT", "Persistent_kidney_dysfunction") else "MAKE30"
    graphics::par(mfrow = if (components) c(2, 2) else c(1, 1),
                  mar = c(5.3, 4.5, 3.0, 0.8))
    for (endpoint in endpoints) {
      s <- rates[rates$endpoint == endpoint, ]
      loc <- graphics::barplot(s$event_pct, names.arg = s$phenotype, col = CLASS_COLORS,
        border = NA, ylim = c(0, max(s$CI95_high_pct, na.rm = TRUE) * 1.25),
        ylab = "Observed events, %", main = paste0(DB, ": ", gsub("_", " ", endpoint, fixed = TRUE)))
      graphics::arrows(loc, s$CI95_low_pct, loc, s$CI95_high_pct, angle = 90, code = 3, length = 0.05)
      graphics::text(loc, s$CI95_high_pct + max(s$CI95_high_pct) * 0.07,
        sprintf("%.1f%%", s$event_pct), cex = 0.9)
      graphics::mtext(paste(sprintf("%s: %d/%d", s$phenotype, s$events, s$known_outcome_n),
                             collapse = "    "), side = 1, line = 3.0, cex = 0.72)
    }
  }
  draw_rr <- function() {
    s <- rr_comparison
    yy <- rev(seq_len(nrow(s)))
    limits <- c(min(0.7, s$CI95_low / 1.3), max(5, s$CI95_high * 1.3))
    graphics::par(mar = c(5.2, 9.0, 3.5, 1.0))
    graphics::plot(NA, NA, xlim = limits, ylim = c(0.5, nrow(s) + 0.8), log = "x",
      xlab = "Adjusted relative risk (95% CI); C1 reference", ylab = "",
      yaxt = "n", xaxt = "n", bty = "n", main = paste0("MAKE30 association: MIMIC and ", DB))
    ticks <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 8, 10, 16)
    ticks <- ticks[ticks >= limits[1L] & ticks <= limits[2L]]
    graphics::axis(1, at = ticks, labels = ticks)
    graphics::axis(2, at = yy, labels = paste0(s$dataset, ": ", s$comparison), las = 1, tick = FALSE, cex.axis = 0.85)
    graphics::abline(v = 1, lty = 2, col = "#999999")
    colors <- ifelse(s$comparison == "C2 vs C1", CLASS_COLORS["C2"], CLASS_COLORS["C3"])
    graphics::segments(s$CI95_low, yy, s$CI95_high, yy, col = colors, lwd = 2)
    graphics::points(s$RR, yy, col = colors, pch = ifelse(s$dataset == "MIMIC", 1, 15), cex = 1.2)
    graphics::text(s$RR, yy + 0.26, labels = sprintf("%.2f (%.2f-%.2f)", s$RR, s$CI95_low, s$CI95_high), cex = 0.85)
    graphics::mtext("Same eight adjustment variables; separate cohort associations, not a transported risk-prediction model.",
      side = 1, line = 3.5, cex = 0.7)
  }
  archive_results <- function() {
    target <- paste0(out_dir, ".zip")
    stage <- paste0(target, ".tmp_", Sys.getpid(), ".zip")
    if (file.exists(stage)) unlink(stage)
    if (.Platform$OS.type == "windows") {
      ps <- Sys.which("powershell.exe")
      if (!nzchar(ps)) stop("PowerShell was not found; results remain in: ", out_dir)
      quote_ps <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
      ps_file <- tempfile(paste0(VERSION, "_zip_"), fileext = ".ps1")
      on.exit(unlink(ps_file), add = TRUE)
      lines <- c("$ErrorActionPreference = 'Stop'",
        paste0("$src = ", quote_ps(out_dir)), paste0("$dst = ", quote_ps(stage)),
        "$items = @(Get-ChildItem -LiteralPath $src -File | Where-Object { $_.Name -notmatch '\\.tmp_[0-9]+$' } | Select-Object -ExpandProperty FullName)",
        "Compress-Archive -LiteralPath $items -DestinationPath $dst -CompressionLevel Optimal -Force")
      con <- file(ps_file, "wb")
      writeBin(c(as.raw(c(239, 187, 191)), charToRaw(enc2utf8(paste(lines, collapse = "\r\n")))), con)
      close(con)
      code <- system2(ps, c("-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass", "-File", shQuote(ps_file)))
      if (!identical(as.integer(code), 0L) || !file.exists(stage)) stop("ZIP creation failed; directory: ", out_dir)
    } else {
      zipper <- Sys.which("zip")
      if (!nzchar(zipper)) stop("zip executable not found; directory: ", out_dir)
      oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE)
      setwd(out_dir)
      entries <- list.files(".", full.names = FALSE)
      entries <- entries[!grepl("\\.tmp_[0-9]+$", entries)]
      utils::zip(stage, files = entries, flags = "-q", zip = zipper)
      setwd(oldwd)
      if (!file.exists(stage)) stop("ZIP creation failed; directory: ", out_dir)
    }
    replace_generated(stage, target)
    target
  }

  progress_log(DB, " external validation started; original V318 classes fixed.")
  status <- tryCatch({
    inputs <- c(
      external_raw = file.path(kit, RAW_NAME),
      external_classification = file.path(scoring, CLASS_NAME),
      frozen_scoring_state = file.path(kit, "V318_frozen_scoring_state_G3.rds"),
      MIMIC_raw_reference = file.path(kit, "MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv"),
      MIMIC_formal_classification = file.path(scoring, "V318_MIMIC_formal_classification_G3.csv"),
      MIMIC_replay_audit = file.path(scoring, "V318_MIMIC_replay_audit_G3.csv"),
      MIMIC_RR_reference = file.path(root, "V319_STEP5_MIMIC_MAKE30_G3",
                                     "V319_MAKE30_modified_Poisson_RR_G3.csv"))
    if (any(!file.exists(inputs))) stop("Missing exact inputs:\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
    inputs <- vapply(inputs, normalizePath, character(1L), winslash = "/", mustWork = TRUE)
    input_hashes <- stats::setNames(unname(tools::md5sum(inputs)), names(inputs))
    expected <- c(external_raw = RAW_MD5, external_classification = CLASS_MD5,
      frozen_scoring_state = "89f684d00680fafc30577d677fee6a30",
      MIMIC_raw_reference = "f1eb2d7fac1303a682b53c679dde9abd",
      MIMIC_formal_classification = "d6824b387dcc2f3cfd5c72cb82d69e94")
    if (!identical(unname(input_hashes[names(expected)]), unname(expected)))
      stop("An exact frozen input differs: ", paste(names(expected)[input_hashes[names(expected)] != expected], collapse = ", "))
    state <- readRDS(inputs["frozen_scoring_state"])
    required_state <- c("schema", "G", "degree", "variables", "source_fit_md5", "source_MIMIC_md5",
      "completion_parameters_iteration", "score_parameters_iteration", "gbmt_version",
      "raw_center", "raw_scale", "archived_class_agreement", "archived_posterior_4dp_max_abs_diff")
    if (!all(required_state %in% names(state)) || state$schema != "V318_recovered_V295_G3_S01" ||
        state$G != 3L || state$degree != 2L || !identical(state$variables, Z) ||
        state$source_fit_md5 != "cf5881479ee206195aa80d59fa3bd17f" ||
        state$source_MIMIC_md5 != expected["MIMIC_raw_reference"] ||
        state$completion_parameters_iteration != 162L || state$score_parameters_iteration != 163L ||
        state$gbmt_version != "0.1.4" || state$archived_class_agreement != 1 ||
        state$archived_posterior_4dp_max_abs_diff != 0 ||
        !identical(names(state$raw_center), Z) || !identical(names(state$raw_scale), Z) ||
        any(!is.finite(state$raw_center)) || any(!is.finite(state$raw_scale) | state$raw_scale <= 0))
      stop("Frozen G3 state identity/standardization differs.")
    audit <- read_csv(inputs["MIMIC_replay_audit"])
    assert_columns(audit, c("n_patients", "n_class_changed", "class_agreement",
      "posterior_max_abs_diff_4dp", "exact_classes", "exact_saved_posterior_4dp",
      "fit_md5", "scoring_state_md5"), "MIMIC exact replay audit")
    if (nrow(audit) != 1L || anyNA(audit) || audit$n_patients != 9804L ||
        audit$n_class_changed != 0L || audit$class_agreement != 1 ||
        audit$posterior_max_abs_diff_4dp != 0 || !isTRUE(audit$exact_classes) ||
        !isTRUE(audit$exact_saved_posterior_4dp) ||
        audit$fit_md5 != state$source_fit_md5 ||
        audit$scoring_state_md5 != expected["frozen_scoring_state"])
      stop("MIMIC audit does not verify the original frozen classification.")
    write_csv(data.frame(role = names(inputs), filename = basename(inputs),
      source_path = unname(inputs), md5_before = unname(input_hashes)), "input_audit")
    reproduce <- c(paste0("# ", VERSION, ": complete independently executable reproduction."),
                   "main_external_validation <-", deparse(SELF, width.cutoff = 500L),
                   "main_external_validation()")
    write_text(reproduce, "reproduce_validation", "R")
    engine_md5 <- unname(tools::md5sum(filename("reproduce_validation", "R")))
    d <- read_csv(inputs["external_raw"], c(ID, PATIENT_ID))
    m <- read_csv(inputs["MIMIC_raw_reference"], c("stay_id", "subject_id"))
    required <- unique(c(ID, PATIENT_ID, unlist(RAW_GROUPS, use.names = FALSE), CORE,
      "sepsis_offset_sec", "landmark_offset_sec", "make30_risk_start_offset_sec",
      "make30_censor_offset_sec", "make30_followup_days_after_landmark", "make30_evaluable_window",
      "make30", "make30_death", "make30_new_krt", "make30_persistent_kidney_dysfunction",
      "prior_krt_before_landmark", "first_new_krt_offset_sec",
      "final_inpatient_creatinine", "final_inpatient_creatinine_offset_sec", "baseline_creatinine"))
    assert_columns(d, required, DB)
    if (anyNA(d[[PATIENT_ID]]) || anyDuplicated(d[[PATIENT_ID]]))
      stop("External cohort differs from the inspected one-patient-per-record export.")
    dd <- read_classes(inputs["external_classification"], d, ID, EXPECTED_N, EXPECTED_CLASS_COUNTS)
    d <- dd$data
    mm <- read_classes(inputs["MIMIC_formal_classification"], m, "stay_id", 9804L,
                       c(4291L, 2841L, 2672L), rounded = TRUE)
    m <- mm$data
    for (v in unlist(RAW_GROUPS, use.names = FALSE)) {
      check_numeric(d[[v]], paste0(DB, ":", v))
      check_numeric(m[[v]], paste0("MIMIC:", v))
    }
    for (v in c("make30", "make30_death", "make30_new_krt", "make30_persistent_kidney_dysfunction"))
      check_binary(d[[v]], v)
    check_binary(d$make30_evaluable_window, "make30_evaluable_window", FALSE)
    positive <- d$make30_death == 1 | d$make30_new_krt == 1 | d$make30_persistent_kidney_dysfunction == 1
    negative <- d$make30_death == 0 & d$make30_new_krt == 0 & d$make30_persistent_kidney_dysfunction == 0
    derived <- ifelse(positive %in% TRUE, 1L, ifelse(negative %in% TRUE, 0L, NA_integer_))
    if (!identical(is.na(derived), is.na(d$make30)) ||
        any(derived[!is.na(derived)] != d$make30[!is.na(d$make30)]))
      stop("Stored composite conflicts with its three stored components.")
    eligible <- d$make30_evaluable_window == 1L
    for (v in c("sepsis_offset_sec", "landmark_offset_sec", "make30_risk_start_offset_sec",
                "make30_censor_offset_sec", "make30_followup_days_after_landmark"))
      check_numeric(d[[v]], v, FALSE)
    if (any(d$landmark_offset_sec - d$sepsis_offset_sec != 86400) ||
        any(d$make30_risk_start_offset_sec != d$landmark_offset_sec) ||
        max(abs((d$make30_censor_offset_sec - d$make30_risk_start_offset_sec) / 86400 -
                  d$make30_followup_days_after_landmark)) > 1e-6 ||
        any(eligible != (d$make30_censor_offset_sec > d$make30_risk_start_offset_sec)) ||
        any(!eligible & !is.na(d$make30)))
      stop("Stored landmark/risk/follow-up/eligibility fields are inconsistent.")
    assert_columns(d, "hospital_discharge_offset_sec", "SICdb stored time interface")
    censor_after_day30 <- eligible & d$make30_censor_offset_sec > 30 * 86400
    window_review_required <- any(censor_after_day30)
    positive_krt <- d$make30_new_krt %in% 1L
    positive_renal <- d$make30_persistent_kidney_dysfunction %in% 1L
    if (anyNA(d$first_new_krt_offset_sec[positive_krt]) ||
        any(d$first_new_krt_offset_sec[positive_krt] <= d$make30_risk_start_offset_sec[positive_krt]) ||
        any(d$first_new_krt_offset_sec[positive_krt] > d$make30_censor_offset_sec[positive_krt]))
      stop("A stored positive SICdb new-KRT event is outside its stored observation window.")
    krt_after_day30 <- positive_krt & d$first_new_krt_offset_sec > 30 * 86400
    renal_measurement_after_day30 <- positive_renal &
      d$final_inpatient_creatinine_offset_sec > 30 * 86400
    write_csv(data.frame(dataset = DB, retained_n = nrow(d), evaluable_n = sum(eligible),
      non_evaluable_n = sum(!eligible), post_sepsis_landmark_hours = 24L,
      time_unit = "seconds", checked_origin = "exported offset zero",
      evaluable_censor_after_offset_day30_n = sum(censor_after_day30),
      known_MAKE30_with_censor_after_offset_day30_n = sum(censor_after_day30 & !is.na(d$make30)),
      new_KRT_event_after_offset_day30_n = sum(krt_after_day30, na.rm = TRUE),
      positive_renal_final_measurement_after_offset_day30_n = sum(renal_measurement_after_day30, na.rm = TRUE),
      exact_ICU_admission_day30_window_established = !window_review_required,
      stored_outcomes_retained = TRUE, source_data_modified = FALSE), "time_axis_audit")
    detail <- data.frame(unit_id = as.character(d[[ID]]),
      sepsis_offset_sec = d$sepsis_offset_sec, landmark_offset_sec = d$landmark_offset_sec,
      censor_offset_sec = d$make30_censor_offset_sec,
      hospital_discharge_offset_sec = d$hospital_discharge_offset_sec,
      followup_days_after_landmark = d$make30_followup_days_after_landmark,
      make30_evaluable_window = d$make30_evaluable_window, stored_make30 = d$make30,
      censor_days_from_offset_zero = d$make30_censor_offset_sec / 86400,
      evaluable_censor_after_offset_day30 = censor_after_day30,
      positive_KRT_after_offset_day30 = krt_after_day30,
      positive_renal_measurement_after_offset_day30 = renal_measurement_after_day30)
    write_csv(detail, "time_window_patient_audit")
    write_csv(detail[censor_after_day30, , drop = FALSE], "time_window_review_records")
    if (window_review_required) progress_log("SICdb stored endpoints retained; ",
      sum(censor_after_day30), " evaluable censor offsets exceed offset day30. See time_axis_audit.")
    keep <- eligible & !is.na(d$make30)
    a <- d[keep, , drop = FALSE]
    a$sex_reported <- as.character(a$sex)
    for (v in setdiff(CORE, "sex")) check_numeric(a[[v]], v, FALSE)
    for (v in c("ckd_or_renal_disease", "diabetes", "mechanical_ventilation_0_24h", "vasopressor_0_24h"))
      check_binary(a[[v]], v, FALSE)
    sex_map <- c(female = "Female", male = "Male")
    if (anyNA(a$sex_reported) || !all(a$sex_reported %in% names(sex_map)))
      stop("Unexpected SICdb sex codes; require the inspected literal female/male codes.")
    a$sex <- factor(unname(sex_map[a$sex_reported]), levels = c("Female", "Male"))
    if (anyNA(a[CORE])) stop("Unexpected missing SICdb core covariates.")
    if (!all(CLASSES %in% as.character(a$phenotype)) || length(unique(a$make30)) != 2L)
      stop("Analysis cohort lacks a phenotype or an outcome category.")
    progress_log("Exact classifications and stored outcomes verified; full n=", nrow(d),
                 "; analysis n=", nrow(a), "; events=", sum(a$make30), ".")
    write_csv(dd$classes, "fixed_classification")
    write_csv(d, "classified_cohort")
    write_csv(a, "analysis_cohort")
    write_csv(data.frame(dataset = DB, full_classified_n = nrow(d),
      outside_evaluable_window_n = sum(!eligible),
      unknown_MAKE30_inside_evaluable_window_n = sum(eligible & is.na(d$make30)),
      total_excluded_n = sum(!keep), analyzed_n = nrow(a), MAKE30_events = sum(a$make30),
      MAKE30_event_pct = 100 * mean(a$make30), ordinary_covariate_imputation_n = sum(is.na(a$sex))),
      "analysis_flow")
    write_csv(do.call(rbind, lapply(CLASSES, function(k) {
      ix <- d$phenotype == k
      data.frame(dataset = DB, phenotype = k, full_n = sum(ix), non_evaluable_n = sum(ix & !eligible),
        unknown_MAKE30_in_evaluable_n = sum(ix & eligible & is.na(d$make30)),
        analyzed_n = sum(ix & keep), MAKE30_events = sum(d$make30[ix & keep]))
    })), "flow_by_class")
    write_csv(data.frame(variable = CORE,
      missing_full_cohort = vapply(d[CORE], function(x) sum(is.na(x)), integer(1L)),
      missing_analysis_model_covariate = vapply(a[CORE], function(x) sum(is.na(x)), integer(1L)),
      imputation_used = vapply(a[CORE], anyNA, logical(1L))), "covariate_missingness")
    p <- as.matrix(d[c("p_class1", "p_class2", "p_class3")])
    lp <- matrix(0, nrow(p), ncol(p)); nz <- p > 0; lp[nz] <- p[nz] * base::log(p[nz])
    entropy <- 1 + sum(lp) / (nrow(p) * base::log(3))
    quality <- do.call(rbind, lapply(CLASSES, function(k) {
      ix <- d$phenotype == k
      data.frame(dataset = DB, phenotype = k, n = sum(ix), pct = 100 * mean(ix),
        APPA = mean(d[[paste0("p_class", match(k, CLASSES))]][ix]),
        mean_max_posterior = mean(d$max_posterior[ix]), median_max_posterior = stats::median(d$max_posterior[ix]),
        pct_max_posterior_ge_070 = 100 * mean(d$max_posterior[ix] >= 0.7),
        pct_max_posterior_ge_090 = 100 * mean(d$max_posterior[ix] >= 0.9),
        normalized_entropy = entropy)
    }))
    write_csv(quality, "classification_quality")
    proportions <- rbind(
      data.frame(dataset = "MIMIC", phenotype = CLASSES, n = as.integer(table(m$phenotype)),
                 pct = 100 * as.integer(table(m$phenotype)) / nrow(m)),
      quality[c("dataset", "phenotype", "n", "pct")])
    write_csv(proportions, "MIMIC_class_proportion_comparison")
    observed <- rbind(observed_summary(m, "MIMIC"), observed_summary(d, DB))
    write_csv(observed, "observed_trajectory_summary")
    write_csv(observed[c("dataset", "phenotype", "indicator", "source_column", "time_hour",
                         "class_n", "observed_n", "missing_n", "missing_pct")], "trajectory_missingness")
    write_csv(data.frame(variable = Z, center = unname(state$raw_center),
                         scale = unname(state$raw_scale)), "frozen_standardization")
    trajectory_comparison <- do.call(rbind, lapply(CLASSES, function(k) {
      do.call(rbind, lapply(INDICATORS, function(indicator) {
        ref <- observed[observed$dataset == "MIMIC" & observed$phenotype == k & observed$indicator == indicator, ]
        ext <- observed[observed$dataset == DB & observed$phenotype == k & observed$indicator == indicator, ]
        stopifnot(identical(ref$time_hour, ext$time_hour))
        delta <- ext$mean_native - ref$mean_native
        data.frame(dataset = DB, phenotype = k, indicator = indicator,
          mean_level_difference_native = mean(delta), RMSE_native = sqrt(mean(delta^2)),
          baseline_centered_RMSE_z = sqrt(mean((ext$baseline_centered_z - ref$baseline_centered_z)^2)),
          MIMIC_end_minus_start_native = tail(ref$mean_native, 1L) - ref$mean_native[1L],
          external_end_minus_start_native = tail(ext$mean_native, 1L) - ext$mean_native[1L])
      }))
    }))
    write_csv(trajectory_comparison, "MIMIC_trajectory_comparison")
    rates <- endpoint_rates(d, DB)
    write_csv(rates, "MAKE30_and_components_rates")
    main_rates <- rates[rates$endpoint == "MAKE30", , drop = FALSE]
    write_csv(rbind(endpoint_rates(m, "MIMIC"), rates), "MIMIC_outcome_rates_comparison")
    write_csv(baseline_table(d), "Table1_full_cohort")
    table_a <- d[keep, , drop = FALSE]
    write_csv(baseline_table(table_a), "Table1_analysis_cohort")
    signature <- list(schema = paste0(VERSION, "_", DB, "_external_G3"),
      engine_md5 = engine_md5, input_md5 = input_hashes,
      analysis_md5 = unname(tools::md5sum(filename("analysis_cohort"))),
      adjusted_covariates = CORE, reference_class = "C1", primary = "modified_Poisson_HC1_RR",
      analysis_n = nrow(a), events = sum(a$make30), R_version = as.character(getRversion()),
      MI = list(m = MI_M, maxit = MI_MAXIT, seed = MI_SEED))
        individual_fits <- list()
    fits <- list()
    for (estimator in c("Poisson_RR", "Logistic_OR")) {
      z <- fit_one(a, estimator, FALSE, "crude")
      fits[[length(fits) + 1L]] <- z
      individual_fits[[length(individual_fits) + 1L]] <- z
    }
    if (anyNA(a$sex)) {
      imputed <- make_imputations()
      for (estimator in c("Poisson_RR", "Logistic_OR")) {
        mi_fits <- lapply(seq_len(MI_M), function(j)
          fit_one(imputed[[j]], estimator, TRUE, "adjusted_common_covariates", imputation = j))
        individual_fits <- c(individual_fits, mi_fits)
        fits[[length(fits) + 1L]] <- pool_models(mi_fits, "adjusted_common_covariates_MICE_20")
        cc <- fit_one(a[!is.na(a$sex), , drop = FALSE], estimator, TRUE,
                      "adjusted_binary_sex_complete_case_sensitivity")
        fits[[length(fits) + 1L]] <- cc
        individual_fits[[length(individual_fits) + 1L]] <- cc
      }
    } else {
      for (estimator in c("Poisson_RR", "Logistic_OR")) {
        z <- fit_one(a, estimator, TRUE, "adjusted_common_covariates")
        fits[[length(fits) + 1L]] <- z
        individual_fits[[length(individual_fits) + 1L]] <- z
      }
    }
    contrasts <- do.call(rbind, lapply(fits, contrast_rows))
    rr <- contrasts[contrasts$estimator == "Poisson_RR", , drop = FALSE]
    names(rr)[names(rr) == "estimate"] <- "RR"
    names(rr)[names(rr) == "SE_log_effect"] <- "robust_SE_logRR"
    ors <- contrasts[contrasts$estimator == "Logistic_OR", , drop = FALSE]
    names(ors)[names(ors) == "estimate"] <- "OR"
    names(ors)[names(ors) == "SE_log_effect"] <- "SE_logOR"
    write_csv(rr, "MAKE30_modified_Poisson_RR")
    write_csv(ors, "MAKE30_supplementary_logistic_OR")
    coefficients <- do.call(rbind, lapply(fits, function(z) {
      b <- z$coefficients; se <- sqrt(diag(z$covariance)); critical <- stats::qt(0.975, df = z$df)
      data.frame(estimator = z$estimator, model = z$model_name, term = names(b),
        log_coefficient = unname(b), SE = unname(se), effect = unname(exp(b)),
        CI95_low = unname(exp(b - critical * se)), CI95_high = unname(exp(b + critical * se)),
        p_value = unname(2 * stats::pt(-abs(b / se), df = z$df)),
        df = unname(z$df), imputations = z$imputations, row.names = NULL)
    }))
    write_csv(coefficients, "regression_coefficients")
    global_tests <- do.call(rbind, lapply(fits, function(z) {
      terms <- c("phenotypeC2", "phenotypeC3")
      b <- z$coefficients[terms]; V <- z$covariance[terms, terms, drop = FALSE]
      chi <- as.numeric(crossprod(b, solve(V, b)))
      data.frame(estimator = z$estimator, model = z$model_name,
        test = if (z$imputations) "Asymptotic joint Wald test using Rubin-pooled covariance" else "Joint phenotype Wald test",
        chi_square = chi, df = 2L, p_value = stats::pchisq(chi, 2L, lower.tail = FALSE))
    }))
    write_csv(global_tests, "global_phenotype_tests")
    diagnostics <- do.call(rbind, lapply(individual_fits, function(z)
      data.frame(estimator = z$estimator, model = z$model_name, imputation = z$imputation,
        n = z$n, model_rank = z$npar, converged = z$converged,
        covariance_method = if (z$estimator == "Poisson_RR") "HC1 sandwich: n/(n-p)" else "Model-based binomial",
        fitted_min = z$fitted_min, fitted_max = z$fitted_max,
        n_fitted_above_one = z$n_fitted_above_one, warnings = paste(z$warnings, collapse = " | "))))
    write_csv(diagnostics, "model_diagnostics")
    adjusted_model <- if (anyNA(a$sex)) "adjusted_common_covariates_MICE_20" else "adjusted_common_covariates"
    adj <- rr[rr$model == adjusted_model, , drop = FALSE]
    ref_rr <- read_csv(inputs["MIMIC_RR_reference"])
    assert_columns(ref_rr, c("dataset", "model", "comparison", "RR", "CI95_low", "CI95_high",
                             "robust_SE_logRR", "p_value", "n", "events"), "MIMIC RR reference")
    ref_rr <- ref_rr[ref_rr$model == "adjusted_common_covariates", , drop = FALSE]
    if (nrow(ref_rr) != 2L || !identical(as.character(ref_rr$comparison), c("C2 vs C1", "C3 vs C1")) ||
        !all(ref_rr$dataset == "MIMIC") || !all(ref_rr$n == 9758L) || !all(ref_rr$events == 2442L) ||
        anyNA(ref_rr[c("RR", "CI95_low", "CI95_high", "p_value")]))
      stop("MIMIC RR reference does not match the reviewed primary analysis.")
    if (max(abs(ref_rr$RR - c(3.355261, 2.357394))) > 1e-5)
      stop("MIMIC reference RR differs materially from the uploaded formal result.")
    comparison_fields <- c("dataset", "model", "comparison", "RR", "CI95_low",
                           "CI95_high", "robust_SE_logRR", "p_value", "n", "events")
    rr_comparison <- rbind(ref_rr[comparison_fields], adj[comparison_fields])
    write_csv(rr_comparison, "MIMIC_adjusted_RR_comparison")
        plot_save("native_trajectories", draw_native, width = 11, height = 7)
    plot_save("joint_standardized_trajectories", function() draw_joint(FALSE), width = 12, height = 5.5)
    plot_save("joint_baseline_centered_trajectories", function() draw_joint(TRUE), width = 12, height = 5.5)
    plot_save("MIMIC_native_trajectory_comparison", function() draw_comparison(FALSE), width = 11, height = 13)
    plot_save("MIMIC_shape_trajectory_comparison", function() draw_comparison(TRUE), width = 11, height = 13)
    plot_save("MAKE30_rates", function() draw_rates(FALSE), width = 7.5, height = 5.2)
    plot_save("MAKE30_component_rates", function() draw_rates(TRUE), width = 10, height = 7.5)
    plot_save("MIMIC_adjusted_RR_comparison", draw_rr, width = 10.5, height = 6.0)
    after <- unname(tools::md5sum(inputs))
    if (!identical(after, unname(input_hashes))) stop("An original input changed during analysis.")
    write_csv(data.frame(role = names(inputs), filename = basename(inputs), source_path = unname(inputs),
      md5_before = unname(input_hashes), md5_after = after, unchanged = TRUE), "input_audit")
    completion_status <- if (window_review_required) "COMPLETE_WITH_TIME_WINDOW_REVIEW" else "COMPLETE"
    summary_lines <- c(
      paste0("Status: ", completion_status), paste0("Dataset: ", DB),
      "Phenotype: frozen V295 G3 S01; exact verified V318 external classifications.",
      "No EM refitting, label permutation, external restandardization or G selection.",
      "All original CSV/RDS inputs and stored outcomes are preserved.",
      "Classification/trajectory descriptions use the entire classified external cohort.",
      "MAKE30 association uses only an evaluable post-sepsis-24h window and a known stored MAKE30 outcome.",
      "Composite: death OR new KRT OR persistent kidney dysfunction; unknown outcomes are never non-events.",
      sprintf("Classified n=%d; analysis n=%d; MAKE30 events=%d.", nrow(d), nrow(a), sum(a$make30)),
      sprintf("Excluded: %d non-evaluable; %d evaluable with unknown MAKE30.", sum(!eligible), sum(eligible & is.na(d$make30))),
      paste0("Adjustment: ", paste(CORE, collapse = ", "), "."),
      "Primary association: modified Poisson log-link; HC1 sandwich n/(n-p); C1 reference.",
      "Supplementary association: logistic OR; distinct from RR.",
      "Existing MIMIC results are read for comparison; MIMIC regression is not refitted.",
      "This validates phenotype transport and its outcome associations, not an individual risk-prediction model.",
      "Component rates have their own known-outcome denominators; overlapping events cannot be added.",
      "Observed graphs use available raw measurements at real 1/3/.../23h or 3/9/15/21h anchors.",
      "Their lines connect group means descriptively; no raw trajectory values are imputed.",
      "Baseline-centered z graphs subtract each dataset/class/indicator's first observed mean z value.",
      "A normalized entropy/APPA value describes classification certainty, not external clinical truth.",
      if (DB == "eICU") paste0("Two reported Unknown/Other sex values: MICE m=", MI_M,
        ", maxit=", MI_MAXIT, ", seed=", MI_SEED, "; Rubin pooling; complete-case sensitivity also included.") else
        "All eight analysis covariates are complete; no ordinary covariate imputation.",
      if (window_review_required) paste0("TIME WINDOW REVIEW: ", sum(censor_after_day30),
        " evaluable stored censor offsets exceed 30*86400 from offset zero. Stored endpoints retained.",
        " This script does not establish that SICdb stored outcomes exactly equal ICU-admission-to-day30 outcomes.") else
        "eICU unified ICU-origin landmark/censor/event conversions and day30 ceiling passed.",
      "Regression and imputation checkpoints are reused only when code, input hashes and cohort match.",
      "Compact checkpoints retain coefficient/covariance/formula/diagnostics without repeated large patient matrices.",
      "", "Observed MAKE30:",
      vapply(seq_len(nrow(main_rates)), function(i) sprintf("%s: %d/%d, %.2f%% (95%% CI %.2f-%.2f%%).",
        main_rates$phenotype[i], main_rates$events[i], main_rates$known_outcome_n[i],
        main_rates$event_pct[i], main_rates$CI95_low_pct[i], main_rates$CI95_high_pct[i]), character(1L)),
      "", "Primary adjusted RR:",
      vapply(seq_len(nrow(adj)), function(i) sprintf("%s: RR %.4f (95%% CI %.4f-%.4f), p=%.6g.",
        adj$comparison[i], adj$RR[i], adj$CI95_low[i], adj$CI95_high[i], adj$p_value[i]), character(1L)))
    write_text(summary_lines, "validation_summary")
    write_text(capture.output(utils::sessionInfo()), "session_info")
    write_csv(data.frame(status = completion_status, dataset = DB, classified_n = nrow(d),
      analysis_n = nrow(a), events = sum(a$make30), raw_md5 = input_hashes["external_raw"],
      class_md5 = input_hashes["external_classification"], engine_md5 = engine_md5,
      completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "completion")
    progress_log("Analysis completed; ", completion_status, "; all original hashes unchanged.")
    print(adj[c("comparison", "RR", "CI95_low", "CI95_high", "p_value")], row.names = FALSE)
    list(ok = TRUE, status = completion_status, adjusted_RR = adj)
  }, error = function(e) {
    write_text(c("Status: ERROR", conditionMessage(e)), "validation_summary")
    write_text(capture.output(utils::sessionInfo()), "session_info")
    write_csv(data.frame(status = "ERROR", dataset = DB, message = conditionMessage(e)), "completion")
    progress_log("ERROR: ", conditionMessage(e))
    list(ok = FALSE, status = "ERROR", error = conditionMessage(e))
  })
  files <- list.files(out_dir, full.names = TRUE)
  files <- files[!grepl("\\.tmp_[0-9]+$", basename(files))]
  write_csv(data.frame(filename = basename(files), size_bytes = file.info(files)$size,
                       md5 = unname(tools::md5sum(files))), "output_manifest")
  zip_path <- archive_results()
  cat("\nOutput directory: ", out_dir, "\nResult ZIP: ", zip_path, "\n", sep = "")
  if (!isTRUE(status$ok)) message("Analysis stopped. Send this ZIP; completed compatible checkpoints are retained.")
  invisible(list(status = status, output_directory = out_dir, zip = zip_path))
}

main_external_validation()
