# V325: MSI-1 independent paths; V320 statistical implementation retained.
# V320: MIMIC MAKE30 / RStudio direct-run edition.
# 使用方法：在 RStudio 中打开本文件，Ctrl+A 全选，然后 Ctrl+Enter 运行。
# 自动读取 Z:/Research/MSI-1 内正式归档输入。
# 输入位置已经与上传的 GBMTM_stage5_output.zip 核对；脚本位置不影响读取。
# 输出文件名、输出目录及 checkpoint 保留 V319 格式，支持继续使用 V319 断点。
# 本文件包含完整分析代码，仅依赖 base R。

main_make30 <- function() {
  previous_options <- options(stringsAsFactors = FALSE, scipen = 999, error = NULL,
                              contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(previous_options), add = TRUE)
  VERSION <- "V319"
  RUN_LOGISTIC_OR <- TRUE  # Supplementary; modified Poisson RR remains primary.
  EXPECTED_MIMIC_MD5 <- "f1eb2d7fac1303a682b53c679dde9abd"
  EXPECTED_CLASS_MD5 <- "d6824b387dcc2f3cfd5c72cb82d69e94"
  EXPECTED_FIT_MD5 <- "cf5881479ee206195aa80d59fa3bd17f"
  EXPECTED_STATE_MD5 <- "89f684d00680fafc30577d677fee6a30"
  CORE <- c("age", "sex", "sepsis_sofa_score", "baseline_egfr_mdrd",
            "ckd_or_renal_disease", "diabetes",
            "mechanical_ventilation_0_24h", "vasopressor_0_24h")
  CLASSES <- c("C1", "C2", "C3")

  archive_root <- "Z:/Research/MSI-1"
  output_root <- file.path(archive_root, "STEP5_冻结分型与预后外部验证", "独立复现")
  dir.create(output_root, recursive=TRUE, showWarnings=FALSE)
  rds_output <- file.path(archive_root, "RDS", "独立复现", "MIMIC")
  dir.create(rds_output, recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(output_root)) stop("STEP5 root does not exist: ", output_root)
  # Exact directories verified in the uploaded current STEP5 archive.
  kit_dir <- archive_root
  scoring_dir <- file.path(archive_root,"STEP5_冻结分型与预后外部验证","冻结分类")
  inputs <- c(
    raw_cohort = file.path(kit_dir, "MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv"),
    formal_classification = file.path(scoring_dir, "V318_MIMIC_formal_classification_G3.csv"),
    replay_audit = file.path(scoring_dir, "V318_MIMIC_replay_audit_G3.csv")
  )
  if (any(!file.exists(inputs)))
    stop("Missing exact input files:\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
  inputs <- vapply(inputs, normalizePath, character(1L), winslash = "/", mustWork = TRUE)
  input_hashes <- unname(tools::md5sum(inputs))
  names(input_hashes) <- names(inputs)
  if (input_hashes["raw_cohort"] != EXPECTED_MIMIC_MD5)
    stop("MIMIC input MD5 differs from the verified V307 cohort.")
  if (input_hashes["formal_classification"] != EXPECTED_CLASS_MD5)
    stop("Formal V318 classification MD5 differs from the reviewed G3 file.")

  out_dir <- file.path(output_root, "V319_STEP5_MIMIC_MAKE30_G3")
  dir.create(out_dir, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Cannot create output directory: ", out_dir)
  filename <- function(content, extension = "csv")
    file.path(if (identical(extension,"rds")) rds_output else out_dir, paste0(VERSION, "_", content, "_G3.", extension))
  log_path <- filename("MAKE30_log", "txt")
  progress_log <- function(...) {
    line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                   paste0(..., collapse = ""))
    cat(line, "\n", sep = "")
    cat(line, "\n", file = log_path, append = TRUE, sep = "")
  }
  replace_generated <- function(tmp, dest) {
    if (file.exists(dest) && unlink(dest) != 0L) stop("Cannot replace generated output: ", dest)
    if (!file.rename(tmp, dest)) stop("Cannot save generated output: ", dest)
    invisible(dest)
  }
  write_csv <- function(d, content) {
    dest <- filename(content)
    tmp <- paste0(dest, ".tmp_", Sys.getpid())
    csv_options <- options(scipen = 0)
    on.exit(options(csv_options), add = TRUE)
    utils::write.csv(d, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
    replace_generated(tmp, dest)
  }
  write_text <- function(lines, content) {
    dest <- filename(content, "txt")
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
  read_exact_csv <- function(path, id = NULL) {
    if (is.null(id)) utils::read.csv(path, check.names = FALSE,
      na.strings = c("", "NA", "NULL")) else
      utils::read.csv(path, check.names = FALSE, colClasses = setNames("character", id),
                      na.strings = c("", "NA", "NULL"))
  }
  assert_columns <- function(d, required, label) {
    if (anyDuplicated(names(d)) || !all(required %in% names(d)))
      stop(label, " schema differs. Required missing fields: ",
           paste(setdiff(required, names(d)), collapse = ", "))
  }
  check_binary <- function(x, label, allow_missing = TRUE) {
    if (!is.numeric(x) || any(!is.na(x) & !x %in% c(0, 1)) ||
        (!allow_missing && anyNA(x))) stop("Invalid binary coding: ", label)
  }
  wilson <- function(events, n) {
    if (n == 0L) return(c(NA_real_, NA_real_))
    z <- stats::qnorm(0.975); p <- events / n
    center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    c(center - half, center + half)
  }
  progress_log("MIMIC MAKE30 started; formal G3 labels fixed; primary modified Poisson RR.")

  status <- tryCatch({
    d <- read_exact_csv(inputs["raw_cohort"], "stay_id")
    cl <- read_exact_csv(inputs["formal_classification"], "unit_id")
    audit <- read_exact_csv(inputs["replay_audit"])
    assert_columns(audit, c("n_patients", "n_class_changed", "class_agreement",
      "posterior_max_abs_diff_4dp", "exact_classes", "exact_saved_posterior_4dp",
      "fit_md5", "scoring_state_md5"), "Replay audit")
    if (nrow(audit) != 1L || anyNA(audit) || audit$n_patients != 9804L ||
        audit$n_class_changed != 0L || audit$class_agreement != 1 ||
        audit$posterior_max_abs_diff_4dp != 0 || !isTRUE(audit$exact_classes) ||
        !isTRUE(audit$exact_saved_posterior_4dp) || audit$fit_md5 != EXPECTED_FIT_MD5 ||
        audit$scoring_state_md5 != EXPECTED_STATE_MD5)
      stop("Replay audit does not establish exact original G3 classification.")
    assert_columns(cl, c("unit_id", "assigned_class", "p_class1", "p_class2", "p_class3"),
                   "Formal G3 classification")
    required <- unique(c("subject_id", "hadm_id", "stay_id", "sepsis_time", "landmark_time",
      "make30_risk_start", "make30_censor_time", "make30_followup_days_after_landmark",
      "make30_evaluable_window", "make30", "make30_death", "make30_new_krt",
      "make30_persistent_kidney_dysfunction", CORE, "race_ethnicity", "weight_kg",
      "bmi", "baseline_creatinine", "baseline_creatinine_source",
      "charlson_comorbidity_index", "congestive_heart_failure", "liver_disease", "malignancy"))
    assert_columns(d, required, "MIMIC")
    if (nrow(d) != 9804L || nrow(cl) != 9804L || anyNA(d$stay_id) ||
        anyNA(cl$unit_id) || anyDuplicated(d$stay_id) || anyDuplicated(cl$unit_id) ||
        anyDuplicated(d$subject_id) || !setequal(d$stay_id, cl$unit_id) ||
        anyNA(cl$assigned_class) || !all(cl$assigned_class %in% 1:3))
      stop("MIMIC patient identities or formal classification differ.")
    pp <- as.matrix(cl[c("p_class1", "p_class2", "p_class3")])
    if (anyNA(pp) || any(!is.finite(pp) | pp < 0 | pp > 1) ||
        max(abs(rowSums(pp) - 1)) > 0.000151)
      stop("Stored four-decimal posterior is invalid.")
    index <- match(d$stay_id, cl$unit_id)
    d$assigned_class <- cl$assigned_class[index]
    d$phenotype <- factor(paste0("C", d$assigned_class), levels = CLASSES)
    for (v in c("p_class1", "p_class2", "p_class3")) d[[v]] <- cl[[v]][index]
    if (!identical(as.integer(table(d$phenotype)), c(4291L, 2841L, 2672L)))
      stop("Formal G3 class counts differ from the reviewed result.")
    for (v in c("make30", "make30_death", "make30_new_krt",
                "make30_persistent_kidney_dysfunction")) check_binary(d[[v]], v)
    check_binary(d$make30_evaluable_window, "make30_evaluable_window", FALSE)
    any_positive <- d$make30_death == 1 | d$make30_new_krt == 1 |
      d$make30_persistent_kidney_dysfunction == 1
    all_negative <- d$make30_death == 0 & d$make30_new_krt == 0 &
      d$make30_persistent_kidney_dysfunction == 0
    derived <- ifelse(any_positive %in% TRUE, 1L,
                      ifelse(all_negative %in% TRUE, 0L, NA_integer_))
    if (!identical(is.na(derived), is.na(d$make30)) ||
        any(derived[!is.na(derived)] != d$make30[!is.na(d$make30)]))
      stop("Stored MAKE30 conflicts with its three component indicators.")
    # Validation only: stored outcome values and dates are never replaced.
    parse_time <- function(x) as.POSIXct(x, format = "%d/%m/%Y %H:%M:%S", tz = "UTC")
    beg <- parse_time(d$sepsis_time); lm <- parse_time(d$landmark_time)
    censor <- parse_time(d$make30_censor_time)
    if (anyNA(beg) || anyNA(lm) || anyNA(censor) ||
        max(abs(as.numeric(difftime(lm, beg, units = "hours")) - 24)) > 1e-5 ||
        any(d$make30_risk_start != d$landmark_time) ||
        any(censor < lm) || anyNA(d$make30_followup_days_after_landmark) ||
        max(abs(as.numeric(difftime(censor, lm, units = "days")) -
                d$make30_followup_days_after_landmark)) > 1e-6)
      stop("Stored landmark/risk/censor dates are inconsistent.")
    eligible <- d$make30_evaluable_window == 1L
    if (any(eligible != (censor > lm)) || any(!eligible & !is.na(d$make30)))
      stop("Stored outcome eligibility is inconsistent.")
    keep <- eligible & !is.na(d$make30)
    if (!all(d$sex %in% c("F", "M"))) stop("Unexpected MIMIC sex codes; require exact F/M.")
    d$sex <- factor(d$sex, levels = c("F", "M"), labels = c("Female", "Male"))
    for (v in setdiff(CORE, "sex"))
      if (!is.numeric(d[[v]]) || any(!is.finite(d[[v]][!is.na(d[[v]])])))
        stop("Invalid numeric core covariate: ", v)
    for (v in c("ckd_or_renal_disease", "diabetes",
                "mechanical_ventilation_0_24h", "vasopressor_0_24h"))
      check_binary(d[[v]], v)
    a <- d[keep, , drop = FALSE]
    if (anyNA(a[CORE]))
      stop("Core covariates are incomplete; do not silently delete or impute records.")
    if (!all(CLASSES %in% as.character(a$phenotype)) || length(unique(a$make30)) != 2L)
      stop("The eligible analysis cohort lacks a phenotype or an outcome category.")
    progress_log("All identities/outcomes verified. Full n=", nrow(d),
      "; MAKE30 analysis n=", nrow(a), "; events=", sum(a$make30), ".")
    write_csv(data.frame(role = names(inputs), filename = basename(inputs),
      source_path = unname(inputs), md5_before = unname(input_hashes)), "MIMIC_input_audit")
    flow <- data.frame(dataset = "MIMIC", full_classified_n = nrow(d),
      outside_evaluable_window_n = sum(!eligible),
      unknown_MAKE30_inside_evaluable_window_n = sum(eligible & is.na(d$make30)),
      total_excluded_n = sum(!keep), analyzed_n = nrow(a),
      MAKE30_events = sum(a$make30), MAKE30_event_pct = 100 * mean(a$make30))
    write_csv(flow, "MIMIC_analysis_flow")
    flow_by_class <- do.call(rbind, lapply(CLASSES, function(k) {
      ix <- d$phenotype == k
      data.frame(dataset = "MIMIC", phenotype = k, full_n = sum(ix),
        non_evaluable_n = sum(ix & !eligible),
        unknown_MAKE30_in_evaluable_n = sum(ix & eligible & is.na(d$make30)),
        analyzed_n = sum(ix & keep), MAKE30_events = sum(d$make30[ix & keep]))
    }))
    write_csv(flow_by_class, "MIMIC_flow_by_class")
    write_csv(data.frame(variable = CORE,
      missing_full_cohort = vapply(d[CORE], function(x) sum(is.na(x)), integer(1L)),
      missing_analysis_cohort = vapply(a[CORE], function(x) sum(is.na(x)), integer(1L)),
      imputation_used = FALSE), "MIMIC_covariate_missingness")
    write_csv(a[unique(c("subject_id", "hadm_id", "stay_id", "assigned_class", "phenotype",
      "p_class1", "p_class2", "p_class3", "make30", "make30_death", "make30_new_krt",
      "make30_persistent_kidney_dysfunction", CORE, "landmark_time", "make30_censor_time"))],
      "MIMIC_analysis_cohort")

    endpoints <- c(MAKE30 = "make30", Death = "make30_death",
      New_KRT = "make30_new_krt", Persistent_kidney_dysfunction = "make30_persistent_kidney_dysfunction")
    rates <- do.call(rbind, lapply(names(endpoints), function(endpoint) {
      v <- endpoints[[endpoint]]
      do.call(rbind, lapply(CLASSES, function(k) {
        pool <- d[d$phenotype == k & eligible, , drop = FALSE]
        x <- pool[[v]][!is.na(pool[[v]])]
        n <- length(x); events <- sum(x); ci <- wilson(events, n)
        data.frame(dataset = "MIMIC", endpoint = endpoint, phenotype = k,
          eligible_n = nrow(pool), known_outcome_n = n,
          missing_outcome_n = sum(is.na(pool[[v]])), events = events,
          event_pct = 100 * events / n, CI95_low_pct = 100 * ci[1L],
          CI95_high_pct = 100 * ci[2L])
      }))
    }))
    write_csv(rates, "MIMIC_MAKE30_and_components_rates")
    main_rates <- rates[rates$endpoint == "MAKE30", , drop = FALSE]

    continuous <- c("age", "weight_kg", "bmi", "sepsis_sofa_score", "baseline_creatinine",
                    "baseline_egfr_mdrd", "charlson_comorbidity_index")
    binary <- c("ckd_or_renal_disease", "diabetes", "congestive_heart_failure",
                "liver_disease", "malignancy", "mechanical_ventilation_0_24h", "vasopressor_0_24h")
    categorical <- c("sex", "race_ethnicity", "baseline_creatinine_source")
    baseline_table <- function(pool) {
      format_cell <- function(x, kind, level = NULL) {
        known <- x[!is.na(x)]; n <- length(known)
        if (!n) return("NA")
        if (kind == "continuous") {
          q <- stats::quantile(known, c(0.25, 0.5, 0.75), names = FALSE)
          return(sprintf("%.2f [%.2f, %.2f]", q[2L], q[1L], q[3L]))
        }
        num <- if (kind == "binary") sum(known == 1) else sum(as.character(known) == level)
        sprintf("%d/%d (%.1f%%)", num, n, 100 * num / n)
      }
      rows <- list(); j <- 0L
      add_row <- function(v, kind, level = "") {
        row <- data.frame(variable = v, statistic = if (kind == "continuous") "median [IQR]" else "n/known (%)",
                          level = level, missing_n = sum(is.na(pool[[v]])))
        row$Overall <- format_cell(pool[[v]], kind, level)
        for (k in CLASSES) row[[k]] <- format_cell(pool[[v]][pool$phenotype == k], kind, level)
        row
      }
      for (v in continuous) { j <- j + 1L; rows[[j]] <- add_row(v, "continuous") }
      for (v in binary) { check_binary(pool[[v]], v); j <- j + 1L; rows[[j]] <- add_row(v, "binary", "1") }
      for (v in categorical) for (level in sort(unique(as.character(pool[[v]][!is.na(pool[[v]])])))) {
        j <- j + 1L; rows[[j]] <- add_row(v, "categorical", level)
      }
      do.call(rbind, rows)
    }
    write_csv(baseline_table(d), "MIMIC_Table1_full_cohort")
    write_csv(baseline_table(a), "MIMIC_Table1_analysis_cohort")

    signature <- list(schema = "V319_MIMIC_MAKE30_G3", input_md5 = input_hashes,
      analysis_md5 = unname(tools::md5sum(filename("MIMIC_analysis_cohort"))),
      adjusted_covariates = CORE, reference_class = "C1", primary = "modified_Poisson_HC1_RR",
      analysis_n = nrow(a), events = sum(a$make30))
    fit_one <- function(estimator, adjusted) {
      model_name <- if (adjusted) "adjusted_common_covariates" else "crude"
      rhs <- c("phenotype", if (adjusted) CORE else character(0L))
      formula <- stats::as.formula(paste("make30 ~", paste(rhs, collapse = " + ")),
                                  env = baseenv())
      content <- paste0("MIMIC_", estimator, "_", model_name, "_checkpoint")
      dest <- filename(content, "rds")
      if (file.exists(dest)) {
        saved <- readRDS(dest)
        if (!identical(saved$signature, signature) || !identical(saved$estimator, estimator) ||
            !identical(saved$model_name, model_name) ||
            !identical(paste(deparse(saved$fit$formula), collapse = ""),
                       paste(deparse(formula), collapse = "")) ||
            !isTRUE(saved$fit$converged) || stats::nobs(saved$fit) != nrow(a))
          stop("Existing regression checkpoint is incompatible: ", dest)
        progress_log("Resumed completed ", estimator, " / ", model_name, ".")
        return(saved)
      }
      warnings <- character(0L)
      fit <- withCallingHandlers(stats::glm(formula, data = a,
        family = if (estimator == "Poisson_RR") stats::poisson(link = "log") else stats::binomial(link = "logit"),
        na.action = stats::na.fail, x = TRUE, y = TRUE), warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
        })
      if (!isTRUE(fit$converged) || anyNA(stats::coef(fit)) ||
          any(!is.finite(stats::coef(fit))) || stats::nobs(fit) != nrow(a) ||
          fit$rank != ncol(fit$x)) stop("Regression convergence/rank failed: ", estimator, " / ", model_name)
      covariance <- stats::vcov(fit)
      if (estimator == "Poisson_RR") {
        # Exactly retain V307's sandwich calculation and n/(n-p) correction.
        scores <- fit$x * as.numeric(fit$y - stats::fitted(fit))
        covariance <- covariance %*% crossprod(scores) %*% covariance *
          nrow(fit$x) / (nrow(fit$x) - ncol(fit$x))
      }
      if (any(!is.finite(covariance)) || any(diag(covariance) <= 0))
        stop("Regression covariance is invalid: ", estimator, " / ", model_name)
      result <- list(signature = signature, estimator = estimator, model_name = model_name,
                     fit = fit, covariance = covariance, warnings = warnings)
      save_checkpoint(result, content)
      progress_log("Finished ", estimator, " / ", model_name, ".")
      result
    }
    fits <- list(fit_one("Poisson_RR", FALSE), fit_one("Poisson_RR", TRUE))
    if (RUN_LOGISTIC_OR) fits <- c(fits, list(fit_one("Logistic_OR", FALSE), fit_one("Logistic_OR", TRUE)))
    contrast_rows <- function(z) {
      terms <- c("phenotypeC2", "phenotypeC3")
      if (!all(terms %in% names(stats::coef(z$fit)))) stop("Expected C2/C3 contrasts unavailable.")
      b <- stats::coef(z$fit)[terms]; se <- sqrt(diag(z$covariance))[terms]
      data.frame(dataset = "MIMIC", estimator = z$estimator, model = z$model_name,
        comparison = c("C2 vs C1", "C3 vs C1"), estimate = exp(b),
        CI95_low = exp(b - stats::qnorm(0.975) * se),
        CI95_high = exp(b + stats::qnorm(0.975) * se), SE_log_effect = se,
        p_value = 2 * stats::pnorm(-abs(b / se)), n = nrow(a), events = sum(a$make30), row.names = NULL)
    }
    contrasts <- do.call(rbind, lapply(fits, contrast_rows))
    rr <- contrasts[contrasts$estimator == "Poisson_RR", , drop = FALSE]
    names(rr)[names(rr) == "estimate"] <- "RR"
    names(rr)[names(rr) == "SE_log_effect"] <- "robust_SE_logRR"
    write_csv(rr, "MAKE30_modified_Poisson_RR")
    if (RUN_LOGISTIC_OR) {
      ors <- contrasts[contrasts$estimator == "Logistic_OR", , drop = FALSE]
      names(ors)[names(ors) == "estimate"] <- "OR"
      names(ors)[names(ors) == "SE_log_effect"] <- "SE_logOR"
      write_csv(ors, "MAKE30_supplementary_logistic_OR")
    }
    coefficients <- do.call(rbind, lapply(fits, function(z) {
      b <- stats::coef(z$fit); se <- sqrt(diag(z$covariance))
      data.frame(estimator = z$estimator, model = z$model_name, term = names(b),
        log_coefficient = as.numeric(b), SE = as.numeric(se), effect = exp(b),
        CI95_low = exp(b - stats::qnorm(0.975) * se), CI95_high = exp(b + stats::qnorm(0.975) * se),
        p_value = 2 * stats::pnorm(-abs(b / se)), row.names = NULL)
    }))
    write_csv(coefficients, "MAKE30_regression_coefficients")
    global_tests <- do.call(rbind, lapply(fits, function(z) {
      terms <- c("phenotypeC2", "phenotypeC3")
      b <- stats::coef(z$fit)[terms]; V <- z$covariance[terms, terms, drop = FALSE]
      chi <- as.numeric(crossprod(b, solve(V, b)))
      data.frame(estimator = z$estimator, model = z$model_name,
                 test = "Joint phenotype Wald test", chi_square = chi, df = 2L,
                 p_value = stats::pchisq(chi, 2L, lower.tail = FALSE))
    }))
    write_csv(global_tests, "MAKE30_global_phenotype_tests")
    diagnostics <- do.call(rbind, lapply(fits, function(z)
      data.frame(estimator = z$estimator, model = z$model_name,
        n = stats::nobs(z$fit), model_rank = z$fit$rank, converged = z$fit$converged,
        covariance_method = if (z$estimator == "Poisson_RR") "HC1 sandwich: n/(n-p)" else "Model-based binomial",
        fitted_min = min(stats::fitted(z$fit)), fitted_max = max(stats::fitted(z$fit)),
        n_fitted_above_one = sum(stats::fitted(z$fit) > 1),
        warnings = paste(z$warnings, collapse = " | "))))
    write_csv(diagnostics, "MAKE30_model_diagnostics")

    plot_save <- function(content, draw, width = 7.2, height = 5.0) {
      grDevices::pdf(filename(content, "pdf"), width = width, height = height, useDingbats = FALSE)
      tryCatch(draw(), finally = grDevices::dev.off())
      if (isTRUE(capabilities("png"))) {
        png_args <- list(filename = filename(content, "png"), width = width * 300,
                         height = height * 300, res = 300)
        if (isTRUE(capabilities("cairo"))) png_args$type <- "cairo"
        do.call(grDevices::png, png_args)
        tryCatch(draw(), finally = grDevices::dev.off())
      }
    }
    palette <- c("#287C8E", "#C5543E", "#6B68A0")
    plot_save("MIMIC_MAKE30_rates", function() {
      graphics::par(mar = c(5.1, 4.8, 3.3, 1.1))
      locations <- graphics::barplot(main_rates$event_pct, names.arg = main_rates$phenotype,
        col = palette, border = NA, ylim = c(0, max(main_rates$CI95_high_pct) * 1.20),
        ylab = "MAKE30, %", main = "MIMIC: observed MAKE30 by frozen G3 phenotype")
      graphics::arrows(locations, main_rates$CI95_low_pct, locations,
                       main_rates$CI95_high_pct, angle = 90, code = 3, length = 0.05)
      graphics::text(locations, main_rates$CI95_high_pct + 2.2,
                    sprintf("%.1f%%", main_rates$event_pct), cex = 1.0)
      graphics::mtext(sprintf("C1: %d/%d    C2: %d/%d    C3: %d/%d; Wilson 95%% CI",
        main_rates$events[1L], main_rates$known_outcome_n[1L],
        main_rates$events[2L], main_rates$known_outcome_n[2L],
        main_rates$events[3L], main_rates$known_outcome_n[3L]), side = 1, line = 3.0, cex = 0.8)
    })
    adj <- rr[rr$model == "adjusted_common_covariates", , drop = FALSE]
    plot_save("MIMIC_MAKE30_adjusted_RR", function() {
      graphics::par(mar = c(6.4, 6.2, 3.3, 1.0))
      xmax <- max(5, adj$CI95_high * 1.4)
      graphics::plot(NA_real_, NA_real_, xlim = c(0.75, xmax), ylim = c(0.5, 2.8), log = "x",
        xlab = "Adjusted relative risk (95% CI)", ylab = "", yaxt = "n", xaxt = "n",
        main = "MIMIC: MAKE30 association, C1 reference", bty = "n")
      graphics::axis(1, at = c(1, 2, 3, 4, 5))
      yy <- c(2, 1)
      graphics::axis(2, at = yy, labels = adj$comparison, las = 1, tick = FALSE)
      graphics::abline(v = 1, lty = 2, col = "#888888")
      graphics::segments(adj$CI95_low, yy, adj$CI95_high, yy, lwd = 2.0, col = palette[2:3])
      graphics::points(adj$RR, yy, pch = 15, cex = 1.3, col = palette[2:3])
      graphics::text(adj$RR, yy + 0.28,
        sprintf("%.2f (%.2f-%.2f)", adj$RR, adj$CI95_low, adj$CI95_high), cex = 0.95)
      graphics::mtext("Age, sex, SOFA, baseline eGFR, CKD, diabetes, ventilation and vasopressors",
                      side = 1, line = 4.6, cex = 0.69)
    })
    after <- unname(tools::md5sum(inputs))
    if (!identical(after, unname(input_hashes))) stop("An original input changed during analysis.")
    write_csv(data.frame(role = names(inputs), filename = basename(inputs),
      source_path = unname(inputs), md5_before = unname(input_hashes),
      md5_after = after, unchanged = TRUE), "MIMIC_input_audit")
    summary_lines <- c(
      "Status: COMPLETE", "Dataset: MIMIC-IV v3.1", "Formal phenotype: original frozen V295 G3 S01; V318 exact replay passed.",
      "Exposure window: sepsis onset to +24 h. Formal classifications are read directly, never reassigned.",
      "MAKE30 is the stored ORAKLE-style composite of death, new KRT and persistent kidney dysfunction.",
      "Outcome window: strictly after sepsis +24 h; earlier of hospital discharge or 30 days after index ICU admission.",
      sprintf("Full classified n=%d; analysis n=%d; MAKE30 events=%d.", nrow(d), nrow(a), sum(a$make30)),
      sprintf("Excluded: %d outside evaluable window; %d unknown MAKE30 within evaluable window.",
        sum(!eligible), sum(eligible & is.na(d$make30))),
      "Unknown outcomes are excluded, never converted to non-events. Existing outcomes and cohort files are unchanged.",
      "Main model: modified Poisson log-link regression, HC1 robust covariance, C1 reference.",
      paste0("Adjustment: ", paste(CORE, collapse = ", "), "."),
      "All eight main covariates are complete in MIMIC; no covariate imputation or complete-case loss is needed.",
      "The common adjustment formula and robust covariance calculation retain V307 exactly.",
      "Supplementary logistic OR, if enabled, is a separate effect measure and is not interpreted as RR.",
      "Components use their own eligible known-outcome denominators; component events may overlap and cannot be added.",
      "Table1_full_cohort uses all 9804 classified patients; Table1_analysis_cohort uses the MAKE30 analysis population.",
      "Completed regression checkpoints are reused only if all input hashes, formula and analysis cohort match.",
      "Poisson models estimate associations; their fitted values are not a calibrated individual prediction tool.",
      "These observational estimates describe adjusted associations, not causal treatment effects.",
      "Reference: Zou G. Am J Epidemiol. 2004;159:702-706. doi:10.1093/aje/kwh090.",
      "", "Observed MAKE30:",
      vapply(seq_len(nrow(main_rates)), function(i) sprintf("%s: %d/%d, %.2f%% (95%% CI %.2f-%.2f%%).",
        main_rates$phenotype[i], main_rates$events[i], main_rates$known_outcome_n[i],
        main_rates$event_pct[i], main_rates$CI95_low_pct[i], main_rates$CI95_high_pct[i]), character(1L)),
      "", "Adjusted modified Poisson RR:",
      vapply(seq_len(nrow(adj)), function(i) sprintf("%s: RR %.4f (95%% CI %.4f-%.4f), p=%.5g.",
        adj$comparison[i], adj$RR[i], adj$CI95_low[i], adj$CI95_high[i], adj$p_value[i]), character(1L)))
    write_text(summary_lines, "MAKE30_summary")
    write_text(capture.output(utils::sessionInfo()), "MAKE30_session_info")
    write_csv(data.frame(status = "COMPLETE", analysis_n = nrow(a), events = sum(a$make30),
      raw_md5 = input_hashes["raw_cohort"], class_md5 = input_hashes["formal_classification"],
      completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "MAKE30_completion")
    progress_log("MAKE30 analysis completed. All original inputs unchanged.")
    console_result <- adj[c("comparison", "RR", "CI95_low", "CI95_high", "p_value")]
    console_result$p_value <- formatC(console_result$p_value, format = "e", digits = 3L)
    print(console_result, row.names = FALSE)
    list(ok = TRUE, flow = flow, adjusted_RR = adj)
  }, error = function(e) {
    write_text(c("Status: ERROR", conditionMessage(e)), "MAKE30_summary")
    progress_log("ERROR: ", conditionMessage(e))
    list(ok = FALSE, error = conditionMessage(e))
  })

  zip_path <- paste0(out_dir, ".zip")
  zip_ok <- FALSE
  if (.Platform$OS.type == "windows") {
    ps <- Sys.which("powershell.exe")
    if (nzchar(ps)) {
      ps_quote <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
      ps_file <- tempfile("V319_zip_", fileext = ".ps1")
      on.exit(unlink(ps_file), add = TRUE)
      lines <- c("$ErrorActionPreference = 'Stop'",
        paste0("$src = ", ps_quote(normalizePath(out_dir, winslash = "/", mustWork = TRUE))),
        paste0("$dst = ", ps_quote(zip_path)),
        "$items = @(Get-ChildItem -LiteralPath $src -File | Select-Object -ExpandProperty FullName)",
        "Compress-Archive -LiteralPath $items -DestinationPath $dst -CompressionLevel Optimal -Force")
      con <- file(ps_file, "wb")
      writeBin(c(as.raw(c(239, 187, 191)), charToRaw(enc2utf8(paste(lines, collapse = "\r\n")))), con)
      close(con)
      code <- tryCatch(system2(ps, c("-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass",
                                    "-File", shQuote(ps_file))), error = function(e) 1L)
      zip_ok <- identical(as.integer(code), 0L) && file.exists(zip_path)
    }
  } else if (nzchar(Sys.which("zip"))) {
    oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE); setwd(out_dir)
    tryCatch(utils::zip(zip_path, list.files(out_dir, full.names = FALSE), flags = "-q",
                        zip = Sys.which("zip")), error = function(e) NULL)
    setwd(oldwd); zip_ok <- file.exists(zip_path)
  }
  cat("\nOutput directory: ", out_dir, "\n", sep = "")
  if (zip_ok) cat("Result ZIP: ", zip_path, "\n", sep = "") else
    message("ZIP unavailable; all results remain in the output directory.")
  invisible(list(status = status, output_directory = out_dir, zip = if (zip_ok) zip_path else NULL))
}

main_make30()

