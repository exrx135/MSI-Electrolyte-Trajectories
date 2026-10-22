#!/usr/bin/env Rscript

# ==============================================================================
# V297_joint_portrait_G3G6.R
# FOURTH STUDY / MSI-Electrolyte GBMTM
#
# Purpose
#   Create an early visual preview of the joint phenotypes for G3-G6 from the
#   completed Stage-3 anchor table. This script does not refit or modify models.
#
# Main outputs for each G
#   1) Joint standardized portrait:
#      - one column per class;
#      - five indicators overlaid in each class panel;
#      - upper row = absolute fitted z values (level + shape);
#      - lower row = baseline-centered fitted z values (shape emphasis).
#
#   2) Native-scale phenotype portrait:
#      - one column per class and one row per indicator;
#      - values remain in their clinical units;
#      - y scales are shared across classes within each indicator row.
#
# Source
#   <USERPROFILE>/Desktop/Research/MSI/提取轨迹文件/
#     GBMTM_stage3_output/tables/V294_class_anchor_G2G20.csv
#
# Output
#   <USERPROFILE>/Desktop/Research/MSI/提取轨迹文件/
#     GBMTM_stage3_output/V297_joint_portrait_G3G6
#
# Important
#   - Uses only the real analysis anchors already present in V294:
#       MSI: 1, 3, 5, ..., 23 h
#       Na/K/Cl/HCO3: 3, 9, 15, 21 h
#   - Does not interpolate or manufacture additional electrolyte time points.
#   - Does not alter Stage 2 or Stage 3 results.
# ==============================================================================

options(stringsAsFactors = FALSE)

VERSION_TAG <- "V297"
TARGET_G <- 3:6

# ---- Resolve Windows user profile without hardcoding a username --------------
resolve_user_profile <- function() {
  user_profile <- Sys.getenv("USERPROFILE", unset = "")
  user_profile <- trimws(user_profile)

  if (nzchar(user_profile) && dir.exists(user_profile)) {
    return(normalizePath(user_profile, winslash = "/", mustWork = TRUE))
  }

  expanded_home <- path.expand("~")
  expanded_home <- normalizePath(expanded_home, winslash = "/", mustWork = FALSE)

  # Some Windows/R configurations resolve ~ to the Documents folder.
  if (tolower(basename(expanded_home)) == "documents") {
    expanded_home <- dirname(expanded_home)
  }

  if (!dir.exists(expanded_home)) {
    stop(
      "Unable to resolve the Windows user profile. ",
      "USERPROFILE is empty/invalid and path.expand('~') is not usable.",
      call. = FALSE
    )
  }

  normalizePath(expanded_home, winslash = "/", mustWork = TRUE)
}

USER_ROOT <- resolve_user_profile()
PROJECT_ROOT <- file.path(
  USER_ROOT, "Desktop", "Research", "MSI", "提取轨迹文件"
)
STAGE3_ROOT <- file.path(PROJECT_ROOT, "GBMTM_stage3_output")
SOURCE_FILE <- file.path(
  STAGE3_ROOT, "tables", "V294_class_anchor_G2G20.csv"
)
OUTPUT_DIR <- file.path(STAGE3_ROOT, "V297_joint_portrait_G3G6")

# ---- Basic checks -------------------------------------------------------------
required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Required R package(s) not installed: ",
    paste(missing_packages, collapse = ", "),
    ". Install with: install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

if (!file.exists(SOURCE_FILE)) {
  stop(
    paste0(
      "Required Stage-3 anchor table was not found.\n",
      "Checked exactly:\n", SOURCE_FILE, "\n\n",
      "Expected project hierarchy:\n",
      "<USERPROFILE>/Desktop/Research/MSI/提取轨迹文件/",
      "GBMTM_stage3_output/tables/V294_class_anchor_G2G20.csv"
    ),
    call. = FALSE
  )
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUTPUT_DIR)) {
  stop("Could not create output directory: ", OUTPUT_DIR, call. = FALSE)
}

LOG_FILE <- file.path(OUTPUT_DIR, "V297_run_log_G3G6.txt")
log_lines <- character()

log_message <- function(...) {
  txt <- paste0(...)
  stamped <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", txt)
  log_lines <<- c(log_lines, stamped)
  message(stamped)
}

write_log <- function() {
  writeLines(log_lines, con = LOG_FILE, useBytes = TRUE)
}

on.exit(write_log(), add = TRUE)

log_message("Script: ", VERSION_TAG, "_joint_portrait_G3G6.R")
log_message("R version: ", R.version.string)
log_message("Resolved user root: ", USER_ROOT)
log_message("Source file: ", SOURCE_FILE)
log_message("Output directory: ", OUTPUT_DIR)
log_message("Requested G values: ", paste(TARGET_G, collapse = ", "))

# ---- Read the single authoritative Stage-3 table ------------------------------
anchor_all <- data.table::fread(
  SOURCE_FILE,
  encoding = "UTF-8",
  showProgress = FALSE
)

required_columns <- c(
  "requested_g",
  "class",
  "class_label",
  "n",
  "class_pct",
  "APPA_model",
  "indicator",
  "time_hour_center",
  "fitted_z",
  "fitted_native",
  "clinical_centered_z",
  "clinical_centered_native",
  "clinical_anchor_type"
)

missing_columns <- setdiff(required_columns, names(anchor_all))
if (length(missing_columns) > 0L) {
  stop(
    "The V294 anchor table is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

anchor <- anchor_all[requested_g %in% TARGET_G]
if (nrow(anchor) == 0L) {
  stop(
    "The source table contains no rows for G3-G6.",
    call. = FALSE
  )
}

# Explicit numeric conversion also catches malformed imported values.
numeric_columns <- c(
  "requested_g", "class", "n", "class_pct", "APPA_model",
  "time_hour_center", "fitted_z", "fitted_native",
  "clinical_centered_z", "clinical_centered_native"
)
anchor[, (numeric_columns) := lapply(.SD, as.numeric), .SDcols = numeric_columns]

if (anyNA(anchor[, ..numeric_columns])) {
  bad_counts <- vapply(
    numeric_columns,
    function(x) sum(is.na(anchor[[x]])),
    integer(1)
  )
  bad_counts <- bad_counts[bad_counts > 0L]
  stop(
    "Missing or non-numeric values detected after import: ",
    paste(names(bad_counts), bad_counts, sep = "=", collapse = ", "),
    call. = FALSE
  )
}

finite_columns <- c(
  "fitted_z", "fitted_native", "clinical_centered_z",
  "clinical_centered_native", "class_pct", "APPA_model"
)
nonfinite_counts <- vapply(
  finite_columns,
  function(x) sum(!is.finite(anchor[[x]])),
  integer(1)
)
nonfinite_counts <- nonfinite_counts[nonfinite_counts > 0L]
if (length(nonfinite_counts) > 0L) {
  stop(
    "Non-finite values detected: ",
    paste(names(nonfinite_counts), nonfinite_counts, sep = "=", collapse = ", "),
    call. = FALSE
  )
}

# ---- Structural validation: G, classes, indicators, and real anchors ----------
indicator_order <- c("MSI", "Na", "K", "Cl", "HCO3")
unknown_indicators <- setdiff(unique(anchor$indicator), indicator_order)
missing_indicators <- setdiff(indicator_order, unique(anchor$indicator))

if (length(unknown_indicators) > 0L || length(missing_indicators) > 0L) {
  stop(
    paste0(
      "Unexpected indicator structure. ",
      "Missing: ", paste(missing_indicators, collapse = ", "),
      "; unknown: ", paste(unknown_indicators, collapse = ", ")
    ),
    call. = FALSE
  )
}

expected_times <- list(
  MSI = seq(1, 23, by = 2),
  Na = c(3, 9, 15, 21),
  K = c(3, 9, 15, 21),
  Cl = c(3, 9, 15, 21),
  HCO3 = c(3, 9, 15, 21)
)

validation_errors <- character()

for (g in TARGET_G) {
  g_data <- anchor[requested_g == g]

  if (nrow(g_data) == 0L) {
    validation_errors <- c(
      validation_errors,
      sprintf("G%d is absent from the source table", g)
    )
    next
  }

  observed_classes <- sort(unique(g_data$class))
  expected_classes <- seq_len(g)
  if (!identical(observed_classes, as.numeric(expected_classes))) {
    validation_errors <- c(
      validation_errors,
      sprintf(
        "G%d classes are {%s}, expected {%s}",
        g,
        paste(observed_classes, collapse = ","),
        paste(expected_classes, collapse = ",")
      )
    )
  }

  for (cl in expected_classes) {
    for (ind in indicator_order) {
      observed_times <- sort(unique(
        g_data[class == cl & indicator == ind, time_hour_center]
      ))
      target_times <- expected_times[[ind]]

      if (!identical(as.numeric(observed_times), as.numeric(target_times))) {
        validation_errors <- c(
          validation_errors,
          sprintf(
            "G%d C%d %s anchors are {%s}, expected {%s}",
            g,
            cl,
            ind,
            paste(observed_times, collapse = ","),
            paste(target_times, collapse = ",")
          )
        )
      }
    }
  }
}

if (length(validation_errors) > 0L) {
  stop(
    paste(
      c(
        "Stage-3 anchor validation failed; no figures were created:",
        paste0("- ", validation_errors)
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

# Verify that class-level metadata is internally constant for every G/class.
metadata_check <- anchor[, .(
  n_n = data.table::uniqueN(n),
  n_pct = data.table::uniqueN(class_pct),
  n_appa = data.table::uniqueN(APPA_model)
), by = .(requested_g, class)]

if (metadata_check[n_n != 1L | n_pct != 1L | n_appa != 1L, .N] > 0L) {
  stop(
    "n, class_pct, or APPA_model is inconsistent within at least one G/class.",
    call. = FALSE
  )
}

data.table::setorder(anchor, requested_g, class, indicator, time_hour_center)
anchor[, indicator := factor(indicator, levels = indicator_order, ordered = TRUE)]

log_message("Input validation passed.")
log_message("Rows retained for G3-G6: ", nrow(anchor))

# ---- Fixed visual encodings ---------------------------------------------------
indicator_colors <- c(
  MSI = "#D55E00",
  Na = "#0072B2",
  K = "#E69F00",
  Cl = "#009E73",
  HCO3 = "#CC79A7"
)

indicator_shapes <- c(
  MSI = 16,
  Na = 17,
  K = 15,
  Cl = 18,
  HCO3 = 8
)

class_colors <- c(
  `1` = "#E64B35",
  `2` = "#4DBBD5",
  `3` = "#00A087",
  `4` = "#3C5488",
  `5` = "#F39B7F",
  `6` = "#8491B4"
)

x_breaks <- c(1, 3, 9, 15, 21, 23)

class_number <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

add_class_panel_label <- function(g_data) {
  class_meta <- unique(g_data[, .(
    class,
    n,
    class_pct,
    APPA_model
  )])
  data.table::setorder(class_meta, class)

  class_meta[, class_panel := sprintf(
    "C%d\nn=%s (%.2f%%) | APPA=%.3f",
    as.integer(class),
    class_number(n),
    class_pct,
    APPA_model
  )]

  g_data <- merge(
    g_data,
    class_meta[, .(class, class_panel)],
    by = "class",
    all.x = TRUE,
    sort = FALSE
  )

  g_data[, class_panel := factor(
    class_panel,
    levels = class_meta$class_panel,
    ordered = TRUE
  )]
  g_data[, class_key := factor(
    as.character(as.integer(class)),
    levels = as.character(seq_len(max(class))),
    ordered = TRUE
  )]

  data.table::setorder(g_data, class, indicator, time_hour_center)
  g_data
}

# ---- Figure 1: joint standardized overlay ------------------------------------
make_joint_overlay <- function(g) {
  g_data <- add_class_panel_label(anchor[requested_g == g])

  absolute <- g_data[, .(
    requested_g,
    class,
    class_panel,
    indicator,
    time_hour_center,
    value = fitted_z,
    view = "Absolute standardized level (level + shape)"
  )]

  centered <- g_data[, .(
    requested_g,
    class,
    class_panel,
    indicator,
    time_hour_center,
    value = clinical_centered_z,
    view = "Baseline-centered temporal shape"
  )]

  plot_data <- data.table::rbindlist(list(absolute, centered), use.names = TRUE)
  plot_data[, view := factor(
    view,
    levels = c(
      "Absolute standardized level (level + shape)",
      "Baseline-centered temporal shape"
    ),
    ordered = TRUE
  )]

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = time_hour_center,
      y = value,
      color = indicator,
      shape = indicator,
      group = indicator
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      linetype = "dashed",
      color = "grey55"
    ) +
    ggplot2::geom_line(linewidth = 0.9, alpha = 0.96) +
    ggplot2::geom_point(size = 2.15, stroke = 0.25) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(view),
      cols = ggplot2::vars(class_panel),
      scales = "free_y",
      switch = "y"
    ) +
    ggplot2::scale_color_manual(
      values = indicator_colors,
      breaks = indicator_order,
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::scale_shape_manual(
      values = indicator_shapes,
      breaks = indicator_order,
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      limits = c(1, 23),
      minor_breaks = NULL,
      expand = ggplot2::expansion(mult = c(0.025, 0.025))
    ) +
    ggplot2::labs(
      title = sprintf("G%d joint MSI-electrolyte phenotype portrait", g),
      subtitle = paste0(
        "Top: absolute fitted z values show level + shape.  ",
        "Bottom: each indicator is centered at its first real anchor to emphasize temporal shape."
      ),
      x = "Hours from sepsis onset",
      y = "Fitted standardized value (z)",
      caption = paste0(
        "Real anchors only: MSI = 1, 3, ..., 23 h; ",
        "Na/K/Cl/HCO3 = 3, 9, 15, 21 h. ",
        "Centered baseline: MSI at 1 h; electrolytes at 3 h."
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey25"),
      plot.caption = ggplot2::element_text(size = 8.5, color = "grey35", hjust = 0),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey55"),
      strip.text.x = ggplot2::element_text(face = "bold", size = 9.2),
      strip.text.y.left = ggplot2::element_text(face = "bold", size = 9.2),
      strip.placement = "outside",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, color = "grey88"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.key.width = grid::unit(1.2, "cm"),
      axis.text.x = ggplot2::element_text(size = 8.5),
      panel.spacing.x = grid::unit(0.55, "lines"),
      panel.spacing.y = grid::unit(0.8, "lines"),
      plot.margin = ggplot2::margin(10, 12, 8, 10)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    )
}

# ---- Figure 2: native-scale class-by-indicator matrix -------------------------
make_native_portrait <- function(g) {
  g_data <- add_class_panel_label(anchor[requested_g == g])

  ggplot2::ggplot(
    g_data,
    ggplot2::aes(
      x = time_hour_center,
      y = fitted_native,
      color = class_key,
      group = class_key
    )
  ) +
    ggplot2::geom_line(linewidth = 1.05, alpha = 0.98) +
    ggplot2::geom_point(size = 2.25, alpha = 0.98) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(indicator),
      cols = ggplot2::vars(class_panel),
      scales = "free_y",
      switch = "y"
    ) +
    ggplot2::scale_color_manual(
      values = class_colors,
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      limits = c(1, 23),
      minor_breaks = NULL,
      expand = ggplot2::expansion(mult = c(0.025, 0.025))
    ) +
    ggplot2::labs(
      title = sprintf("G%d native-scale joint phenotype portrait", g),
      subtitle = paste0(
        "Each column is one phenotype; each row is one indicator. ",
        "The y scale is shared across classes within an indicator row, allowing direct clinical-level comparison."
      ),
      x = "Hours from sepsis onset",
      y = "Fitted value (indicator-specific clinical unit)",
      caption = paste0(
        "Real anchors only: MSI = 1, 3, ..., 23 h; ",
        "Na/K/Cl/HCO3 = 3, 9, 15, 21 h. No added electrolyte time points."
      )
    ) +
    ggplot2::theme_bw(base_size = 10.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey25"),
      plot.caption = ggplot2::element_text(size = 8.5, color = "grey35", hjust = 0),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey55"),
      strip.text.x = ggplot2::element_text(face = "bold", size = 9),
      strip.text.y.left = ggplot2::element_text(face = "bold", size = 10),
      strip.placement = "outside",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, color = "grey88"),
      axis.text.x = ggplot2::element_text(size = 8.2),
      panel.spacing.x = grid::unit(0.48, "lines"),
      panel.spacing.y = grid::unit(0.55, "lines"),
      plot.margin = ggplot2::margin(10, 12, 8, 10)
    )
}

# ---- Generate all per-G PNG files and multi-page PDFs -------------------------
plot_list <- vector("list", length(TARGET_G))
names(plot_list) <- paste0("G", TARGET_G)

manifest <- data.table::data.table(
  version = character(),
  requested_g = integer(),
  figure_type = character(),
  file_name = character(),
  full_path = character(),
  width_in = numeric(),
  height_in = numeric(),
  dpi = integer()
)

for (g in TARGET_G) {
  log_message("Creating G", g, " figures...")

  joint_plot <- make_joint_overlay(g)
  native_plot <- make_native_portrait(g)

  joint_width <- max(11.5, 2.50 * g + 3.2)
  native_width <- max(11.5, 2.35 * g + 3.2)
  joint_height <- 8.2
  native_height <- 11.2

  joint_file <- file.path(
    OUTPUT_DIR,
    sprintf("V297_joint_overlay_G%d.png", g)
  )
  native_file <- file.path(
    OUTPUT_DIR,
    sprintf("V297_native_portrait_G%d.png", g)
  )

  ggplot2::ggsave(
    filename = joint_file,
    plot = joint_plot,
    width = joint_width,
    height = joint_height,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    filename = native_file,
    plot = native_plot,
    width = native_width,
    height = native_height,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )

  manifest <- data.table::rbindlist(
    list(
      manifest,
      data.table::data.table(
        version = VERSION_TAG,
        requested_g = g,
        figure_type = c("joint_standardized_overlay", "native_scale_portrait"),
        file_name = basename(c(joint_file, native_file)),
        full_path = normalizePath(
          c(joint_file, native_file),
          winslash = "/",
          mustWork = TRUE
        ),
        width_in = c(joint_width, native_width),
        height_in = c(joint_height, native_height),
        dpi = c(300L, 300L)
      )
    ),
    use.names = TRUE
  )

  plot_list[[paste0("G", g)]] <- list(
    joint = joint_plot,
    native = native_plot
  )
}

joint_pdf <- file.path(OUTPUT_DIR, "V297_joint_overlay_G3G6.pdf")
grDevices::pdf(
  joint_pdf,
  width = 18,
  height = 8.5,
  onefile = TRUE,
  useDingbats = FALSE
)
for (g in TARGET_G) {
  print(plot_list[[paste0("G", g)]]$joint)
}
grDevices::dev.off()

native_pdf <- file.path(OUTPUT_DIR, "V297_native_portrait_G3G6.pdf")
grDevices::pdf(
  native_pdf,
  width = 18,
  height = 11,
  onefile = TRUE,
  useDingbats = FALSE
)
for (g in TARGET_G) {
  print(plot_list[[paste0("G", g)]]$native)
}
grDevices::dev.off()

atlas_pdf <- file.path(OUTPUT_DIR, "V297_joint_atlas_G3G6.pdf")
grDevices::pdf(
  atlas_pdf,
  width = 18,
  height = 11,
  onefile = TRUE,
  useDingbats = FALSE
)
for (g in TARGET_G) {
  print(plot_list[[paste0("G", g)]]$joint)
  print(plot_list[[paste0("G", g)]]$native)
}
grDevices::dev.off()

pdf_rows <- data.table::data.table(
  version = VERSION_TAG,
  requested_g = NA_integer_,
  figure_type = c(
    "joint_standardized_overlay_multipage_pdf",
    "native_scale_portrait_multipage_pdf",
    "joint_atlas_alternating_multipage_pdf"
  ),
  file_name = basename(c(joint_pdf, native_pdf, atlas_pdf)),
  full_path = normalizePath(
    c(joint_pdf, native_pdf, atlas_pdf),
    winslash = "/",
    mustWork = TRUE
  ),
  width_in = c(18, 18, 18),
  height_in = c(8.5, 11, 11),
  dpi = NA_integer_
)
manifest <- data.table::rbindlist(list(manifest, pdf_rows), use.names = TRUE)

# ---- Export exact plotting data, manifest, README, and session information ----
plot_data_file <- file.path(OUTPUT_DIR, "V297_plot_data_G3G6.csv")
manifest_file <- file.path(OUTPUT_DIR, "V297_figure_manifest_G3G6.csv")
readme_file <- file.path(OUTPUT_DIR, "V297_README_G3G6.txt")
session_file <- file.path(OUTPUT_DIR, "V297_sessionInfo_G3G6.txt")

data.table::fwrite(anchor, plot_data_file, bom = TRUE)
data.table::fwrite(manifest, manifest_file, bom = TRUE)

readme_lines <- c(
  "V297: joint MSI-electrolyte phenotype portraits for G3-G6",
  "",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Source: ", SOURCE_FILE),
  paste0("Output: ", OUTPUT_DIR),
  "",
  "WHAT THIS SCRIPT DOES",
  "1. Reads the completed V294 Stage-3 class-anchor table.",
  "2. Strictly verifies G3-G6, class counts, five indicators, and real anchors.",
  "3. Creates two complementary figure types for every G.",
  "4. Does not refit, overwrite, or modify any Stage-2/Stage-3 model.",
  "",
  "FIGURE TYPE 1: V297_joint_overlay_Gx.png",
  "- One column is one class/phenotype.",
  "- Five colored trajectories are overlaid because all use standardized z units.",
  "- Top row retains absolute level plus temporal shape.",
  "- Bottom row removes each indicator's fitted first-anchor level and emphasizes shape.",
  "- If classes separate above but substantially overlap below, their distinction is",
  "  mainly persistent level rather than temporal shape.",
  "",
  "FIGURE TYPE 2: V297_native_portrait_Gx.png",
  "- One column is one class/phenotype; one row is one indicator.",
  "- The five indicators are not overlaid because they have different clinical units.",
  "- The y scale is shared across class columns within each indicator row.",
  "- Use this figure to judge the clinical magnitude of level and change differences.",
  "",
  "REAL TIME ANCHORS",
  "- MSI: 1, 3, 5, ..., 23 hours.",
  "- Na, K, Cl, HCO3: 3, 9, 15, 21 hours.",
  "- No extra electrolyte time points are interpolated or displayed.",
  "",
  "MULTI-PAGE PDF FILES",
  "- V297_joint_overlay_G3G6.pdf: four joint standardized pages, G3 through G6.",
  "- V297_native_portrait_G3G6.pdf: four native-scale pages, G3 through G6.",
  "- V297_joint_atlas_G3G6.pdf: eight pages, alternating the two views for each G.",
  "",
  "INTERPRETATION RULE",
  "The standardized overlay is for joint pattern/coupling; the native portrait is for",
  "clinical magnitude. Neither figure alone selects the final G. They should be read",
  "together with class size, APPA, stability, and the Stage-3 split/refinement evidence."
)
writeLines(readme_lines, con = readme_file, useBytes = TRUE)
capture.output(sessionInfo(), file = session_file)

log_message("Created per-G PNG files for G3-G6.")
log_message("Created joint PDF: ", joint_pdf)
log_message("Created native PDF: ", native_pdf)
log_message("Created combined atlas PDF: ", atlas_pdf)
log_message("Created plotting-data CSV: ", plot_data_file)
log_message("Created figure manifest: ", manifest_file)
log_message("Completed successfully.")

write_log()

cat(
  paste0(
    "\n",
    "============================================================\n",
    "V297 completed successfully.\n",
    "Output directory:\n",
    normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE), "\n",
    "\nRecommended first file to open:\n",
    normalizePath(atlas_pdf, winslash = "/", mustWork = TRUE), "\n",
    "============================================================\n"
  )
)

