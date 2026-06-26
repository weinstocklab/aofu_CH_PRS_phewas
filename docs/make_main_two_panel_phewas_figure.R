#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
})

script_file <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])
if (is.na(script_file) || script_file == "") {
  script_file <- "docs/make_main_two_panel_phewas_figure.R"
}
root <- normalizePath(file.path(dirname(script_file), ".."))
if (!dir.exists(file.path(root, "results"))) {
  root <- getwd()
}
out_dir <- file.path(root, "docs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

format_count <- function(x, noun) {
  paste0(formatC(x, format = "d", big.mark = ","), " ", noun)
}

format_p_plotmath <- function(p) {
  vapply(p, function(value) {
    if (is.na(value)) {
      return("")
    }
    exponent <- floor(log10(value))
    mantissa <- value / (10^exponent)
    if (mantissa >= 9.95) {
      mantissa <- 1
      exponent <- exponent + 1
    }
    mantissa_text <- formatC(mantissa, format = "fg", digits = 2)
    paste0("italic(P)==", mantissa_text, "%*%10^{", exponent, "}")
  }, character(1))
}

shorten_label <- function(x) {
  recode(
    x,
    "Chronic liver disease and cirrhosis" = "Chronic liver disease",
    "Acute upper respiratory infections" = "Respiratory infections",
    .default = x
  )
}

normalize_phecode <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "\\.0$", "")
  x <- if_else(str_detect(x, "^\\d+$") & nchar(x) < 3, str_pad(x, width = 3, pad = "0"), x)
  x
}

download_pheinfo <- function() {
  cache <- file.path(root, "docs", "pheinfo.rda")
  if (!file.exists(cache)) {
    download.file(
      "https://raw.githubusercontent.com/weinstocklab/gwasplot/main/data/pheinfo.rda",
      cache,
      quiet = TRUE,
      mode = "wb"
    )
  }
  env <- new.env(parent = emptyenv())
  load(cache, envir = env)
  env$pheinfo %>%
    mutate(phecode_chr = as.character(phecode)) %>%
    select(phecode_chr, phewas_group = group, phewas_color = color)
}

with_phewas_groups <- function(df, pheinfo) {
  df %>%
    mutate(phecode_chr = normalize_phecode(phecodes)) %>%
    left_join(pheinfo, by = "phecode_chr") %>%
    mutate(
      phewas_group = case_when(
        phenotype_id == "mortality" ~ "mortality",
        str_detect(category %||% "", "heme") ~ "neoplasms",
        TRUE ~ coalesce(phewas_group, category, "other")
      ),
      phewas_color = case_when(
        phenotype_id == "mortality" ~ "gray20",
        str_detect(category %||% "", "heme") ~ "darkcyan",
        TRUE ~ coalesce(phewas_color, "gray40")
      )
    )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

make_panel <- function(
    df,
    effect_col,
    lower_col,
    upper_col,
    count_col,
    count_noun,
    xlab,
    title,
    overlay_df = NULL,
    overlay_effect_col = NULL,
    overlay_lower_col = NULL,
    overlay_upper_col = NULL,
    overlay_count_col = NULL) {
  effect_col <- rlang::ensym(effect_col)
  lower_col <- rlang::ensym(lower_col)
  upper_col <- rlang::ensym(upper_col)
  count_col <- rlang::ensym(count_col)

  plot_df <- df %>%
    filter(status == "ok", !is.na(!!effect_col), !is.na(!!lower_col), !is.na(!!upper_col)) %>%
    mutate(
      effect = !!effect_col,
      lower = !!lower_col,
      upper = !!upper_col,
      count = !!count_col,
      p_label = format_p_plotmath(p),
      label_plot = shorten_label(label),
      label_with_n = paste0(label_plot, " (", format_count(count, count_noun), ")"),
      label_with_n = fct_reorder(label_with_n, effect, .desc = FALSE)
    )

  overlay_plot_df <- NULL
  if (!is.null(overlay_df)) {
    overlay_effect_col <- rlang::ensym(overlay_effect_col)
    overlay_lower_col <- rlang::ensym(overlay_lower_col)
    overlay_upper_col <- rlang::ensym(overlay_upper_col)
    overlay_count_col <- rlang::ensym(overlay_count_col)

    overlay_plot_df <- overlay_df %>%
      filter(status == "ok", !is.na(!!overlay_effect_col), !is.na(!!overlay_lower_col), !is.na(!!overlay_upper_col)) %>%
      transmute(
        phenotype_id,
        overlay_count = !!overlay_count_col,
        overlay_effect = !!overlay_effect_col,
        overlay_lower = !!overlay_lower_col,
        overlay_upper = !!overlay_upper_col,
        overlay_p = p,
        covariate_model = covariate_model %||% NA_character_
      ) %>%
      inner_join(
        plot_df %>%
          mutate(label_with_n = as.character(label_with_n)) %>%
          select(phenotype_id, label_with_n, phewas_group, phewas_color),
        by = "phenotype_id"
      ) %>%
      mutate(label_with_n = factor(label_with_n, levels = levels(plot_df$label_with_n)))
  }

  color_values <- plot_df %>%
    distinct(phewas_group, phewas_color)
  color_values <- stats::setNames(color_values$phewas_color, color_values$phewas_group)

  effect_min <- min(
    c(plot_df$lower, if (!is.null(overlay_plot_df)) overlay_plot_df$overlay_lower else NA_real_),
    na.rm = TRUE
  )
  effect_max <- max(
    c(plot_df$upper, if (!is.null(overlay_plot_df)) overlay_plot_df$overlay_upper else NA_real_),
    na.rm = TRUE
  )
  effect_range <- effect_max - effect_min
  p_x <- effect_max + 0.035 * effect_range
  x_limit_max <- effect_max + 0.22 * effect_range
  main_position <- if (!is.null(overlay_plot_df) && nrow(overlay_plot_df) > 0) {
    position_nudge(y = -0.12)
  } else {
    "identity"
  }
  overlay_position <- position_nudge(y = 0.12)

  p <- ggplot(plot_df, aes(x = effect, y = label_with_n, color = phewas_group)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "gray35") +
    geom_errorbar(
      aes(xmin = lower, xmax = upper),
      orientation = "y",
      width = 0,
      linewidth = 0.55,
      position = main_position
    ) +
    geom_point(aes(shape = "Germline CH PRS (P shown)"), size = 1.8, position = main_position) +
    geom_text(
      aes(x = p_x, label = p_label),
      inherit.aes = FALSE,
      y = plot_df$label_with_n,
      parse = TRUE,
      hjust = 0,
      size = 2.0,
      color = "gray20"
    ) +
    annotate(
      "segment",
      x = 1,
      xend = effect_min,
      y = 0.42,
      yend = 0.42,
      arrow = arrow(length = unit(0.065, "inches")),
      color = "gray35",
      linewidth = 0.3
    ) +
    annotate(
      "segment",
      x = 1,
      xend = effect_max,
      y = 0.42,
      yend = 0.42,
      arrow = arrow(length = unit(0.065, "inches")),
      color = "gray35",
      linewidth = 0.3
    ) +
    annotate("text", x = effect_min, y = 0.08, label = "protective", hjust = 0, size = 2.1, color = "gray25") +
    annotate("text", x = effect_max, y = 0.08, label = "increased risk", hjust = 1, size = 2.1, color = "gray25") +
    scale_color_manual(values = color_values, name = "PheWAS group") +
    scale_shape_manual(
      values = c("Germline CH PRS (P shown)" = 16, "Observed passenger mutation burden" = 17),
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(effect_min, x_limit_max),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_y_discrete(expand = expansion(add = c(1.65, 0.45))) +
    coord_cartesian(clip = "off") +
    labs(title = title, x = xlab, y = NULL) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 8.5),
      plot.margin = margin(4, 12, 4, 20),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 5.9),
      axis.text.x = element_text(size = 6.5),
      axis.title.x = element_text(size = 6.6),
      legend.position = "none"
    )

  if (!is.null(overlay_plot_df) && nrow(overlay_plot_df) > 0) {
    p <- p +
      geom_errorbar(
        data = overlay_plot_df,
        aes(x = overlay_effect, xmin = overlay_lower, xmax = overlay_upper, y = label_with_n, color = phewas_group),
        inherit.aes = FALSE,
        orientation = "y",
        width = 0,
        linewidth = 0.5,
        alpha = 0.35,
        position = overlay_position
      ) +
      geom_point(
        data = overlay_plot_df,
        aes(x = overlay_effect, y = label_with_n, color = phewas_group, shape = "Observed passenger mutation burden"),
        inherit.aes = FALSE,
        size = 1.75,
        alpha = 0.45,
        position = overlay_position
      )
  }

  p
}

make_bottom_legend <- function(group_df) {
  group_df <- group_df %>%
    distinct(phewas_group, phewas_color) %>%
    arrange(phewas_group) %>%
    mutate(
      group_label = str_to_sentence(str_replace_all(phewas_group, "_", " ")),
      group_label = recode(
        group_label,
        "Endocrine/metabolic" = "Endocrine/metabolic",
        "Genitourinary" = "Genitourinary",
        "Mental disorders" = "Mental disorders",
        .default = group_label
      ),
      idx = row_number()
    )

  first_row_n <- ceiling(nrow(group_df) / 2)
  group_df <- group_df %>%
    mutate(
      legend_row = if_else(idx <= first_row_n, 1L, 2L),
      row_index = if_else(legend_row == 1L, idx, idx - first_row_n),
      row_n = if_else(legend_row == 1L, first_row_n, nrow(group_df) - first_row_n),
      x = 0.20 + (row_index - 1) * (0.74 / pmax(row_n - 1, 1)),
      y = if_else(legend_row == 1L, 0.43, 0.18)
    )

  exposure_df <- data.frame(
    x = c(0.23, 0.49),
    y = c(0.78, 0.78),
    label = c("Germline CH PRS", "Observed passenger mutation burden"),
    shape = c(16, 17)
  )

  ggplot() +
    annotate("text", x = 0.04, y = 0.78, label = "Exposure", hjust = 0, size = 2.3, fontface = "bold") +
    geom_point(
      data = exposure_df,
      aes(x = x, y = y, shape = label),
      color = "gray25",
      size = 2.0
    ) +
    geom_text(
      data = exposure_df,
      aes(x = x + 0.025, y = y, label = label),
      hjust = 0,
      size = 2.25,
      color = "gray20"
    ) +
    annotate("text", x = 0.04, y = 0.43, label = "PheWAS category", hjust = 0, size = 2.3, fontface = "bold") +
    geom_point(
      data = group_df,
      aes(x = x, y = y, color = phewas_group),
      size = 2.0
    ) +
    geom_text(
      data = group_df,
      aes(x = x + 0.018, y = y, label = group_label),
      hjust = 0,
      size = 2.0,
      color = "gray20"
    ) +
    scale_shape_manual(values = c("Germline CH PRS" = 16, "Observed passenger mutation burden" = 17)) +
    scale_color_manual(values = stats::setNames(group_df$phewas_color, group_df$phewas_group)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    cowplot::theme_cowplot(font_size = 8) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(0, 10, 0, 10)
    )
}

panel_divider <- ggplot() +
  geom_vline(xintercept = 0.5, color = "gray80", linewidth = 0.45) +
  xlim(0, 1) +
  cowplot::theme_cowplot(font_size = 8) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

pheinfo <- download_pheinfo()

logistic <- read_tsv(
  file.path(root, "results/phewas/PRSFNN_out_final.eur.phewas.tsv"),
  col_types = cols(.default = col_guess(), phecodes = col_character())
) %>%
  with_phewas_groups(pheinfo)

survival <- read_tsv(
  file.path(root, "results/phewas/PRSFNN_out_final.eur.survival_phewas.tsv"),
  show_col_types = FALSE
) %>%
  left_join(logistic %>% select(phenotype_id, phecodes, category), by = "phenotype_id") %>%
  with_phewas_groups(pheinfo)

somatic_survival <- read_tsv(
  file.path(root, "results/phewas/somatic_burden_chr21_pos12mb_ad3.post_draw_survival_phewas.tsv"),
  show_col_types = FALSE
) %>%
  left_join(logistic %>% select(phenotype_id, phecodes, category), by = "phenotype_id") %>%
  with_phewas_groups(pheinfo)

bottom_legend <- make_bottom_legend(bind_rows(logistic, survival, somatic_survival))

logistic_panel <- make_panel(
  logistic,
  effect_col = or,
  lower_col = ci_lower,
  upper_col = ci_upper,
  count_col = cases,
  count_noun = "cases",
  xlab = "OR per 1 s.d. CH PRS",
  title = "a  Logistic regression PheWAS"
)

survival_panel <- make_panel(
  survival,
  effect_col = hr,
  lower_col = ci_lower,
  upper_col = ci_upper,
  count_col = events,
  count_noun = "events",
  xlab = "HR per 1 s.d. exposure",
  title = "b  Survival Analysis PheWAS",
  overlay_df = somatic_survival,
  overlay_effect_col = hr,
  overlay_lower_col = ci_lower,
  overlay_upper_col = ci_upper,
  overlay_count_col = events
)

panel_plot <- logistic_panel + panel_divider + survival_panel +
  plot_layout(ncol = 3, widths = c(1, 0.035, 1)) +
  plot_annotation(
    title = "Germline CH PRS and observed passenger mutation burden across targeted EHR phenotypes",
    theme = theme(plot.title = element_text(face = "bold", size = 10.5))
  )

main_plot <- cowplot::plot_grid(
  cowplot::ggdraw(panel_plot),
  bottom_legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

main_plot <- cowplot::ggdraw(main_plot) +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(
  file.path(out_dir, "main_two_panel_logistic_survival_phewas.pdf"),
  main_plot,
  width = 9,
  height = 9,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(out_dir, "main_two_panel_logistic_survival_phewas.png"),
  main_plot,
  width = 9,
  height = 9,
  units = "in",
  dpi = 360,
  bg = "white"
)

bind_rows(
  logistic %>%
    transmute(
      analysis = "logistic",
      exposure = "Germline CH PRS",
      phenotype_id,
      label,
      phewas_group,
      phewas_color,
      n,
      count = cases,
      effect = or,
      ci_lower,
      ci_upper,
      p,
      p_label = format_p_plotmath(p),
      fdr_bh,
      covariate_model = NA_character_
    ),
  survival %>%
    transmute(
      analysis = "survival_germline_prs",
      exposure = "Germline CH PRS",
      phenotype_id,
      label,
      phewas_group,
      phewas_color,
      n,
      count = events,
      effect = hr,
      ci_lower,
      ci_upper,
      p,
      p_label = format_p_plotmath(p),
      fdr_bh,
      covariate_model = NA_character_
    ),
  somatic_survival %>%
    transmute(
      analysis = "survival_observed_somatic_burden_post_draw",
      exposure = "Observed passenger mutation burden",
      phenotype_id,
      label,
      phewas_group,
      phewas_color,
      n,
      count = events,
      effect = hr,
      ci_lower,
      ci_upper,
      p,
      p_label = format_p_plotmath(p),
      fdr_bh,
      covariate_model = covariate_model %||% NA_character_
    )
) %>%
  write_tsv(file.path(out_dir, "main_two_panel_logistic_survival_source_data.tsv"))
