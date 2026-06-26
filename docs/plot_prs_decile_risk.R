#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

script_file <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])
if (is.na(script_file) || script_file == "") {
  script_file <- "docs/plot_prs_decile_risk.R"
}
root <- normalizePath(file.path(dirname(script_file), ".."))
if (!dir.exists(file.path(root, "results"))) {
  root <- getwd()
}
out_dir <- file.path(root, "docs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ci_path <- file.path(root, "results/phewas/supplementary/prs_decile_risk_ci_obesity_alzheimers.tsv")
fallback_path <- file.path(root, "results/phewas/supplementary/prs_decile_risk_obesity_alzheimers.tsv")
input_path <- if (file.exists(ci_path)) ci_path else fallback_path

deciles <- read_tsv(input_path, show_col_types = FALSE) %>%
  mutate(
    decile = as.integer(decile),
    label = recode(label, "Alzheimer's disease" = "Alzheimer's disease"),
    observed_percent = 100 * observed_risk,
    adjusted_percent = 100 * adjusted_risk,
    observed_ci_lower_percent = 100 * coalesce(observed_ci_lower, observed_risk),
    observed_ci_upper_percent = 100 * coalesce(observed_ci_upper, observed_risk),
    adjusted_ci_lower_percent = 100 * coalesce(adjusted_ci_lower, adjusted_risk),
    adjusted_ci_upper_percent = 100 * coalesce(adjusted_ci_upper, adjusted_risk)
  )

plot_df <- bind_rows(
  deciles %>%
    transmute(
      phenotype_id, label, decile, n, cases,
      risk_type = "Observed",
      risk_percent = observed_percent,
      risk_lower = observed_ci_lower_percent,
      risk_upper = observed_ci_upper_percent
    ),
  deciles %>%
    transmute(
      phenotype_id, label, decile, n, cases,
      risk_type = "Adjusted",
      risk_percent = adjusted_percent,
      risk_lower = adjusted_ci_lower_percent,
      risk_upper = adjusted_ci_upper_percent
    )
)

p <- ggplot(plot_df, aes(x = decile, y = risk_percent, color = risk_type, group = risk_type)) +
  geom_ribbon(aes(ymin = risk_lower, ymax = risk_upper, fill = risk_type), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.8) +
  facet_wrap(~ label, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = c(Observed = "gray45", Adjusted = "#007C89")) +
  scale_fill_manual(values = c(Observed = "gray45", Adjusted = "#007C89")) +
  labs(
    title = "Risk by CH PRS decile with 95% confidence intervals",
    x = "CH PRS decile",
    y = "Risk (%)",
    color = NULL
  ) +
  guides(fill = "none") +
  cowplot::theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.25),
    panel.grid.major.x = element_blank()
  )

ggsave(
  file.path(out_dir, "prs_decile_risk_obesity_alzheimers.pdf"),
  p,
  width = 7.2,
  height = 3.6,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(out_dir, "prs_decile_risk_obesity_alzheimers.png"),
  p,
  width = 7.2,
  height = 3.6,
  units = "in",
  dpi = 360
)
