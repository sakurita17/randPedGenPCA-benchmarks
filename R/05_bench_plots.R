# Generates scalability and memory plots from benchmark_results.rds.
# Usage: Rscript R/05_bench_plots.R [path/to/benchmark_results.rds]
# Without an argument, uses the most recent output directory.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# --- Locate the RDS ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && file.exists(args[1])) {
  rds_path <- args[1]
  out_dir  <- dirname(rds_path)
} else {
  base <- "/Users/s2745345/Documents/GitHub/randPedGenPCA-benchmarks/output"
  dirs <- sort(list.dirs(base, recursive = FALSE, full.names = TRUE), decreasing = TRUE)
  stopifnot("No directories found in output/" = length(dirs) > 0)
  rds_path <- file.path(dirs[1], "benchmark_results.rds")
  out_dir  <- dirs[1]
}
stopifnot("RDS file not found" = file.exists(rds_path))
results_all <- readRDS(rds_path)
message("Reading: ", rds_path)

#results <- results_all #plot all scenarios

results <- results_all %>%
  filter(scenario == "hybrid_two_metafounders")

presentation_theme <- theme_bw(base_size = 26) +
  theme(
    plot.title = element_text(size = 32, face = "bold"),
    plot.subtitle = element_text(size = 26),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 24),
    strip.text = element_text(size = 26, face = "bold"),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 24),
    legend.position = "bottom"
  )

# n_markers label for grouping/faceting (nChr * segSites, always available from params)
results$n_markers_label <- paste0(results$nChr * results$segSites, " markers")

# Memory column: mem_mb (peak RAM from bench) or object_mb (result object size, legacy)
if ("mem_mb" %in% names(results) && any(!is.na(results$mem_mb))) {
  results$mem_plot <- results$mem_mb
  mem_label <- "Peak RAM (MB)"
} else {
  results$mem_plot <- results$object_mb
  mem_label <- "Result object size (MB)"
}

save_plot <- function(p, name, w = 10, h = 6) {
  path <- file.path(out_dir, paste0(name, ".png"))
  ggsave(path, plot = p, width = w, height = h, dpi = 300)
  message("Saved: ", path)
}

fmt_num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

step_colors <- c(
  "gammaL"      = "#E69F00",
  "getLInv"     = "#56B4E9",
  "rhpca"       = "#009E73",
  "getH_prcomp" = "#D55E00"
)

step_labels <- c(
  "gammaL"      = "getGammaL",
  "getLInv"     = "getLInv",
  "rhpca"       = "rhpca",
  "getH_prcomp" = "getH + prcomp (exact)"
)

# ============================================================
# Plot 1: Log-log — time vs n_pedigree (all steps)
# Faceted by scenario (cols) and n_markers (rows).
# ============================================================
p1_data <- results %>%
  filter(
    step %in% c("getLInv", "rhpca", "getH_prcomp", "gammaL"),
    !is.na(elapsed_sec),
    is.na(centred) | centred == "0.5",
    is.na(x_type)  | x_type  == "ped"
  ) %>%
  group_by(step, n_pedigree, scenario, n_markers_label) %>%
  summarise(elapsed_sec = mean(elapsed_sec, na.rm = TRUE), .groups = "drop")

p1 <- ggplot(p1_data, aes(x = n_pedigree, y = elapsed_sec,
                            color = step, group = step)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.5) +
  scale_x_log10(labels = fmt_num) +
  scale_y_log10(labels = function(x) signif(x, 2)) +
  scale_color_manual(values = step_colors, labels = step_labels) +
  facet_grid(n_markers_label ~ scenario) +
  labs(
    title = "Scalability: time vs pedigree size",
    x     = "n_pedigree (log scale)",
    y     = "elapsed_sec (log scale)",
    color = "Step"
  ) +
  presentation_theme

save_plot(p1, "plot1_scalability_time", w = 14, h = 8)

# ============================================================
# Plot 2: rhpca — time vs n_genotyped, color = popSize, linetype = n_markers
# Shows how both pedigree size and marker count affect rhpca time.
# ============================================================
p2_data <- results %>%
  filter(step == "rhpca", !is.na(elapsed_sec), centred == "0.5", x_type == "ped") %>%
  group_by(popSize, n_genotyped, n_markers_label) %>%
  summarise(elapsed_sec = mean(elapsed_sec, na.rm = TRUE), .groups = "drop") %>%
  mutate(popSize = factor(popSize))

p2 <- ggplot(p2_data, aes(x = n_genotyped, y = elapsed_sec,
                            color = popSize,
                            linetype = n_markers_label,
                            group = interaction(popSize, n_markers_label))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.5) +
  scale_x_log10(labels = fmt_num) +
  scale_y_log10(labels = function(x) signif(x, 2)) +
  labs(
    title    = "rhpca: time vs genotyped animals",
    subtitle = "Color = popSize (pedigree size)  |  Line type = number of markers",
    x        = "n_genotyped (log scale)",
    y        = "elapsed_sec (log scale)",
    color    = "popSize",
    linetype = "Markers"
  ) +
  presentation_theme

save_plot(p2, "plot2_rhpca_vs_ngenotyped", w = 12, h = 7)

# ============================================================
# Plot 3: Memory — rhpca vs exact PCA, linetype = n_markers
# getH does not use markers; rhpca memory grows with n_markers.
# ============================================================
p3_data <- results %>%
  filter(
    step %in% c("rhpca", "getH_prcomp"),
    !is.na(mem_plot), mem_plot > 0,
    is.na(centred) | centred == "0.5",
    x_type == "ped"
  ) %>%
  group_by(step, scenario, n_pedigree, n_markers_label) %>%
  summarise(mem_plot = mean(mem_plot, na.rm = TRUE), .groups = "drop")

p3 <- ggplot(p3_data, aes(x = n_pedigree, y = mem_plot,
                            color = step,
                            linetype = n_markers_label,
                            group = interaction(step, n_markers_label))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.5) +
  scale_x_log10(labels = fmt_num) +
  scale_y_log10() +
  scale_color_manual(values = step_colors, labels = step_labels) +
  facet_wrap(~scenario, ncol = 3) +
  labs(
    title    = paste("Memory comparison: rhpca vs exact PCA —", mem_label),
    subtitle = "getH memory is independent of markers; rhpca grows with marker count",
    x        = "n_pedigree (log scale)",
    y        = paste0(mem_label, " (log scale)"),
    color    = "Step",
    linetype = "Markers"
  ) +
  presentation_theme

save_plot(p3, "plot3_memory", w = 12, h = 7)

# ============================================================
# Plot 4: Stacked bars — time breakdown by step
# Faceted by n_markers to show how marker count shifts the bottleneck.
# ============================================================
p4_data <- results %>%
  filter(
    step %in% c("gammaL", "getLInv", "rhpca", "getH_prcomp"),
    !is.na(elapsed_sec),
    is.na(centred) | centred == "0.5",
    is.na(x_type)  | x_type  == "ped",
    scenario == "hybrid_two_metafounders"
  ) %>%
  group_by(step, popSize, n_markers_label) %>%
  summarise(elapsed_sec = mean(elapsed_sec, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    popSize = factor(popSize, levels = sort(unique(popSize))),
    step    = factor(step, levels = c("getH_prcomp", "rhpca", "getLInv", "gammaL"))
  )

p4 <- ggplot(p4_data, aes(x = popSize, y = elapsed_sec, fill = step)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = step_colors, labels = step_labels) +
  facet_wrap(~n_markers_label) +
  labs(
    title    = "Total time breakdown by step",
    subtitle = "Scenario: two metafounders | centred = '0.5' | x = ped",
    x        = "popSize",
    y        = "Cumulative time (s)",
    fill     = "Step"
  ) +
  presentation_theme

save_plot(p4, "plot4_time_breakdown", w = 12, h = 7)

# ============================================================
# Plot 5: Speedup — rhpca vs exact PCA (where both ran)
# n_markers added to join key to avoid many-to-many.
# ============================================================
rhpca_t <- results %>%
  filter(step == "rhpca", x_type == "ped", centred == "0.5", !is.na(elapsed_sec)) %>%
  select(scenario, n_pedigree, geno_pct, n_markers_label, elapsed_rhpca = elapsed_sec)

exact_t <- results %>%
  filter(step == "getH_prcomp", exact_status == "ok", centred == "0.5", !is.na(elapsed_sec)) %>%
  select(scenario, n_pedigree, geno_pct, n_markers_label, elapsed_exact = elapsed_sec)

p5_data <- inner_join(rhpca_t, exact_t,
                       by = c("scenario", "n_pedigree", "geno_pct", "n_markers_label")) %>%
  mutate(
    speedup  = elapsed_exact / elapsed_rhpca,
    geno_pct = paste0(geno_pct, "% genotyped")
  )

if (nrow(p5_data) > 0) {
  p5 <- ggplot(p5_data, aes(x = n_pedigree, y = speedup,
                              color = scenario,
                              linetype = n_markers_label,
                              group = interaction(scenario, n_markers_label))) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3.5) +
    facet_wrap(~geno_pct) +
    labs(
      title    = "Speedup: exact PCA / rhpca",
      subtitle = "Values > 1 mean rhpca is faster",
      x        = "n_pedigree",
      y        = "Times faster",
      color    = "Scenario",
      linetype = "Markers"
    ) +
      presentation_theme
    save_plot(p5, "plot5_speedup", w = 12, h = 7)
} else {
  message("Plot 5 (speedup): no data — needs popSize <= 50 where both exact and rhpca ran.")
}

# ============================================================
# Plot 6: Direct comparison rhpca vs exact PCA.
# Compare the exact same PCA rank (k) used by randomised PCA, and select the
# largest marker series available for both methods.
# ============================================================
if (!"pca_rank" %in% names(results_all)) {
  results_all$pca_rank <- NA_integer_
}

p6_base <- results_all %>%
  filter(
    scenario == "hybrid_no_metafounder",
    centred == "0.5",
    step %in% c("rhpca", "getH_prcomp"),
    !is.na(elapsed_sec)
  )

p6_ranks <- p6_base %>%
  filter(step == "rhpca", x_type == "ped", !is.na(pca_rank)) %>%
  distinct(pca_rank) %>%
  pull(pca_rank)

if (length(p6_ranks) == 0) {
  message("Plot 6: no pca_rank column/data found. Re-run the benchmark to compare exact PCA with the same k as rhpca.")
} else {
  p6_rank <- min(p6_ranks)

  p6_common_markers <- p6_base %>%
    filter(
      pca_rank == p6_rank,
      (step == "rhpca" & x_type == "ped") |
        (step == "getH_prcomp" & (is.na(exact_status) | exact_status == "ok"))
    ) %>%
    distinct(step, n_markers) %>%
    group_by(n_markers) %>%
    summarise(n_methods = n_distinct(step), .groups = "drop") %>%
    filter(n_methods == 2, !is.na(n_markers)) %>%
    pull(n_markers)

  if (length(p6_common_markers) == 0) {
    message("Plot 6: no common marker series found for rhpca and exact PCA at k = ", p6_rank, ".")
  } else {
    p6_markers <- max(p6_common_markers)

    p6_rhpca <- p6_base %>%
      filter(
        step == "rhpca",
        x_type == "ped",
        pca_rank == p6_rank,
        n_markers == p6_markers
      ) %>%
      group_by(n_pedigree) %>%
      summarise(elapsed_sec = mean(elapsed_sec, na.rm = TRUE), .groups = "drop") %>%
      mutate(method = paste0("Randomised PCA (k = ", p6_rank, ")"))

    p6_exact <- p6_base %>%
      filter(
        step == "getH_prcomp",
        pca_rank == p6_rank,
        n_markers == p6_markers,
        is.na(exact_status) | exact_status == "ok"
      ) %>%
      group_by(n_pedigree) %>%
      summarise(elapsed_sec = mean(elapsed_sec, na.rm = TRUE), .groups = "drop") %>%
      mutate(method = paste0("Standard PCA (k = ", p6_rank, ")"))

    p6_common_n <- intersect(p6_rhpca$n_pedigree, p6_exact$n_pedigree)
    p6_rhpca <- p6_rhpca %>% filter(n_pedigree %in% p6_common_n)
    p6_exact <- p6_exact %>% filter(n_pedigree %in% p6_common_n)
    p6_data <- bind_rows(p6_rhpca, p6_exact)

    if (length(p6_common_n) == 0) {
      message("Plot 6: no paired pedigree sizes found for rhpca and exact PCA at k = ", p6_rank, ".")
    } else {
      method_colors <- setNames(
        c("#4CC9F0", "#F72585"),
        c(
          paste0("Randomised PCA (k = ", p6_rank, ")"),
          paste0("Standard PCA (k = ", p6_rank, ")")
        )
      )

      p6 <- ggplot(p6_data, aes(x = n_pedigree, y = elapsed_sec,
                                  color = method,
                                  group = method)) +
        geom_line(linewidth = 1.8) +
        geom_point(size = 6.0, stroke = 0) +
        scale_x_log10(labels = fmt_num) +
        scale_y_log10(labels = function(x) signif(x, 2)) +
        scale_color_manual(values = method_colors) +
        labs(
          subtitle = paste0(fmt_num(p6_markers), " markers; same PCA rank for both methods"),
          x        = "Number of individuals",
          y        = "Execution time (s, log)",
          color    = NULL
        ) +
        theme_bw(base_size = 22, base_family = "Lato") +
        theme(
          text = element_text(family = "Lato"),
          plot.subtitle = element_text(size = 24, color = "black"),
          axis.title = element_text(size = 34, face = "plain"),
          axis.text = element_text(size = 28, color = "black"),
          axis.ticks = element_line(color = "grey35", linewidth = 0.8),
          axis.ticks.length = grid::unit(5, "pt"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.9),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.55),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 27, color = "black"),
          legend.key = element_blank(),
          legend.key.width = grid::unit(1.5, "cm"),
          legend.spacing.x = grid::unit(0.6, "cm"),
          plot.margin = margin(10, 12, 10, 10)
        )

      save_plot(p6, "plot6_direct_comparison", w = 12, h = 7)
    }
  }
}

message("\nAll plots saved in: ", out_dir)
