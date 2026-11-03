# Export the summary CSVs to Excel and scientific plots. Run after summarise_eddie.py.
args <- commandArgs(TRUE)
folder <- if (length(args)) args[1] else "reports/eddie_comparison"
english <- length(args) < 2L || args[2] == "en"
tr <- function(es, en) if (english) en else es
stopifnot(dir.exists(folder))
tables_dir <- if (dir.exists(file.path(folder, "data"))) file.path(folder, "data") else folder
figures_dir <- file.path(folder, "figures")
dir.create(figures_dir, showWarnings = FALSE)
library(openxlsx)
library(ggplot2)
tables <- c("escalabilidad", "precision_resumen", "genotipados", "precision_semillas",
            "referencia_exacta", "inventario", "fallos", "duplicados", "auditoria_anidacion")
sheet_names <- c("scalability", "accuracy_summary", "genotyped_animals", "accuracy_by_seed",
                 "exact_reference", "run_inventory", "failures", "duplicates", "nested_subset_audit")
headers <- c(
  fuentes_id="source_signature", genotipados="genotyped_animals", marcadores="markers",
  operador="operator", semilla="seed", tiempo_pca_s="pca_time_s",
  tiempo_pca_min="pca_time_min", tiempo_pca_h="pca_time_h", proceso_pca_s="pca_process_time_s",
  rss_proceso_mib="process_peak_rss_mib", rss_proceso_gib="process_peak_rss_gib",
  diagnostico_s="diagnostic_time_s", residuo_max_relativo="max_relative_residual",
  cambio_eigenvalores_vs_base="eigenvalue_change_vs_base",
  cambio_subespacio_vs_base="subspace_change_vs_base", error_ortogonalidad="orthogonality_error",
  advertencias="warnings", diagnostico_estado="diagnostic_status", archivo_resultados="results_file",
  porcentaje_genotipado="genotyped_percentage", error_eigenvalores_pct="eigenvalue_error_pct",
  error_subespacio="subspace_error", error_scores_alineados="aligned_scores_error",
  correlacion_pc1="abs_pc1_correlation", correlacion_pc2="abs_pc2_correlation",
  error_traza_pct="trace_error_pct", error_proporcion_variacion="variance_proportion_error",
  semillas="pca_seeds", tiempo_pca_mediana_s="median_pca_time_s",
  tiempo_pca_min_s="min_pca_time_s", tiempo_pca_max_s="max_pca_time_s",
  referencia_H_mas_eigen_s="reference_H_plus_eigen_s", rss_pca_max_mib="max_pca_process_rss_mib",
  error_eigenvalores_peor_pct="worst_eigenvalue_error_pct", error_subespacio_peor="worst_subspace_error",
  error_scores_alineados_peor="worst_aligned_scores_error", residuo_peor="worst_relative_residual",
  correlacion_pc1_min="min_abs_pc1_correlation", correlacion_pc2_min="min_abs_pc2_correlation",
  error_traza_abs_peor_pct="worst_abs_trace_error_pct",
  error_proporcion_variacion_peor="worst_variance_proportion_error",
  construir_H_s="build_H_time_s", eigen_exacto_s="exact_eigen_time_s",
  H_mas_eigen_s="H_plus_eigen_time_s", proceso_referencia_s="reference_process_time_s",
  rss_referencia_mib="reference_peak_rss_mib", eigenvalor_minimo="min_eigenvalue",
  separacion_relativa_PC5_PC6="relative_PC5_PC6_eigengap", experimento="experiment",
  pca_esperadas="expected_pca_runs", pca_ok="successful_pca_runs",
  diagnosticos_ok="successful_diagnostics", etapas_con_fallo="failed_stages", cobertura="coverage",
  numero_ids="animal_id_count", anidacion_correcta="nested_subsets_verified", nota="note"
)
wb <- createWorkbook()
for (i in seq_along(tables)) {
  name <- tables[i]
  data <- read.csv(file.path(tables_dir, paste0(name, ".csv")), check.names = FALSE,
                   na.strings = c("", "NA"))
  sheet <- if (english) sheet_names[i] else name
  if (english) {
    if ("cobertura" %in% names(data)) data$cobertura <- ifelse(data$cobertura == "completa", "complete", "incomplete")
    if ("anidacion_correcta" %in% names(data)) data$anidacion_correcta <- as.character(data$anidacion_correcta)
    translated <- unname(headers[names(data)])
    names(data) <- ifelse(is.na(translated), names(data), translated)
  }
  addWorksheet(wb, sheet)
  if (nrow(data)) {
    writeDataTable(wb, sheet, data, tableStyle = "TableStyleMedium2")
    freezePane(wb, sheet, firstRow = TRUE)
    widths <- vapply(seq_along(data), function(i) {
      x <- data[[i]]
      min(42, max(14, nchar(names(data)[i]), nchar(as.character(x[!is.na(x)]))))
    }, 1.0)
    setColWidths(wb, sheet, cols = seq_along(data), widths = widths)
  } else {
    writeData(wb, sheet, tr("Sin registros", "No records"))
  }
}
notes <- c(
  "Leer comparacion.md para interpretación y límites.",
  "Tiempo PCA = cálculo interno; proceso_pca_s incluye arranque/carga/guardado.",
  "RSS = máximo del proceso completo, no memoria incremental de PCA.",
  "Precisión: mediana de tiempos y peor error entre tres semillas de PCA, sobre un dataset.",
  "Errores de eigenvalores y traza en columnas *_pct están expresados en porcentaje.",
  "Error de subespacio y residuos son adimensionales; no son porcentajes de error de eigenvalores.",
  "Los residuos de poblaciones grandes no se comparan con una referencia exacta.",
  "En genotipados se usa una sola semilla; no se puede concluir una ley de rendimiento.",
  "Sin datos descargados de CHC con 100.000: no se han inventado valores faltantes.",
  "Los registros fallidos y ejecuciones incompletas se muestran aparte.",
  "No se han descargado coordenadas .rds: estos gráficos no son PC1 frente a PC2."
)
if (english) notes <- c(
  "Read summary_en.html or data/comparison_en.md for interpretation and limitations.",
  "PCA time is internal computation; process time includes startup, loading and saving.",
  "RSS is peak resident memory of the entire process, not incremental PCA-only memory.",
  "Accuracy summary: median runtime and worst error across three PCA seeds on one dataset.",
  "Eigenvalue and trace error columns ending in _pct are expressed as percentages.",
  "Subspace errors and residuals are dimensionless; they are not eigenvalue error percentages.",
  "Large-population residuals are not comparisons against an exact reference.",
  "The genotyped-animal experiment uses one PCA seed; a general performance law cannot be inferred.",
  "No downloaded CHC results at 100,000 individuals; missing values have not been invented.",
  "Failed stages and incomplete run coverage are listed separately.",
  "Coordinate .rds files were not downloaded; the report figures are not PC1 versus PC2 plots."
)
if (english) notes <- c(notes,
  "The algorithm uses X=U*diag(d). Historical aligned-score errors used a sqrt(d) reference and are not valid coordinate errors; see score_comparison_status.",
  "Historical raw score errors are retained for audit. Corrected coordinate diagnostics require the saved RDS files; PCA runtimes/eigenvalue comparisons are unaffected."
)
readme <- tr("LEEME", "README")
addWorksheet(wb, readme)
writeData(wb, readme, setNames(data.frame(notes), tr("nota", "note")))
setColWidths(wb, readme, cols = 1, widths = 110)
saveWorkbook(wb, file.path(folder, tr("cuadros_comparativos.xlsx", "comparison_tables_en.xlsx")), overwrite = TRUE)

scale <- read.csv(file.path(tables_dir, "escalabilidad.csv"))
geno <- read.csv(file.path(tables_dir, "genotipados.csv"))
accuracy <- read.csv(file.path(tables_dir, "precision_resumen.csv"))
theme_set(theme_bw(base_size = 12))
save_plot <- function(name, plot) {
  if (english) name <- paste0(name, "_en")
  ggsave(file.path(figures_dir, paste0(name, ".png")), plot, width = 9, height = 5.5, dpi = 180)
  ggsave(file.path(figures_dir, paste0(name, ".pdf")), plot, width = 9, height = 5.5)
}
scale$ajuste <- paste0("q=", scale$q, ", l=", scale$l)
save_plot("escalabilidad_tiempo",
  ggplot(scale, aes(n, tiempo_pca_h, colour = ajuste, group = ajuste)) +
    geom_line() + geom_point(size = 2) + facet_wrap(~operador) +
    scale_x_log10(breaks = c(10000, 50000, 100000)) + scale_y_log10() +
    labs(x = tr("Individuos totales (escala log)", "Total individuals (log scale)"),
         y = tr("Tiempo interno PCA (h, escala log)", "Internal PCA time (h, log scale)"),
         colour = tr("Ajuste", "Setting"),
         caption = tr("Una semilla. CHC de 100.000 sin resultados descargados.",
                      "One PCA seed. No downloaded CHC results at 100,000 individuals.")))
save_plot("genotipados_tiempo",
  ggplot(geno, aes(genotipados, tiempo_pca_min, colour = operador)) +
    geom_line() + geom_point(size = 2) + scale_x_continuous(breaks = c(500, 1000, 2000, 4000)) +
    labs(x = tr("Individuos genotipados", "Genotyped animals"),
         y = tr("Tiempo interno PCA (min)", "Internal PCA time (min)"),
         colour = tr("Operador", "Operator"),
         caption = tr("Pedigrí fijo de 20.000; q=4, l=20. Una semilla; resultado descriptivo.",
                      "Fixed pedigree: 20,000; q=4, l=20. One PCA seed; descriptive results.")))
accuracy$ajuste <- paste0("q=", accuracy$q, ", l=", accuracy$l)
save_plot("precision_eigenvalores",
  ggplot(accuracy, aes(n, error_eigenvalores_peor_pct, colour = ajuste, group = ajuste)) +
    geom_line() + geom_point(size = 2) + facet_wrap(~operador) + scale_y_log10() +
    scale_x_continuous(breaks = c(200, 500, 1000)) +
    labs(x = tr("Individuos totales", "Total individuals"),
         y = tr("Peor error relativo eigenvalores (%, escala log)", "Worst relative eigenvalue error (%, log scale)"),
         colour = tr("Ajuste", "Setting"),
         caption = tr("Máximo entre cinco eigenvalores y tres semillas de PCA por configuración.",
                      "Maximum across five eigenvalues and three PCA seeds per configuration.")))
cat("Excel and plots saved to", folder, "\n")
