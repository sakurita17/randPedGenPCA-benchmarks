# Generate exact TOTAL pedigree sizes, with every individual genotyped.
# Run from the repository root: Rscript R/generate_pedigrees_and_markers.R
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)][1L])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root, "R", "simulate_abc_population.R"))

sizes <- c(2000L, 5000L, 10000L)
markers <- 1000L
generations <- 10L
base_seed <- 12345L
output_dir <- file.path(root, "output", paste0(
  format(Sys.time(), "%Y%m%d_%H%M%S"), "_ABC_populations"))
if (!dir.create(output_dir, recursive = TRUE)) stop("Cannot create output directory")

manifest <- data.frame()
for (n in sizes) {
  prefix <- paste0("population_", n)
  dataset <- file.path(output_dir, paste0(prefix, ".rds"))
  seed <- base_seed + n
  message("Generating ", n, " total individuals with ", markers, " markers")
  data <- simulate_abc_population(n, markers, seed)
  ped <- data$ped_df
  ids <- ped$id
  M <- data$M[ids, , drop = FALSE]
  colnames(M) <- paste0("SNP_", seq_len(ncol(M)))
  stopifnot(nrow(ped) == n, !anyDuplicated(ids),
            identical(rownames(M), ids), nrow(M) == n, ncol(M) == markers,
            all(M %in% 0:2))
  offspring <- which(!is.na(ped$fid) & !is.na(ped$mid))
  fathers <- match(ped$fid[offspring], ids)
  mothers <- match(ped$mid[offspring], ids)
  stopifnot(!anyNA(fathers), !anyNA(mothers),
            all(fathers < offspring), all(mothers < offspring))
  # A/B breed within their group; C is A x B in generation five, then C x C.
  for (group in c("A", "B", "C")) {
    rows <- which(ped$population == group & ped$generation > 0L)
    sire_group <- ped$population[match(ped$fid[rows], ids)]
    dam_group <- ped$population[match(ped$mid[rows], ids)]
    initial_cross <- group == "C" & ped$generation[rows] == data$hybrid_generation
    stopifnot(all(sire_group[!initial_cross] == group),
              all(dam_group[!initial_cross] == group))
    if (any(initial_cross)) {
      stopifnot(all(sire_group[initial_cross] == "B"),
                all(dam_group[initial_cross] == "A"))
    }
  }
  # Check that each observed dosage is possible given its two parents.
  for (start in seq.int(1L, length(offspring), by = 500L)) {
    block <- start:min(start + 499L, length(offspring))
    paternal <- M[fathers[block], , drop = FALSE]
    maternal <- M[mothers[block], , drop = FALSE]
    child <- M[offspring[block], , drop = FALSE]
    stopifnot(all(child >= (paternal == 2L) + (maternal == 2L)),
              all(child <= (paternal > 0L) + (maternal > 0L)))
  }

  data$ped_df <- ped
  data$M <- M
  saveRDS(data, dataset)
  pedigree_file <- paste0(prefix, "_pedigree.csv")
  marker_file <- paste0(prefix, "_markers.csv")
  write.csv(ped, file.path(output_dir, pedigree_file), row.names = FALSE, na = "0")
  write.csv(data.frame(id = ids, M, check.names = FALSE),
            file.path(output_dir, marker_file), row.names = FALSE)
  manifest <- rbind(manifest, data.frame(
    n_pedigree = n, n_genotyped = nrow(M), n_markers = ncol(M),
    n_A = sum(ped$population == "A"), n_B = sum(ped$population == "B"),
    n_C = sum(ped$population == "C"),
    generations = generations, seed = seed, pedigree_file = pedigree_file,
    marker_file = marker_file, rds_file = basename(dataset)
  ))
  write.csv(manifest, file.path(output_dir, "population_manifest.csv"), row.names = FALSE)
  message("Validated and saved ", n, " individuals")
}

writeLines(c(
  "Tres poblaciones simuladas independientes: 2000, 5000 y 10000 individuos totales.",
  "Cada pedigree contiene 10 generaciones, incluidos los fundadores (generacion 0).",
  "Todos los individuos tienen 1000 marcadores, codificados como dosis alelicas 0/1/2.",
  "Modelo AlphaSimR: runMacs, GENERIC, split=20, dos cromosomas de 500 sitios cada uno.",
  "A y B son grupos fundadores separados que se reproducen dentro de su grupo.",
  "C (equivalente al grupo AB) nace de cruces madre A x padre B en la generacion 5;",
  "despues C se reproduce mediante C x C. Todos los padres son de la generacion anterior.",
  "Hay A/B en generaciones 0-4 y A/B/C en generaciones 5-9, con igual tamano por grupo/generacion.",
  "El total acumulado por pedigree se reparte en 40% A, 40% B y 20% C.",
  "Se simulan herencia mendeliana, ligamiento y recombinacion; no se aplica seleccion.",
  "CSV pedigree: id, population (A/B/C), fid (padre), mid (madre), generation;",
  "0 en las columnas fid/mid significa progenitor desconocido.",
  "CSV marcadores: id, SNP_1, ..., SNP_1000; el orden coincide con el pedigree.",
  "Los identificadores son locales a cada poblacion.",
  "RDS: lista con ped (pedigreeTools), ped_df (data.frame), M (matriz con IDs en rownames),",
  "y metadatos de simulacion. Para cargar: datos <- readRDS('population_2000.rds').",
  "Verificaciones: tamanos exactos, IDs alineados, padres anteriores a hijos,",
  "dosis 0/1/2, compatibilidad mendeliana y origen parental de cada grupo A/B/C."
), file.path(output_dir, "README.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
message("Done. Populations saved in: ", output_dir)
