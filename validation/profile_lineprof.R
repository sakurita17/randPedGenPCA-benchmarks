# Use Hadley Wickham's legacy lineprof method, as in Advanced R first edition.
# Rscript profile_lineprof.R DATASET PACKAGE_ROOT OUTPUT [CENTER Q L TRACE SEED]
args <- commandArgs(TRUE)
if (length(args) < 3L) stop("Specify DATASET PACKAGE_ROOT OUTPUT [CENTER Q L TRACE SEED]")
if (!requireNamespace("lineprof", quietly = TRUE))
  stop("Install hadley/lineprof in the project R library before submitting this job.")
dataset <- normalizePath(args[1], mustWork = TRUE)
package_root <- normalizePath(args[2], mustWork = TRUE)
output <- args[3]
center <- if (length(args) >= 4L) as.logical(args[4]) else FALSE
q <- if (length(args) >= 5L) as.integer(args[5]) else 4L
l <- if (length(args) >= 6L) as.integer(args[6]) else 20L
trace_vectors <- if (length(args) >= 7L) as.integer(args[7]) else 50L
seed <- if (length(args) >= 8L) as.integer(args[8]) else 123L
stopifnot(!is.na(center), q >= 0L, l >= 5L,
          trace_vectors == 0L || trace_vectors >= 2L, seed > 0L)
if (dir.exists(output)) stop("Output folder must be new: ", output)
dir.create(file.path(output, "source"), recursive = TRUE)
output <- normalizePath(output, mustWork = TRUE)
files <- file.path(package_root, "R", c("utils.R", "rhpca.R", "hutchH.R"))
stopifnot(all(file.copy(files, file.path(output, "source"))))
files <- file.path(output, "source", basename(files))
md5 <- tools::md5sum(files)
write.csv(data.frame(file = names(md5), md5 = unname(md5)),
          file.path(output, "sources_md5.csv"), row.names = FALSE)
pkg <- new.env(parent = asNamespace("pedigreeTools"))
for (file in files) sys.source(file, envir = pkg, keep.source = TRUE)
writeLines(unlist(lapply(c("getLInv", "getHx", "implicitMasuda", "oraculumLi"), function(name) {
  c(paste0("# ", name), deparse(get(name, envir = asNamespace("pedigreeTools"))), "")
})), file.path(output, "pedigreeTools_functions.R"))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
data <- readRDS(dataset)
cat("lineprof:", data$n, "pedigree individuals;", nrow(data$M), "genotyped\n")

warnings <- character()
failure <- NULL
fit <- NULL
completed <- FALSE
run_pca <- function() {
  # lineprof catches errors internally, so explicitly record whether PCA finished.
  tryCatch({
    fit <<- withCallingHandlers(pkg$rhpca.pedigree(
      data$ped, data$M, rho = 5L, q = q, l = l, centredH = center,
      centred = "0.5", epsilon = 0.01,
      numVectorsTrace = if (trace_vectors == 0L) NULL else trace_vectors),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    completed <<- TRUE
  }, error = function(e) {failure <<- conditionMessage(e)})
  invisible(NULL)
}
gc()
set.seed(seed)
timing <- system.time(
  profile <- lineprof::lineprof(run_pca(), interval = 0.02, torture = FALSE)
)
saveRDS(profile, file.path(output, "lineprof.rds"))
# format.lineprof flattens the reference column while retaining time/memory columns.
write.csv(format(profile, digits = 6), file.path(output, "lineprof.csv"), row.names = FALSE)
writeLines(capture.output(print(profile, depth = 5)), file.path(output, "summary.txt"))
write.csv(data.frame(n = data$n, n_genotyped = nrow(data$M), markers = ncol(data$M),
  centredH = center, q = q, l = l, rho = 5L, trace_vectors = trace_vectors, seed = seed,
  profiled_elapsed_sec = unname(timing["elapsed"]),
  status = if (completed) "ok" else "error", warnings = paste(warnings, collapse = "; "),
  error = if (is.null(failure)) "" else failure),
  file.path(output, "timing.csv"), row.names = FALSE)
if (!completed) stop("PCA did not finish: ", failure)
saveRDS(fit, file.path(output, "fit.rds"), compress = FALSE)
cat("Finished. lineprof results:", output, "\n")
