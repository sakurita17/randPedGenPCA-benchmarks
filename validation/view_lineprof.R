# Source this in a local R session, then call view_lineprof("path/to/g500").
view_lineprof <- function(folder) {
  if (!requireNamespace("lineprof", quietly = TRUE)) stop("Install lineprof first.")
  folder <- normalizePath(folder, mustWork = TRUE)
  profile <- readRDS(file.path(folder, "lineprof.rds"))
  # The profile was collected on Eddie; point source references at downloaded copies.
  for (i in seq_len(nrow(profile))) {
    refs <- profile$ref[[i]]
    for (j in which(!is.na(refs$path))) {
      candidate <- file.path(folder, "source", basename(refs$path[j]))
      if (file.exists(candidate)) refs$path[j] <- normalizePath(candidate, mustWork = TRUE)
    }
    profile$ref[[i]] <- refs
  }
  lineprof::shine(profile)
}
