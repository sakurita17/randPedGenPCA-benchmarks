# Semantic checks after the small smoke commands documented in followup_experiments.md.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
exact_dir <- args[1]
genotyped_dir <- args[2]
package_root <- args[3]
pkg <- new.env(parent = asNamespace("pedigreeTools"))
for (file in c("utils.R", "rhpca.R", "hutchH.R"))
  sys.source(file.path(package_root, "R", file), envir = pkg)

# Genotyping changes membership only; the pedigree and shared genotypes stay fixed.
master <- readRDS(file.path(genotyped_dir, "master.rds"))
previous_ids <- character()
for (count in c(8L, 16L, 32L, 40L)) {
  data <- readRDS(file.path(genotyped_dir, paste0("data_g", count, ".rds")))
  ids <- rownames(data$M)
  stopifnot(identical(data$ped, master$ped), nrow(data$M) == count,
            all(previous_ids %in% ids),
            identical(data$M, master$M[ids, , drop = FALSE]))
  previous_ids <- ids
}

# The dense reference must represent exactly the implicit operator used by rhpca.
data <- readRDS(file.path(exact_dir, "data_80.rds"))
Li <- pedigreeTools::getLInv(data$ped)
for (center in c(FALSE, TRUE)) {
  ref <- readRDS(file.path(exact_dir, paste0("n80_c", center, "_exact.rds")))
  op <- pkg$makeHxOperator(pkg$inputHx(Li, data$M), center, "0.5", NULL,
                           0.01, TRUE, 1000L)
  set.seed(345)
  V <- matrix(rnorm(80 * 7), 80, 7)
  stopifnot(norm(op(V) - ref$K %*% V, "F") / norm(ref$K %*% V, "F") < 1e-10)
  # With a full-size sketch, the leading eigenvalues should match the reference.
  set.seed(123)
  fit <- pkg$rhpca.pedigree(data$ped, data$M, rho = 3, q = 0, l = 80,
                           centredH = center, centred = "0.5", epsilon = 0.01)
  stopifnot(max(abs(fit$d - ref$d) / abs(ref$d)) < 1e-8)
}
for (folder in c(exact_dir, genotyped_dir)) {
  rows <- read.csv(file.path(folder, "results.csv"))
  stopifnot(all(rows$status == "ok"), all(rows$process_wall_sec >= 0),
            all(rows$peak_process_rss_mib > 0))
}
cat("Passed: nested datasets, fixed pedigree/genotypes, exact/implicit operator, full-sketch eigenvalues\n")
