build_no_metafounder <- function(sim, config = NULL) {
  list(
    name = "hybrid_no_metafounder",
    nmf = 0L,
    n_pedigree = length(sim$ped@label),
    ped = sim$ped,
    gamma = NULL,
    M = sim$M,
    notes = "Full A/B/AB pedigree without metafounders."
  )
}

build_single_metafounder <- function(sim, config) {
  mf <- "-1"
  ped_df <- sim$ped_df
  sire <- ped_df$fid
  dam <- ped_df$mid
  sire[is.na(sire)] <- mf
  dam[is.na(dam)] <- mf

  ped_mf <- pedigreeTools::pedigree(
    sire = c(NA, sire),
    dam = c(NA, dam),
    label = c(mf, ped_df$id)
  )

  gamma_t <- estimate_gamma(ped_mf, sim$M, nmf = 1L, config)

  list(
    name = "hybrid_single_metafounder",
    nmf = 1L,
    n_pedigree = length(ped_mf@label),
    ped = ped_mf,
    gamma = gamma_t$value$gamma,
    gamma_timing = gamma_t,
    M = sim$M,
    notes = "Full A/B/AB pedigree; all unknown parents assigned to one metafounder."
  )
}

build_two_metafounders <- function(sim, config) {
  ped_df <- sim$ped_df
  built <- build_two_mf_pedigree(ped_df)
  gamma_t <- estimate_gamma(built$ped, sim$M, nmf = 2L, config)

  list(
    name = "hybrid_two_metafounders",
    nmf = 2L,
    n_pedigree = length(built$ped@label),
    ped = built$ped,
    gamma = gamma_t$value$gamma,
    gamma_timing = gamma_t,
    M = sim$M,
    notes = "Full A/B/AB pedigree; base A assigned to MF1 and base B assigned to MF2."
  )
}

build_two_mf_pedigree <- function(ped_df) {
  mf1 <- "-1"
  mf2 <- "-2"
  sire <- ped_df$fid
  dam <- ped_df$mid
  pop <- as.character(ped_df$population)

  missing_sire <- is.na(sire)
  missing_dam <- is.na(dam)

  sire[missing_sire & pop == "A"] <- mf1
  dam[missing_dam & pop == "A"] <- mf1
  sire[missing_sire & pop == "B"] <- mf2
  dam[missing_dam & pop == "B"] <- mf2

  sire[missing_sire & !(pop %in% c("A", "B"))] <- mf1
  dam[missing_dam & !(pop %in% c("A", "B"))]  <- mf2

  ped_mf <- pedigreeTools::pedigree(
    sire = c(NA, NA, sire),
    dam = c(NA, NA, dam),
    label = c(mf1, mf2, ped_df$id)
  )

  list(ped = ped_mf, mf1 = mf1, mf2 = mf2)
}

estimate_gamma <- function(ped, M, nmf, config) {
  timed(pedigreeTools::getGammaL(
    ped = ped,
    M = M,
    nmf = nmf,
    epsilon = config$gammaL$epsilon
  ))
}

scenario_builders <- list(
  hybrid_no_metafounder = build_no_metafounder,
  hybrid_single_metafounder = build_single_metafounder,
  hybrid_two_metafounders = build_two_metafounders
)
