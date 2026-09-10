library(tidyverse)

# Build data/levendusky_replication.rda from the CGGK (2017) replication archive.
#
# Run with: source("data-raw/levendusky_replication.R")

# Fetch from the Harvard Dataverse ----

# The replication archive for Coppock, Gerber, Green, and Kern (2017) is
# doi:10.7910/DVN/AQB4MP. Files there are addressed by a numeric id that is
# stable across dataset versions, so this script can name the exact file it
# read: 2887314 is levendusky_mturk_clean.csv, the cleaned subject-level file
# every published result runs on.
#
# `?format=original` matters: without it Dataverse serves its own ingested TSV
# rather than the CSV as deposited.

cache_dir <- "data-raw/cache"

dataverse_file <- function(file_id) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(cache_dir, file_id)
  if (!file.exists(dest)) {
    download.file(
      paste0("https://dataverse.harvard.edu/api/access/datafile/", file_id, "?format=original"),
      dest, mode = "wb", quiet = TRUE
    )
  }
  dest
}

# levendusky_replication ----

# The archive file holds all three conditions. Every result in the paper is the
# polarized-versus-moderate contrast, so the placebo condition goes here rather
# than in each analysis, and the two contrasts defined only against the placebo
# (Z2 and Z3 in the archive) go with it. Z_Levendusky is a lower-case duplicate
# of Z_lev and is dropped as well.

levendusky_replication <-
  read_csv(dataverse_file(2887314), show_col_types = FALSE) |>
  filter(Z_lev %in% c("Moderate", "Polarized")) |>
  transmute(
    X_party_id = factor(pid_3_recoded, levels = c("Dem", "Ind", "Rep")),
    Z_condition = factor(Z_lev, levels = c("Moderate", "Polarized")),
    Z = Z1,
    R1 = R1,
    Attempt = Attempt,
    R2 = R2,
    Y_polarization_w1 = L_dif,
    Y_polarization_w2 = L_dif_w2,
    Y_extremity_w1 = L_ex,
    Y_extremity_w2 = L_ex_w2
  )

usethis::use_data(levendusky_replication, overwrite = TRUE, compress = "xz")
