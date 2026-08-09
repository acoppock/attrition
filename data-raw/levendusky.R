# Build data/levendusky.rda from the CGGK (2017) replication archive.
#
# Source: Coppock, Gerber, Green, and Kern (2016). Replication Data for:
# Combining double sampling and bounds to address non-ignorable missing
# outcomes in randomized experiments. Harvard Dataverse.
# doi:10.7910/DVN/AQB4MP
#
# Run with: source("data-raw/levendusky.R")

levendusky <- read.csv("data-raw/levendusky_mturk_clean.csv", stringsAsFactors = FALSE)

levendusky$Z_lev         <- factor(levendusky$Z_lev, levels = c("Placebo", "Moderate", "Polarized"))
levendusky$pid_3_recoded <- factor(levendusky$pid_3_recoded, levels = c("Dem", "Ind", "Rep"))
levendusky$Z_Levendusky  <- NULL   # duplicate of Z_lev in lower case

usethis::use_data(levendusky, overwrite = TRUE, compress = "xz")
