options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
root_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input_path <- file.path(root_dir, "data", "inputs", "perioperative_source.tsv")
results_dir <- file.path(root_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

source_peri <- read.delim(input_path, check.names = FALSE)

md_effect <- function(athermal_mean, athermal_sd, athermal_n,
                      energy_mean, energy_sd, energy_n) {
  md <- athermal_mean - energy_mean
  se <- sqrt(athermal_sd^2 / athermal_n + energy_sd^2 / energy_n)
  c(effect = md, low = md - qnorm(0.975) * se,
    high = md + qnorm(0.975) * se, variance = se^2)
}

peri_values <- t(mapply(
  md_effect,
  source_peri$athermal_mean, source_peri$athermal_sd, source_peri$athermal_n,
  source_peri$energy_mean, source_peri$energy_sd, source_peri$energy_n
))
source_peri[, c("effect", "low", "high", "variance")] <- peri_values

write.table(
  source_peri,
  file.path(results_dir, "perioperative_effects.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

message("Wrote study-specific perioperative mean differences to results/perioperative_effects.tsv")
