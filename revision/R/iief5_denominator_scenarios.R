# Post hoc hypothetical scenarios for the Fuschi IIEF-5 mean differences.
# These scenarios describe how the review's exploratory intervals respond to
# assumed follow-up counts. They do not estimate the unreported counts, do not
# replace any analysis, and do not model the missing-data mechanism.
options(stringsAsFactors = FALSE)
dir.create('results', showWarnings = FALSE)
d <- read.delim('data/inputs/functional_trajectory_source.tsv', na.strings = c('', 'NA'))
d <- d[d$panel == 'iief5', ]
N1 <- 97; N2 <- 98
z <- qnorm(.975)
scen <- data.frame(scenario = c('80% in both arms', '50% in both arms'),
                   f1 = c(.8, .5), f2 = c(.8, .5))
out <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
  md <- d$athermal_mean[i] - d$energy_mean[i]
  do.call(rbind, lapply(seq_len(nrow(scen)), function(j) {
    n1 <- round(N1 * scen$f1[j]); n2 <- round(N2 * scen$f2[j])
    se <- sqrt(d$athermal_sd[i]^2 / n1 + d$energy_sd[i]^2 / n2)
    data.frame(visit = d$time_label[i], scenario = scen$scenario[j], n_cold = n1, n_thermal = n2,
               md = md, low = md - z * se, high = md + z * se)
  }))
}))
write.table(out, 'results/iief5_denominator_scenarios.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

m <- read.delim('data/inputs/fuschi_source_reported_iief5_models.tsv')
write.table(m, 'results/iief5_source_reported_models.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
cat('Wrote hypothetical IIEF-5 denominator scenarios; no primary result changed.\n')
