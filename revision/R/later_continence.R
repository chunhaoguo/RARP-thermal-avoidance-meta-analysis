# Basourakos 2022 Table 4, footnote c: any report of 0-pad continence at a follow-up visit before 12 months.
# This is not a fixed 12-month assessment or the adjusted >12-month EPIC-CP score model.
# Unadjusted, non-randomised; reported descriptively and not pooled with the randomised trial.
options(stringsAsFactors = FALSE)
dir.create('results', showWarnings = FALSE)
a <- 61; n1 <- 159; b <- 67; n2 <- 133
rr <- (a / n1) / (b / n2)
se <- sqrt(1 / a - 1 / n1 + 1 / b - 1 / n2)
out <- data.frame(study = 'Basourakos 2022', outcome = 'Any report of 0 pads before 12 months', clips = paste0(a, '/', n1), bipolar = paste0(b, '/', n2),
                  rr_clips_vs_bipolar = rr, low = exp(log(rr) - qnorm(.975) * se), high = exp(log(rr) + qnorm(.975) * se))
write.table(out, 'results/later_continence_basourakos.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
cat(sprintf('Clips vs bipolar RR %.2f (%.2f-%.2f)\n', rr, out$low, out$high))
