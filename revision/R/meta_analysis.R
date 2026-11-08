options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
root_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input_dir <- file.path(root_dir, "data", "inputs")
results_dir <- file.path(root_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

read_input <- function(filename) {
  read.delim(file.path(input_dir, filename), check.names = FALSE)
}

meta_rr <- function(dat) {
  dat$rr <- (dat$athermal_events / dat$athermal_total) /
    (dat$energy_events / dat$energy_total)
  dat$yi <- log(dat$rr)
  dat$vi <- 1 / dat$athermal_events - 1 / dat$athermal_total +
    1 / dat$energy_events - 1 / dat$energy_total
  dat$se <- sqrt(dat$vi)
  dat$ci_low <- exp(dat$yi - qnorm(0.975) * dat$se)
  dat$ci_high <- exp(dat$yi + qnorm(0.975) * dat$se)

  w_fixed <- 1 / dat$vi
  mu_fixed <- sum(w_fixed * dat$yi) / sum(w_fixed)
  se_fixed <- sqrt(1 / sum(w_fixed))
  q <- sum(w_fixed * (dat$yi - mu_fixed)^2)
  df <- nrow(dat) - 1
  c_value <- sum(w_fixed) - sum(w_fixed^2) / sum(w_fixed)
  tau2_dl <- max(0, (q - df) / c_value)
  w_random <- 1 / (dat$vi + tau2_dl)
  mu_random <- sum(w_random * dat$yi) / sum(w_random)
  se_random <- sqrt(1 / sum(w_random))
  i2 <- if (q > 0) max(0, (q - df) / q) * 100 else 0

  reml_objective <- function(tau2) {
    w <- 1 / (dat$vi + tau2)
    mu <- sum(w * dat$yi) / sum(w)
    0.5 * (sum(log(dat$vi + tau2)) + log(sum(w)) +
      sum(w * (dat$yi - mu)^2))
  }
  tau2_reml <- optimize(reml_objective, interval = c(0, 5))$minimum
  w_reml <- 1 / (dat$vi + tau2_reml)
  mu_reml <- sum(w_reml * dat$yi) / sum(w_reml)
  se_reml_wald <- sqrt(1 / sum(w_reml))
  hk_scale <- sum(w_reml * (dat$yi - mu_reml)^2) / df
  se_reml_hk <- sqrt(hk_scale / sum(w_reml))
  hk_critical <- qt(0.975, df = df)

  dat$fixed_weight_pct <- 100 * w_fixed / sum(w_fixed)
  dat$dl_weight_pct <- 100 * w_random / sum(w_random)
  dat$reml_weight_pct <- 100 * w_reml / sum(w_reml)
  dat$random_weight_pct <- dat$reml_weight_pct

  prediction <- c(NA_real_, NA_real_)
  if (nrow(dat) >= 3) {
    critical <- qt(0.975, df = df)
    prediction <- exp(mu_reml + c(-1, 1) * critical *
      sqrt(tau2_reml + se_reml_hk^2))
  }

  summary <- data.frame(
    k = nrow(dat),
    fixed_rr = exp(mu_fixed),
    fixed_ci_low = exp(mu_fixed - qnorm(0.975) * se_fixed),
    fixed_ci_high = exp(mu_fixed + qnorm(0.975) * se_fixed),
    random_rr_dl = exp(mu_random),
    random_ci_low = exp(mu_random - qnorm(0.975) * se_random),
    random_ci_high = exp(mu_random + qnorm(0.975) * se_random),
    random_rr_reml_hk = exp(mu_reml),
    random_reml_hk_ci_low = exp(mu_reml - hk_critical * se_reml_hk),
    random_reml_hk_ci_high = exp(mu_reml + hk_critical * se_reml_hk),
    random_rr_reml_wald = exp(mu_reml),
    random_reml_wald_ci_low = exp(mu_reml - qnorm(0.975) * se_reml_wald),
    random_reml_wald_ci_high = exp(mu_reml + qnorm(0.975) * se_reml_wald),
    q = q,
    q_df = df,
    q_p = pchisq(q, df = df, lower.tail = FALSE),
    i2_pct = i2,
    tau2_dl = tau2_dl,
    tau2_reml = tau2_reml,
    hk_scale = hk_scale,
    prediction_low = prediction[1],
    prediction_high = prediction[2]
  )
  list(studies = dat, summary = summary)
}

leave_one_out <- function(dat) {
  rows <- lapply(seq_len(nrow(dat)), function(i) {
    fit <- meta_rr(dat[-i, , drop = FALSE])$summary
    data.frame(
      omitted_study_family = dat$study_family[i],
      omitted_study_label = dat$study_label[i],
      k = fit$k,
      random_rr_reml_hk = fit$random_rr_reml_hk,
      random_reml_hk_ci_low = fit$random_reml_hk_ci_low,
      random_reml_hk_ci_high = fit$random_reml_hk_ci_high,
      random_rr_reml_wald = fit$random_rr_reml_wald,
      random_reml_wald_ci_low = fit$random_reml_wald_ci_low,
      random_reml_wald_ci_high = fit$random_reml_wald_ci_high,
      i2_pct = fit$i2_pct,
      tau2_reml = fit$tau2_reml
    )
  })
  do.call(rbind, rows)
}

write_meta <- function(result, prefix) {
  write.table(
    result$studies,
    file.path(results_dir, paste0(prefix, "_study_effects.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  write.table(
    result$summary,
    file.path(results_dir, paste0(prefix, "_summary.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

draw_forest <- function(result, filename, title_text, subtitle_text) {
  dat <- result$studies
  pooled <- result$summary
  labels <- c(dat$study_label, "REML random-effects (HK)", "Common-effect sensitivity")
  effects <- c(dat$rr, pooled$random_rr_reml_hk, pooled$fixed_rr)
  lows <- c(dat$ci_low, pooled$random_reml_hk_ci_low, pooled$fixed_ci_low)
  highs <- c(dat$ci_high, pooled$random_reml_hk_ci_high, pooled$fixed_ci_high)
  y <- rev(seq_along(labels))
  x_min <- max(0.1, min(lows, na.rm = TRUE) * 0.75)
  x_max <- max(highs, na.rm = TRUE) * 1.25

  png(file.path(results_dir, filename), width = 1800, height = 1100, res = 180)
  par(mar = c(5, 12, 5, 4))
  plot(
    effects, y, log = "x", xlim = c(x_min, x_max), ylim = c(0.5, max(y) + 0.8),
    xlab = "Risk ratio for potency recovery (athermal / energy)", ylab = "",
    yaxt = "n", pch = c(rep(15, nrow(dat)), 18, 18), cex = c(rep(1.2, nrow(dat)), 1.6, 1.6),
    main = title_text, sub = subtitle_text
  )
  axis(2, at = y, labels = labels, las = 1, tick = FALSE)
  abline(v = 1, lty = 2, col = "#6B7280")
  segments(lows, y, highs, y, lwd = 2, col = c(rep("#334155", nrow(dat)), "#0F766E", "#B45309"))
  points(effects, y, pch = c(rep(15, nrow(dat)), 18, 18), cex = c(rep(1.2, nrow(dat)), 1.6, 1.6),
         col = c(rep("#334155", nrow(dat)), "#0F766E", "#B45309"))
  dev.off()
}

long_term <- meta_rr(read_input("meta_long_term_potency_input.tsv"))
early <- meta_rr(read_input("meta_3month_potency_input.tsv"))

write_meta(long_term, "meta_long_term_potency")
write_meta(early, "meta_3month_potency_exploratory")
write.table(
  leave_one_out(read_input("meta_long_term_potency_input.tsv")),
  file.path(results_dir, "meta_long_term_potency_leave_one_out.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

draw_forest(
  long_term,
  "forest_long_term_potency.png",
  "Long-term potency recovery after nerve-sparing RARP",
  sprintf("REML-HK RR %.2f (95%% CI %.2f–%.2f); I² %.1f%%; observational evidence",
          long_term$summary$random_rr_reml_hk,
          long_term$summary$random_reml_hk_ci_low,
          long_term$summary$random_reml_hk_ci_high,
          long_term$summary$i2_pct)
)

feasibility <- data.frame(
  outcome = c("Long-term potency recovery", "Three-month potency recovery"),
  independent_families = c(long_term$summary$k, early$summary$k),
  decision = c("Pool with prespecified REML-Hartung-Knapp primary and common-effect sensitivity",
               "Do not pool as a main result; use SWiM direction-of-effect synthesis"),
  rationale = c(
    sprintf("Three compatible observational families; I2=%.1f%% and the REML-Hartung-Knapp CI crosses 1, so certainty is very low. The prediction interval is descriptive only because k=3.", long_term$summary$i2_pct),
    sprintf("Only two small families with opposite effects; I2=%.1f%%, indicating an unstable pooled estimate.", early$summary$i2_pct)
  )
)
write.table(feasibility, file.path(results_dir, "meta_feasibility.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Long-term fixed RR:", long_term$summary$fixed_rr,
    "CI", long_term$summary$fixed_ci_low, long_term$summary$fixed_ci_high,
    "I2", long_term$summary$i2_pct, "\n")
cat("Long-term DL random RR:", long_term$summary$random_rr_dl,
    "CI", long_term$summary$random_ci_low, long_term$summary$random_ci_high, "\n")
cat("Long-term REML-HK random RR:", long_term$summary$random_rr_reml_hk,
    "CI", long_term$summary$random_reml_hk_ci_low,
    long_term$summary$random_reml_hk_ci_high, "\n")
cat("Long-term REML-Wald sensitivity RR:", long_term$summary$random_rr_reml_wald,
    "CI", long_term$summary$random_reml_wald_ci_low,
    long_term$summary$random_reml_wald_ci_high, "\n")
cat("Three-month I2:", early$summary$i2_pct, "- main pooling not recommended\n")
