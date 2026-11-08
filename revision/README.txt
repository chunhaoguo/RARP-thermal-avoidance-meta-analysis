RARP reproducibility files

Run from this folder with base R 4.6.0 or a compatible installation:
Rscript R/meta_analysis.R
Rscript R/perioperative_effects.R
Rscript R/functional_reproduction.R
Rscript R/iief5_denominator_scenarios.R
Rscript R/later_continence.R

Output mapping
Main Figure 3 and Table 3: results/functional_effects_verified.tsv.
Supplementary Tables S7A-S7B: results/meta_long_term_potency_summary.tsv and results/meta_long_term_potency_leave_one_out.tsv.
Main Table 2 and Supplementary Table S11: results/perioperative_effects.tsv.
Exploratory observational blood loss: results/blood_loss_exploratory.tsv.
Supplementary Table S9 source-reported coefficients: results/iief5_source_reported_models.tsv.
Supplementary Table S10 illustrative scenarios: results/iief5_denominator_scenarios.tsv.
Basourakos 0-pad continence before 12 months: results/later_continence_basourakos.tsv.

The primary analysis is REML with Hartung-Knapp intervals. Model and leave-one-study-out analyses are sensitivities. Three-month pooled outputs are protocol-transparency calculations, not a reported primary pooled result.
IIEF-5 intervals assume 97 and 98 observations per visit unless a hypothetical scenario is stated. Scenarios hold published means and SDs fixed and do not model differential attrition. Adjusted coefficients are source-reported.
No meta-regression was fitted; Supplementary Table S8 records feasibility.
