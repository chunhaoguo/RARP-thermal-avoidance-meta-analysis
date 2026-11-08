# Thermal avoidance versus energy use in nerve-sparing RARP

This repository contains the redistributable derived data and analysis code supporting the systematic review and meta-analysis:

> Thermal avoidance versus energy use at the prostatic pedicle and neurovascular bundle during nerve-sparing robot-assisted radical prostatectomy: a systematic review and meta-analysis.

The review was prospectively registered in PROSPERO (CRD420261488794).

## Revised analyses

The `revision/` directory contains the analysis scripts, aggregate inputs and result tables for the revised manuscript. The original-submission files remain at the repository root for version history.

Run the revised calculations from that directory using base R 4.6.0 or a compatible installation:

```bash
cd revision
Rscript R/meta_analysis.R
Rscript R/perioperative_effects.R
Rscript R/functional_reproduction.R
Rscript R/iief5_denominator_scenarios.R
Rscript R/later_continence.R
```

See [`revision/README.txt`](revision/README.txt) for the result-to-manuscript map and analysis roles. [`revision/MANIFEST_SHA256.tsv`](revision/MANIFEST_SHA256.tsv) records the archived file hashes. The primary analysis is REML with Hartung-Knapp intervals. The IIEF-5 reconstructions and hypothetical completion scenarios are explicitly exploratory, and no meta-regression was fitted because the available evidence did not support it.

## Original submission

### Contents

- `data/inputs/`: study-level inputs used for the quantitative syntheses.
- `data/derived/`: derived tables underlying the main and supplementary results and figures.
- `R/meta_analysis.R`: primary and sensitivity meta-analyses for potency recovery, including Hartung-Knapp and common-effect estimates and leave-one-family-out analyses.
- `R/perioperative_effects.R`: study-specific mean differences for perioperative outcomes.
- `results/`: output directory created when the scripts are run.

### Reproduce the original analyses

The scripts use base R only. From the repository root, run:

```bash
Rscript R/meta_analysis.R
Rscript R/perioperative_effects.R
```

The scripts write tab-separated result tables and the long-term potency forest plot to `results/`.

## Data scope

Only redistributable study-level inputs and author-derived data are included. Copyrighted full-text articles and licensed bibliographic database exports are not redistributed. Source publications are cited in the manuscript.

## Licences

- Analysis code: MIT License (`LICENSE-CODE`).
- Derived data: Creative Commons Attribution 4.0 International (`LICENSE-DATA`).

## Contact

Questions about the repository may be directed to the corresponding authors listed in the manuscript.
