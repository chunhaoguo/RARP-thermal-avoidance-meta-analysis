# Thermal avoidance versus energy use in nerve-sparing RARP

This repository contains the redistributable derived data and analysis code supporting the systematic review and meta-analysis:

> Thermal avoidance versus energy use at the prostatic pedicle and neurovascular bundle during nerve-sparing robot-assisted radical prostatectomy: a systematic review and meta-analysis.

The review was prospectively registered in PROSPERO (CRD420261488794).

## Contents

- `data/inputs/`: study-level inputs used for the quantitative syntheses.
- `data/derived/`: derived tables underlying the main and supplementary results and figures.
- `R/meta_analysis.R`: primary and sensitivity meta-analyses for potency recovery, including Hartung-Knapp and common-effect estimates and leave-one-family-out analyses.
- `R/perioperative_effects.R`: study-specific mean differences for perioperative outcomes.
- `results/`: output directory created when the scripts are run.

## Reproduce the analyses

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
