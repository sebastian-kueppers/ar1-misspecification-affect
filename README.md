# Beyond Linear Dynamics: AR(1) Misspecification and Complexity in Daily Life Affect Data

Analysis code for Kueppers, S., Holtmann, J., Koslowski, K., Lange, J., Kalokerinos, E. K., & Koval, P. (2026). _Beyond Linear Dynamics: AR(1) Misspecification and Complexity in Daily Life Affect Data_.

- Preregistration: [https://doi.org/10.17605/OSF.IO/3CWXH](https://doi.org/10.17605/OSF.IO/3CWXH)
- Archived version of this repository: 

## Data availability

This repository contains **no raw or processed participant data**. We conducted secondary analyses of 18 experience sampling datasets from the EMOTE database (Kalokerinos et al., 2026; [https://emotedatabase.com](https://emotedatabase.com/)). 

To reproduce the analyses, request the data from EMOTE and quote **request ID `VRXLEPA17Y`**. This request ID specifies the exact datasets and variables we used. Our data were downloaded on 2026-06-11. The GEMS 2024 dataset was downloaded again on 2026-08-18. 

The Monte Carlo power and Type I error simulation (`power_sim/`) does not require EMOTE data. 

## Repository structure

```
.
├── _data preprocessing.R      # Step 1: raw EMOTE data -> processed data
├── _analysis.R                # Step 2: main analyses (RQ1-RQ3), main-text figures
├── supplementary analyses/    # Step 3: supplementary analyses (Supplementary Materials)
├── helpers/                   # functions sourced by the scripts (tests, refit bootstrap, signatures)
├── power_sim/                 # power and Type I error simulation (self-contained)
│   ├── power_sim.R
│   ├── Failenschmid-et-al-2025/   # third-party code (GPL-3), see README there
│   ├── data/                  # simulated data
│   ├── results/               # simulation results
│   └── figures/
├── data/
│   ├── raw/                   # <- place the EMOTE dataset CSVs here (not included)
│   │   └── summary/           # <- place the EMOTE request summary CSV here (not included)
│   └── processed/             # written by _data preprocessing.R
├── inputs/                    # variable coding sheets used by preprocessing
├── results/                   # written by _analysis.R and the supplementary scripts
└── figures/                   # figures as reported in the manuscript
```

The `data/` and `results/` folders are empty in this repository.

## How to reproduce

All scripts assume that the **working directory is the repository root**. 

1. **Get the data.** Request the EMOTE data with ID `VRXLEPA17Y`. Put the 18 dataset CSVs (`data_downloads_VRXLEPA17Y_*.csv`) in `data/raw/` and the request summary CSV (`data_reports_VRXLEPA17Y_*_request_summary_*.csv`) in `data/raw/summary/`.
2. **Preprocess.** Run `_data preprocessing.R`. It writes `data/processed/data_processed.Rdata` and the exclusion and reliability files to `results/`.
3. **Main analyses.** Run `_analysis.R`. RQ1 (AR(1) misspecification tests with the refit bootstrap null, B = 499) and RQ2 (complex-system signatures) are computationally intensive. Results are saved with a timestamp in `results/` (e.g., `rq1_results_<timestamp>.Rdata`).
4. **Supplementary analyses.** Run the scripts in `supplementary analyses/`. They load the RQ1/RQ2 results files from step 3. Because your files carry a new timestamp, update the file names in the `load()` calls at the top of each script.
5. **Simulation (optional, independent of steps 1–4).** `power_sim/power_sim.R` reproduces Tables S4–S6 and Figure S1 from the archived simulation files. Set the `RUN_*` flags at the top to `TRUE` to regenerate data or results.

## Software

Analyses were run in R 4.6.0 (R Core Team, 2026). Required packages:

ComplexUpset, data.table, dHSIC, dplyr, ecp, energy, furrr, future, future.apply, ggplot2, ggtext, glmmTMB, lubridate, mgcv, ordinal, parallelly, progressr, psych, purrr, readxl, rEDM, rlang, tidyr, tidyverse, tseries, tseriesEntropy, writexl

`rEDM` has changed its interface between versions, so use the version listed in `sessionInfo.txt`.

## License

The code in this repository is licensed under the GNU General Public License v3.0 (see `LICENSE`). `power_sim/Failenschmid-et-al-2025/` contains unmodified code by Failenschmid et al. (2025), which is also licensed under GPL-3.0.
