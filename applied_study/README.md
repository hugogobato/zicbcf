# MEPS 2023 dental-expenditure applied study

This directory reproduces the health-economics application of ZIC-BCF-Smear:
the causal effect of dental-insurance coverage on annual dental expenditure in
the 2023 Medical Expenditure Panel Survey (MEPS) Full-Year Consolidated file,
HC-251.

Everything here reproduces from a clean checkout. The only file that is not
committed is the raw AHRQ release, which is about 60 MB as CSV and is fetched
by the download script below.

## Data and conditions of use

The MEPS public-use files are released by the Agency for Healthcare Research
and Quality and are freely redistributable. Two restrictions attach to them and
survive redistribution: the data may be used for statistical reporting and
analysis only, and they may not be linked with individually identifiable
records from any source other than MEPS itself and the National Health
Interview Survey from which the MEPS sample is drawn.

Release page: <https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-251>

## How to reproduce

```sh
Rscript applied_study/download_meps_h251.R      # fetch and prepare h251.csv
Rscript applied_study/run_study.R               # ZIC-BCF-Smear fit, 4 chains
Rscript applied_study/run_cate_subgroups.R      # subgroup tables and figures
Rscript applied_study/run_gamma_hurdle.R        # Gamma-hurdle benchmark
Rscript applied_study/run_gamma_hurdle_subgroups.R
```

The two fitting scripts each run four chains sequentially and take roughly an
hour in total on a recent laptop. They discard each chain after summarizing it,
so peak memory stays near the size of a single fit.

## Files

| File | Role |
| --- | --- |
| `meps_common.R` | Cleaning rules, design matrix, subgroup definitions, survey-design helpers. Every other script sources it. |
| `download_meps_h251.R` | Fetches HC-251 from AHRQ and writes `h251.csv`. |
| `run_study.R` | ZIC-BCF-Smear fit; writes posterior summaries, convergence diagnostics and smearing diagnostics. |
| `run_cate_subgroups.R` | Subgroup effects, pairwise contrasts and figures, from the saved summaries. |
| `run_gamma_hurdle.R` | Gamma-hurdle benchmark fit. |
| `run_gamma_hurdle_subgroups.R` | Benchmark subgroup effects and contrasts. |

## Two estimands are reported for every quantity

MEPS is a stratified, clustered, unequal-probability sample, so a sample
average over its records is not a United States population quantity. Every
table therefore carries two columns.

**Analytic-sample average.** The unweighted average over the 18,763 analytic
records. It describes the sample, which oversamples several subpopulations.

**Survey-weighted population target.** The average weighted by `PERWT23F`, with
a design-based variance component obtained by Taylor linearization over the
strata (`VARSTR`) and primary sampling units (`VARPSU`) and added to the
posterior variance. Only this column supports a national reading. Income
tertile cutpoints are weighted quantiles for the same reason.

Clustering usually widens an interval and stratification usually narrows it, so
which effect dominates is not knowable without computing both, and the two
columns are reported side by side rather than one being presented as a
correction to the other.

## Note on the earlier version of this analysis

Three defects in the first version of these scripts are corrected here, and are
documented in the header of `meps_common.R`: nominal covariates were coerced to
numeric category codes instead of being expanded into indicators; race was taken
from `RACEV2X`, which carries no ethnicity dimension and placed 3,758 Hispanic
respondents inside the "White" category; and the survey design was ignored
entirely.
