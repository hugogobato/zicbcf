# Alternative endpoint analysis: distinct severe-AE calendar days

## Purpose

This analysis replaces the earlier outcome, the sum of durations across individual grade 3+ adverse-event (AE) records, with the number of distinct on-study calendar days on which a patient experienced at least one grade 3+ AE. The new outcome avoids counting a date twice when two qualifying events overlap. It is a sensitivity analysis of outcome definition, not a correction to an error in the prior analysis.

## Data linkage and endpoint construction

The supplied `ae.sas7bdat` file contains 9,934 AE records for 514 patients. Its original `SUBJID` field links exactly to `../headneck_265_severe_ae_patient_level.csv`, which contains the 520 randomized patients, their treatment arm, and the prespecified baseline covariates. Recomputing the prior summed-duration outcome from the AE file reproduced the existing patient-level value for every one of the 520 patients.

The alternative endpoint uses on-study records (`AEPRIOR != "Y"`, or missing) with `AESEVCD >= 3`. For each qualifying record, the imputed study-day interval `[AESTDYI, AEENDYI]` is inclusive, as specified by the supplied data dictionary: `AEDURI = (AEENDTI - AESTDTI) + 1`. Within each patient, overlapping or adjacent intervals are merged and counted once. Of 1,671 qualifying records, 1,547 (92.6%) have an imputed start and end day. The 124 records with an unresolved boundary are not extended beyond the observed data because the supplied files do not contain patient-specific follow-up end dates. They also contribute zero days to the existing duration endpoint. Complete date intervals reproduce that endpoint exactly.

## Descriptive comparison

| Endpoint | Total days | Chemotherapy mean | Panitumumab plus chemotherapy mean | Unadjusted arm difference | Zero outcomes |
|---|---:|---:|---:|---:|---:|
| Summed AE-record days | 17,988 | 26.75 | 42.43 | +15.68 | 123 |
| Distinct severe-AE calendar days | 13,107 | 19.18 | 31.23 | +12.06 | 123 |

The interval union removes 4,881 record-days (27.1% of the original total) across 217 patients. The unadjusted arm difference decreases by 3.62 days (23.1%), but remains positive. The zero-outcome count is unchanged because every positive recorded duration is associated with a complete date interval in this extract.

## Adjusted causal estimates

Both analyses use the same 520 patients, treatment coding, baseline adjustment set, and posterior settings. The estimand is the average effect of panitumumab plus chemotherapy relative to chemotherapy.

| Model | Outcome | Record-duration ATE, 95% CrI | Distinct-calendar-day ATE, 95% CrI | Posterior probability of a positive distinct-day ATE |
|---|---|---:|---:|---:|
| ZIC-BCF-Smear | Total severe-AE burden | +13.99 [4.10, 24.00] | +9.95 [3.33, 16.82] | 0.999 |
| Gamma hurdle | Total severe-AE burden | +16.06 [6.25, 27.37] | +12.31 [6.14, 19.25] | 0.999 |

Both models retain a credible positive effect after overlapping events are counted once. The expected reduction in magnitude is consistent with the endpoint change rather than a reversal of the safety conclusion. The Gamma-hurdle result is numerically closer to the unadjusted distinct-day contrast (+12.31 versus +12.06 days) than the ZIC-BCF estimate (+9.95 days), but raw-arm proximity is not causal-model validation.

The hurdle-margin conclusions are essentially unchanged because the zero indicator does not change under the interval union. For ZIC-BCF, the distinct-day probability effect is +2.79 percentage points (95% CrI -3.81 to +10.25; posterior probability positive 0.785); it was +2.68 points (-3.95 to +9.84) for record-days. For Gamma hurdle, it is +1.76 points (-5.56 to +9.13; posterior probability positive 0.675), exactly as in the record-duration analysis. Neither model supports a credible treatment effect on the probability of any severe-AE day.

An independent Gamma-hurdle run with seed 2026 gives a distinct-day ATE of +12.36 days (95% CrI 5.30 to 19.74; posterior probability positive 1.000), supporting stability of the overall Gamma result.

## Conditional effects

ZIC-BCF does not produce a credible distinct-day CATE contrast for any prespecified subgroup. Its ECOG 0 minus ECOG 1 contrast is +6.79 days (95% CrI -0.84 to 20.87). Gamma hurdle retains a larger ECOG 0 minus ECOG 1 contrast of +23.54 days (7.28 to 41.67), while the corresponding ZIC-BCF interval includes zero. Therefore this apparent heterogeneity is model-dependent and should remain exploratory, not a substantive subgroup conclusion.

## Interpretation and literature convergence

The new endpoint improves the clinical interpretation of a day of toxicity: concurrent severe events count as one day experienced with severe toxicity. It nevertheless remains an observed-complete-interval endpoint, so ongoing or administratively censored events with unresolved end dates are not assigned unobserved days. The primary conclusion is robust: adding panitumumab is associated with more severe-toxicity burden, whether burden is represented by summed event durations or by distinct calendar days.

The SPECTRUM literature reports grade 3/4 AE incidence and specific toxicity categories rather than either duration endpoint. It therefore supports only the direction of the overall finding, not the numerical ATE. Gamma hurdle does not converge more strongly with the literature than ZIC-BCF on that basis: both estimates are directionally consistent with the reported excess severe toxicity. The distinct-day outcome is even more clearly a new estimand, because it explicitly handles co-occurring AEs. No supplied publication externally corroborates the Gamma-specific ECOG signal.

## Potential contribution and remaining literature review

This analysis does not newly establish that panitumumab causes severe toxicity. That clinical direction is already reported in SPECTRUM and related panitumumab studies. The potential contribution is instead a new secondary analysis of the existing randomized-trial data: 1. a patient-centred endpoint defined as distinct calendar days with at least one grade 3+ AE, rather than AE counts or summed record durations; 2. a randomized-treatment ATE estimated on that day scale; 3. a direct sensitivity analysis showing that the positive overall conclusion persists after concurrent events are not double-counted; and 4. flexible, uncertainty-aware subgroup-effect estimation for a semicontinuous outcome.

The appropriate manuscript wording at this stage is therefore cautious: "To our knowledge, this is a new secondary causal analysis of the SPECTRUM data using an overlap-adjusted severe-AE-day endpoint." It should not be described as a new clinical discovery, and it should not yet claim that no previous publication has used an equivalent endpoint or analysis. A targeted additional literature review is still required to verify that originality claim, including searches for SPECTRUM secondary safety analyses, clinical-study reports, regulatory materials, and studies using AE-duration, time-with-toxicity, overlap-adjusted, or burden-based endpoints. Until that review is complete, the work is best positioned as a promising applied-study contribution within a methodological manuscript.

## Reproducibility

The endpoint is constructed by `../build_nonoverlap_calendar_day_endpoint.R`. The matched ZIC-BCF and Gamma-hurdle fits are run by `../run_nonoverlap_calendar_day_analysis.R`; the independent Gamma diagnostic is run by `../run_nonoverlap_gamma_seed_sensitivity.R`. Machine-readable data, posterior draws, subgroup results, and figures are in this directory.
