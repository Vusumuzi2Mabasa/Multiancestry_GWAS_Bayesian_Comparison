# Synthetic multiancestry pilot results

Status: **COMPLETE_WITH_MANC**

## Design and interpretation

Eight SNPs, 2,000 EUR-labelled individuals and nested AFR-labelled samples of 250, 500, 1,000 and 2,000. Three scenarios; 20 replicates per setting. These are stylised labels, not real ancestry data. See README.md for the complete model and limitations.

## Primary PIP >= 0.95 results

| Scenario | AFR N | Method | Exact recall | AFR-only recall | AFR RMSE | Mean FDP |
|---|---:|---|---:|---:|---:|---:|
| heterogeneous | 250 | Bayes_rho0_PIP0.95 | 0.275 | 0.100 | 0.078 | 0.025 |
| heterogeneous | 250 | Bayes_rho0.8_PIP0.95 | 0.300 | 0.000 | 0.086 | 0.050 |
| heterogeneous | 250 | Manc_COJO | 0.200 | 0.000 | 0.105 | 0.000 |
| heterogeneous | 500 | Bayes_rho0_PIP0.95 | 0.500 | 0.450 | 0.054 | 0.000 |
| heterogeneous | 500 | Bayes_rho0.8_PIP0.95 | 0.550 | 0.400 | 0.062 | 0.000 |
| heterogeneous | 500 | Manc_COJO | 0.200 | 0.000 | 0.104 | 0.000 |
| heterogeneous | 1000 | Bayes_rho0_PIP0.95 | 0.750 | 0.900 | 0.030 | 0.000 |
| heterogeneous | 1000 | Bayes_rho0.8_PIP0.95 | 0.800 | 0.900 | 0.035 | 0.000 |
| heterogeneous | 1000 | Manc_COJO | 0.225 | 0.050 | 0.100 | 0.000 |
| heterogeneous | 2000 | Bayes_rho0_PIP0.95 | 0.875 | 1.000 | 0.018 | 0.000 |
| heterogeneous | 2000 | Bayes_rho0.8_PIP0.95 | 0.900 | 1.000 | 0.020 | 0.000 |
| heterogeneous | 2000 | Manc_COJO | 0.425 | 0.550 | 0.072 | 0.000 |
| null | 250 | Bayes_rho0_PIP0.95 | NA | NA | 0.012 | 0.000 |
| null | 250 | Bayes_rho0.8_PIP0.95 | NA | NA | 0.007 | 0.000 |
| null | 250 | Manc_COJO | NA | NA | 0.000 | 0.000 |
| null | 500 | Bayes_rho0_PIP0.95 | NA | NA | 0.018 | 0.000 |
| null | 500 | Bayes_rho0.8_PIP0.95 | NA | NA | 0.011 | 0.000 |
| null | 500 | Manc_COJO | NA | NA | 0.000 | 0.000 |
| null | 1000 | Bayes_rho0_PIP0.95 | NA | NA | 0.013 | 0.000 |
| null | 1000 | Bayes_rho0.8_PIP0.95 | NA | NA | 0.011 | 0.000 |
| null | 1000 | Manc_COJO | NA | NA | 0.000 | 0.000 |
| null | 2000 | Bayes_rho0_PIP0.95 | NA | NA | 0.008 | 0.000 |
| null | 2000 | Bayes_rho0.8_PIP0.95 | NA | NA | 0.008 | 0.000 |
| null | 2000 | Manc_COJO | NA | NA | 0.000 | 0.000 |
| shared | 250 | Bayes_rho0_PIP0.95 | 0.600 | NA | 0.050 | 0.000 |
| shared | 250 | Bayes_rho0.8_PIP0.95 | 0.675 | NA | 0.039 | 0.000 |
| shared | 250 | Manc_COJO | 0.450 | NA | 0.077 | 0.000 |
| shared | 500 | Bayes_rho0_PIP0.95 | 0.675 | NA | 0.040 | 0.000 |
| shared | 500 | Bayes_rho0.8_PIP0.95 | 0.675 | NA | 0.037 | 0.000 |
| shared | 500 | Manc_COJO | 0.550 | NA | 0.065 | 0.000 |
| shared | 1000 | Bayes_rho0_PIP0.95 | 0.875 | NA | 0.028 | 0.000 |
| shared | 1000 | Bayes_rho0.8_PIP0.95 | 0.875 | NA | 0.026 | 0.000 |
| shared | 1000 | Manc_COJO | 0.675 | NA | 0.054 | 0.025 |
| shared | 2000 | Bayes_rho0_PIP0.95 | 1.000 | NA | 0.017 | 0.000 |
| shared | 2000 | Bayes_rho0.8_PIP0.95 | 1.000 | NA | 0.016 | 0.000 |
| shared | 2000 | Manc_COJO | 0.975 | NA | 0.021 | 0.000 |

## Limits

- Recall counts exact planted variants. Correlated proxies may recover a region without recovering the exact SNP; block recall is also tabulated.
- PIP and p-value thresholds are not matched for error rate. Do not infer method superiority directly from a recall difference.
- Zero false selections in 20 null replicates has a one-sided 95% upper bound of about 13.9% on the per-replicate probability; it does not establish genome-wide calibration.
- Shrinkage can improve total effect RMSE while biasing ancestry-specific effects. Inspect per-variant estimates, not only the aggregate score.
- The rho=0 model retains shared variant inclusion; it is not a fully independent-ancestry baseline.
- In-sample LD, fixed effect priors and a small quantitative locus favour clean calculations; external LD, binary/survival outcomes and real population structure remain untested.

Tables: summary.tsv (means and Monte Carlo SEs), paired_differences.tsv, replicate_metrics.tsv and variant_results.tsv. Figure: comparison.png. No real participant data were used.
