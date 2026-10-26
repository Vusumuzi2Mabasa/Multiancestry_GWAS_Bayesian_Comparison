# Multiancestry simulation pilot

**Completed synthetic pilot: 240 paired datasets analysed with the official
Manc-COJO v1.1.0 executable and a custom Bayesian hierarchical selection model.**

- [Detailed Word report with methods before results](report/Multiancestry_GWAS_Methodology_and_Results.docx)
- [Compact results](results/RESULTS.md)
- [R report source](scripts/build_word_report.R)
- [Signal recovery figure](results/figures/01_signal_recovery.png)
- [Effect estimation figure](results/figures/02_effect_error.png)

Strong effect borrowing helped in the shared-effect simulation but increased
error relative to the uncorrelated-effect prior when effects differed. These
are pilot findings, not proof that a Bayesian method is superior: decision
thresholds are not error-rate-matched, only 20 replicates were used per setting,
and the cohorts and LD are deliberately simplified.

This project compares the actual **Manc-COJO executable** with a small,
custom Bayesian hierarchical variable-selection model. If that executable
is unavailable, the program clearly records that Manc-COJO was **not run**;
it never labels a home-made stepwise algorithm as Manc-COJO.

## Scientific question

When one ancestry cohort is smaller, how do cross-population borrowing and
different LD patterns affect recovery of known signals and effect estimation?
The simulation is deliberately small enough for a laptop. It is not a
genome-wide benchmark or a new validated fine-mapping method.

## Design fixed before examining results

- Quantitative phenotype; eight biallelic SNPs in two independent four-SNP blocks.
- Larger cohort labelled EUR: 2,000 people. Smaller cohort labelled AFR:
  250, 500, 1,000 or 2,000 people. Labels describe artificial scenarios only.
- Two Markov-generated haplotypes per person; allele frequency 0.3.
  Within-block haplotype persistence is 0.8 in EUR and 0.35 in AFR.
  Both cohorts use the same MAF to isolate sample-size and LD effects.
  This does not reproduce real African diversity, admixture, demography or
  the full range of LD found within any ancestry grouping.
- Additive per-allele quantitative effects and independent N(0,1) noise.
- Shared scenario: rsSIM002 = +0.18 and rsSIM006 = -0.16 in both cohorts.
- Heterogeneous scenario: rsSIM002 = +0.18 EUR / +0.09 AFR;
  rsSIM006 = 0 EUR / +0.28 AFR. The latter is AFR-specific by construction.
- Null scenario: all effects zero.
- Twenty replicates per scenario by default: 240 paired datasets. Seed 20261026.
- Within a replicate, the larger cohort is unchanged across sample-size
  comparisons and the smaller samples are nested subsets of the same 2,000.
- Both methods receive exactly the same marginal summary statistics and
  in-sample LD. This is a favourable LD setting, not a test of external-reference
  mismatch. No relatives, missing genotypes, population stratification or
  sample overlap between cohorts are simulated.

## Bayesian hierarchy and likelihood

For SNP j, gamma_j ~ Bernoulli(0.2). When gamma_j=0, both ancestry effects
are zero. Otherwise beta_jk = mu_j + delta_jk, with
mu_j ~ N(0, rho*tau^2) and delta_jk ~ N(0, (1-rho)*tau^2).
tau=0.15 per allele. Primary borrowing rho=0.8; sensitivity rho=0.
The latter removes effect-size correlation but still shares SNP inclusion;
it is **not** a fully separate-ancestry analysis.

The likelihood uses X'X and X'y reconstructed from marginal estimates,
LD and exact in-sample genotype variances. Residual variance is estimated
from the full eight-SNP regression using those sufficient statistics and
then held fixed. Thus this is a plug-in / empirical-Bayes analysis; it does
not integrate uncertainty in residual variance or learn the prior parameters.
It is valid for this small quantitative-trait linear-model simulation,
not a drop-in likelihood for binary/survival GWAS.

All 256 inclusion configurations are enumerated. Gaussian effects are
integrated analytically, yielding exact model probabilities conditional
on the plug-in variance and fixed priors. Reported effects are model-averaged
posterior means, not effects conditional only on the winning model.
No MCMC is required. SDs are posterior standard deviations, not credible
interval endpoints.

The shared inclusion indicator identifies the union of signals. It cannot
assign separate probabilities of a zero effect within each ancestry.
An AFR-only signal can be selected but its effect may be over-shrunk toward
the larger cohort, especially with rho=0.8. That is a limitation to measure.

## Manc-COJO comparison

The real command is saved alongside each dataset. It uses two PLINK BED
reference prefixes and two `.ma` files, ordered EUR then AFR, with
`--cojo-slct --cojo-p 5e-8 --thread-num 1`. Other settings remain the
executable's defaults. All BIM A1 alleles and summary effect alleles are A.
The locally used executable's SHA256 and full configuration are recorded.
Its source/version provenance belongs in `tools/provenance.json`.

Manc-COJO's stepwise meta-analysis criterion and Bayesian shared-inclusion
criterion differ. Neither a PIP threshold nor p<5e-8 guarantees equal
frequentist error rates. Results at PIP >=0.95 and >=0.5 are both reported;
they must not be presented as error-rate-matched comparisons.

## Outcomes

Per replicate: union-of-causal-SNP recall; false discovery proportion
(zero when nothing is selected); probability of any false selection;
AFR-only SNP recovery; mean-squared effect error per ancestry; and AFR
effect error at causal SNPs. For Manc-COJO, unselected SNP effects are set
to zero to define a whole-locus estimator. Bayesian effect MSE uses all
posterior means regardless of the chosen PIP reporting threshold.
Exact causal-SNP recovery is distinct from recovering an LD-equivalent
association signal: COJO can choose a correlated proxy. Therefore this
metric alone is not a fair complete assessment of COJO's stated objective.
Coarse causal-block recall (any selected variant in a causal four-SNP block)
is also reported. It is deliberately permissive and is not a formal
LD-thresholded signal-matching criterion or a credible-set coverage metric.

Report paired differences across methods because inputs are identical.
Replicates, rather than SNPs, are the units for Monte Carlo uncertainty.
Twenty null replicates per sample-size setting cannot estimate extremely
small genome-wide false-positive rates. Zero observed false positives is
not proof of calibration. The prior scale is not tuned to favour results.
Runtimes describe this tiny implementation and input/output overhead,
not a scalable performance comparison.

## Files and execution

`scripts/simulate_compare.py` needs existing Python plus NumPy; it installs
nothing. `scripts/validate.py` checks BED encoding, summary-data likelihood
reconstruction, an analytic posterior calculation and ancestry symmetry.
`scripts/report.R` uses existing base R for tables and figures.
`scripts/build_word_report.R` builds the detailed Word report, tables and
ggplot2 figures. It requires R packages officer, flextable, ggplot2, jsonlite,
xml2 and zip. Install dependencies in your own chosen R library if needed;
the scripts themselves do not install software. The original run used R
4.5.1 for Word authoring and R 4.3.3 for the initial base-R summaries.

From this project directory, with your available interpreter:

```bash
python scripts/validate.py
python scripts/simulate_compare.py --replicates 20 --manc tools/manc_cojo
Rscript scripts/report.R
Rscript scripts/build_word_report.R
python scripts/audit_results.py
```

Use the Linux executable with Linux Python, or Windows executable with
Windows Python. Do not mix the two. For a Bayesian-only preparation run,
omit `--manc`; the status will explicitly read `BAYES_ONLY_MANC_NOT_RUN`.
Do not call that run a completed method comparison.
After obtaining an approved executable, `scripts/run_manc_existing.py --manc
tools/manc_cojo_win.exe` adds Manc-COJO to the already generated paired datasets
without rerunning the Bayesian model. Rerun `scripts/report.R` afterwards.

Results are synthetic. No AWI-Gen or other real participant data are used.
This is a new project; earlier assessment and GATE files are left unchanged.

## Sources

- Manc-COJO paper: https://doi.org/10.1038/s41588-026-02799-1
- Official tutorial: https://light156.github.io/multi-ancestry-COJO-docs/tutorial/
- Official options: https://light156.github.io/multi-ancestry-COJO-docs/options/
- Official code: https://github.com/light156/multi-ancestry-COJO
- Related established Bayesian method MESuSiE:
  https://doi.org/10.1038/s41588-023-01604-7

The custom model here is not MESuSiE and is not claimed to reproduce it.

## Runtime used for this run

The approved official Windows Manc-COJO v1.1.0 executable is stored under
`tools/manc_cojo_win.exe` (about 2.1 MB); `tools/provenance.json` records
the official URL and SHA256. Python/NumPy already existed in the bundled
Windows runtime. Base R already existed in WSL Ubuntu-24.04. Nothing was
installed into a system environment and no cluster jobs were submitted.

Use Python 3.12 with the NumPy version in requirements.txt. On Windows use
`--manc tools/manc_cojo_win.exe` instead of the Linux executable name shown
above. Invoke your installed Rscript for the R report scripts.
Re-running simulation commands
replaces this project's generated datasets/results; it does not modify
earlier GATE or assessment projects.

## Repository contents

The repository includes source code, a portable study configuration, compact
replicate and variant results, R-generated figures, software provenance and
the Word report. Simulated genotype inputs, detailed per-run logs, downloaded
executables and machine-specific paths are excluded by .gitignore. They remain
available locally and the inputs can be regenerated using the fixed seed.
`checksums_public.sha256` records the published source and artifact files.
No real participant data or credentials are included.
