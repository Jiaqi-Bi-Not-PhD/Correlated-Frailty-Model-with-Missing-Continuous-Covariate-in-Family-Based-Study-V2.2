# Continuous-covariate correlated-frailty simulations

This repository contains the R code needed to reproduce the synthetic
continuous-covariate simulation study in *Correlated Frailty Model with Missing
Continuous Covariate in Family-Based Study*.

The publication scope is limited to the case in which the continuous polygenic
risk score (PRS) is missing. It includes the corresponding full-data benchmark,
continuous-case complete-case analysis (CCA), continuous-case MI-SMCFCS,
C-O-PDMI, and C-R-PDMI. 

Reviewers generate all results locally with the R scripts described below.

## Simulation grid

Every published method/scenario cell uses 1,000 attempted simulation
replicates.

| Analysis | Replicates | Imputations |
|---|---:|---:|
| No-missing full-data benchmark | 1,000 per frailty variance | not applicable |
| Complete-case analysis (CCA) | 1,000 per scenario | not applicable |
| MI-SMCFCS | 1,000 per scenario | 20 |
| C-O-PDMI | 1,000 per scenario | 20 |
| C-R-PDMI | 1,000 per scenario | 20 |

The scenario grid uses:

- 498 generated families;
- frailty variances \(\sigma_u^2\in\{0.2,0.5\}\);
- continuous-PRS missingness rates of 20%, 50%, and 80%;
- true parameters
  \((\log\rho,\log\lambda,\beta_b,\beta_c)=(0.804,4.71,2.2,1.0)\);
- PRS variance 0.1;
- 10 PDMI iterations and 20 retained PDMI imputations;
- 20 MI-SMCFCS imputations.

Seeds are deterministic functions of the replicate number, scenario, and
method. A failed model fit is retained with its convergence flag and failure
reason.

## Repository contents

```text
.
├── code/
│   ├── pdmi-v2.2/
│   │   ├── Continuous Only/
│   │   │   └── run_continuous_pdmi_simulation.R
│   │   ├── Shared/
│   │   │   ├── config.R
│   │   │   ├── continuous_missingness.R
│   │   │   ├── data_generation.R
│   │   │   ├── frailtypack_analysis.R
│   │   │   ├── numerics.R
│   │   │   ├── packages_and_sources.R
│   │   │   ├── pdmi_continuous.R
│   │   │   ├── pdmi_diagnostics.R
│   │   │   ├── pdmi_exact_slice_continuous.R
│   │   │   ├── results_and_summaries.R
│   │   │   └── source_all.R
│   │   ├── dependencies/          # required family-generation R helpers
│   │   └── Misc/
│   │       ├── check_source_all.R
│   │       ├── summarize_results.R
│   │       ├── test_exact_slice_pdmi.R
│   │       ├── test_pdmi_diagnostics.R
│   │       └── test_pdmi_pooling.R
│   └── comparators-may23/
│       ├── Continuous Only/
│       │   └── run_continuous_simulation.R
│       ├── No Missing/
│       │   └── run_no_missing_benchmark.R
│       ├── Shared/
│       │   ├── config.R
│       │   ├── continuous_missingness_and_cca.R
│       │   ├── data_generation.R
│       │   ├── frailtypack_analysis.R
│       │   ├── numerics.R
│       │   ├── packages_and_sources.R
│       │   ├── results_and_summaries.R
│       │   ├── smcfcs_continuous_comparator.R
│       │   └── source_all.R
│       ├── dependencies/          # required family-generation R helpers
│       └── Misc/
│           ├── check_source_all.R
│           └── summarize_results.R
├── scripts/
│   ├── install_dependencies.R
│   ├── validate_repository.R
│   └── run_simulation_grid.R
├── .gitignore
└── README.md
```

Only the files shown above, the R helper files below the two `dependencies/`
directories, and the two repository-level documentation/control files are
part of the GitHub publication set. Dedicated binary-missingness, joint-
missingness, legacy MI-Cong, cluster, and Slurm modules are excluded.

## R dependencies

The code was developed with R 4.5.0 and later validated with R 4.5.3. The two
critical package versions were:

- `frailtypack` 2.13;
- `smcfcs` 2.0.2.

The scripts also use:

- `Matrix`;
- `parallel`;
- `survival`;
- `kinship2`;
- `MASS`;
- `truncnorm`.

Their required transitive dependencies are installed automatically by R.
Installing `frailtypack` may require C, C++, and Fortran compilers together
with BLAS and LAPACK development libraries. Please contact the maintainer of `frailtypack`
for any download specifications. You may find the following link useful: 
https://cran.r-project.org/web/packages/frailtypack/index.html.

Parallel execution is controlled by `SIM_CORES`. On Windows, the scripts use
one core because `parallel::mclapply()` does not provide multicore execution
there.

## Step 1: obtain the code

Download or clone the repository, open a terminal, and change to the repository
root—the directory containing this `README.md`.

Multi-line commands below use the macOS/Linux continuation character `\` for
readability. In Windows PowerShell, enter the same command and options on one
line.

## Step 2: install the R packages

The following creates a project-local R library and installs missing packages
from CRAN:

```bash
export R_LIBS_USER="$PWD/.R-library"
Rscript scripts/install_dependencies.R
```

In Windows PowerShell, set the same library with:

```powershell
$env:R_LIBS_USER = "$PWD/.R-library"
Rscript scripts/install_dependencies.R
```

The installer reports the installed package versions. If CRAN supplies newer
versions, it warns that the study was validated with `frailtypack` 2.13 and
`smcfcs` 2.0.2.

## Step 3: validate the R code

Run the structural and configuration check:

```bash
Rscript scripts/validate_repository.R
```

To require all packages and the two tested method-critical versions:

```bash
Rscript scripts/validate_repository.R --strict-packages
```

The following optional checks load the source modules and run the focused PDMI
tests without running the simulation grid:

```bash
Rscript code/pdmi-v2.2/Misc/check_source_all.R
Rscript code/pdmi-v2.2/Misc/test_exact_slice_pdmi.R
Rscript code/pdmi-v2.2/Misc/test_pdmi_pooling.R
Rscript code/pdmi-v2.2/Misc/test_pdmi_diagnostics.R
Rscript code/comparators-may23/Misc/check_source_all.R
```

## Step 4: run the full simulations

The full study is computationally expensive. Run the three components
separately so that progress and failures are easy to inspect. Replace
`--cores=4` with a suitable number for the computer being used.

### 4a. No-missing benchmark

This runs 1,000 replicates for each of the two frailty variances:

```bash
Rscript scripts/run_simulation_grid.R \
  --component=no-missing \
  --n=1000 \
  --families=498 \
  --cores=4 \
  --run-prefix=reviewer
```

### 4b. CCA and MI-SMCFCS

This runs CCA and MI-SMCFCS for all six combinations of frailty variance and
missingness rate. MI-SMCFCS uses 20 imputations:

```bash
Rscript scripts/run_simulation_grid.R \
  --component=comparators \
  --n=1000 \
  --families=498 \
  --cores=4 \
  --m-smcfcs=20 \
  --run-prefix=reviewer
```

### 4c. C-O-PDMI and C-R-PDMI

This runs both V2.2 PDMI methods for all six combinations of frailty variance
and missingness rate. Each cell uses 20 imputations and 10 PDMI iterations:

```bash
Rscript scripts/run_simulation_grid.R \
  --component=pdmi \
  --n=1000 \
  --families=498 \
  --cores=4 \
  --m-pdmi=20 \
  --pdmi-numit=10 \
  --run-prefix=reviewer
```

The driver processes scenario cells sequentially. Within each cell, replicate
computations use the requested number of cores. Existing replicate files are
skipped, so rerunning the same command resumes an interrupted run. Add
`--force=true` only when results should be recomputed deliberately.

## Step 5: find and inspect the outputs

By default, locally generated files are written below:

```text
reproduced-results/
├── no-missing/
├── comparators/
└── pdmi/
```

Each component has a run-label directory containing one RDS file per replicate.
After all scenario cells for a component finish, the driver automatically runs
the corresponding `Misc/summarize_results.R` script. Its `summary_tables/`
directory includes:

- parameter estimates and performance summaries;
- penetrance estimates and performance summaries;
- convergence and failure diagnostics;
- manuscript-oriented summary tables.

Penetrance figures are written under the same run-label directory. All
generated output is excluded from version control.

To write output somewhere else, add an absolute path:

```bash
Rscript scripts/run_simulation_grid.R \
  --component=pdmi \
  --n=1000 \
  --families=498 \
  --cores=4 \
  --output=/absolute/path/to/simulation-output
```

## Driver options

| Option | Meaning | Default |
|---|---|---:|
| `--component` | `no-missing`, continuous `comparators`, `pdmi`, or `all` | `all` |
| `--n` | attempted replicates per cell | `1000` |
| `--families` | generated families per replicate | `498` |
| `--cores` | parallel replicate workers | `1` |
| `--m-smcfcs` | MI-SMCFCS imputations | `20` |
| `--m-pdmi` | retained PDMI imputations | `20` |
| `--pdmi-numit` | PDMI iterations | `10` |
| `--run-prefix` | prefix for output run labels | `reviewer` |
| `--output` | output root | `reproduced-results/` |
| `--maxit` | optional frailtypack iteration limit | configuration default |
| `--force` | recompute existing replicates | `false` |
| `--dry-run` | print planned R calls without running them | `false` |

Arguments must use the `--name=value` form shown above.
