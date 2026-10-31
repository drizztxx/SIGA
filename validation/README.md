# Validation of the SIGA package

This directory is excluded from the package build. It contains the scripts
that validate the package against the original SIGA code (`../legacy`) and
against the results published in the SIGA paper. Run the scripts from this
directory with the installed SIGA package (parallel runs start new R
processes that load the installed package).

| Step | Script | What is checked | Runtime |
|---|---|---|---|
| V1 | `../tests/testthat` | unit and regression tests (`devtools::test()`) | < 1 min |
| V2 | `01_component_equivalence.R` | decoding, objective values, prescreening, the complete genetic algorithm trajectory and the simulation program reproduce the original code on identical inputs | 2 min |
| V3 | `02_ga_comparison.R [n]` | search quality of one genetic algorithm run, original code vs package, on the same simulated B3 datasets | 15 min |
| V3 | `02b_siga_comparison.R [workers] [b1]` | the complete SIGA procedure, original code vs package, on the same B1 and B4 datasets (needs the run `04_simulation.R endpoint=tte R=24 s=5 q=10 method=bernoulli sim=legacy` first) | 25 min per dataset (32 datasets: 3.5 h on 4 workers) |
| V4 | `03_radiant_case_study.R <RData> [cores]` | the RADIANT case study of the paper (restricted data, read from a local file; aggregate output only) | 8 min on 8 cores |
| V5 | `04_simulation.R endpoint=tte R=100 cores=8` | simulation study of the paper (Tables 3 to 6), 100 trials per scenario | 1.5 h on 8 cores |
| V6 | `04_simulation.R endpoint=continuous R=100 cores=8`, `04_simulation.R endpoint=binary R=100 cores=8` | the same models for the binary and continuous endpoints | 3 h on 8 cores (continuous 50 min, binary 2 h) |
| | `05_summarise.R` | tables of V5/V6 compared with the published values (`paper_tables.csv`); for runs with 100 resamples, the decisions with fewer resamples | < 1 min |
| | `06_report.Rmd` | validation report (text and tables in `validation-chapter.Rmd`) | < 1 min |

`04_simulation.R` uses the package defaults (`s=5` sub-populations per
resample, `q=10` resamples, `method=exact` sampling of 63.2% of the patients)
and time-to-event data with the observed time min(T, C) (`sim=proper`). Other
settings are run with, e.g., `s=1`, `q=100`, `method=bernoulli` (sampling of the
original code) or `sim=legacy` (original simulation program) and are stored
under their own label (e.g. `tte_s1_q100_legacydata`). The comparison of
settings in the report uses

    Rscript 04_simulation.R endpoint=tte R=50 method=bernoulli
    Rscript 04_simulation.R endpoint=tte R=50 method=bernoulli sim=legacy
    Rscript 04_simulation.R endpoint=tte R=100 s=1 sim=legacy
    Rscript 04_simulation.R endpoint=tte R=100 s=1 q=100 sim=legacy
    Rscript 04_simulation.R endpoint=tte R=50 q=100 method=bernoulli sim=legacy models=A1,A2 sizes=300

(add `cores=8`; about 40 min, 40 min, 25 min, 3 h and 2 h on 8 cores).

`paper_tables.csv` contains the SIGA rows of Tables 3 to 6 of the paper,
parsed from the final manuscript. The simulation trials are written to
`%LOCALAPPDATA%/SIGA-validation/<label>/` (one file per trial, so that a run
can be interrupted and resumed); summaries are written to `results/`.

The RADIANT data were obtained from ClinicalStudyDataRequest.com under a data
sharing agreement. They must not be copied into this repository.
