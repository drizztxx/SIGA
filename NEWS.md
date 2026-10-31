# SIGA 0.1.0

First release as an R package. The original scripts (now in `legacy/` of the
source repository) were rewritten; the genetic algorithm reproduces the
original one operator by operator (verified generation by generation).

## New features

* Formula interface `outcome ~ treatment | predictors` for `siga()`,
  `siga_search()`, `siga_screen()` and `siga_cv()`.
* Continuous endpoints (linear model interaction test) in addition to
  time-to-event and binary endpoints; `direction` sets which direction of the
  treatment effect means benefit.
* `update()` re-applies the decision rule with other `eta` / `theta` without
  repeating the resampling.
* `subgroup_effect()`, `predict()`, `print()`, `summary()` and `plot()`
  methods; `siga_cv()` for repeated k-fold cross-validation.
* `siga_simulate()` generates data from the models of the SIGA paper
  (`legacy = TRUE` reproduces the original simulation program exactly).
* Reproducible parallel computing: results are identical for any number of
  cores.
* About 35 times faster than the original code per genetic algorithm run
  (fast model fitters, no string evaluation, caching of repeated subgroups);
  a complete analysis of 300 patients with the default settings takes about
  20 seconds on one core.
* The resampling defaults are those of the original `SIGA()` function, with
  which the package reproduces the false discovery rates of the simulation
  study of the SIGA paper: 10 resamples (`n_resample`), each searched with the
  5 sub-populations of the final search (`resample_subpop`). Exactly 63.2% of
  the patients are drawn per resample, as described in the paper
  (`resample_method = "bernoulli"` gives the sampling of the original code).
* `eta` defaults to 0.2 as recommended in the paper, and to 0.25 for binary
  endpoints: in the validation simulations the false discovery rate of binary
  endpoints was 8.5% with 0.2 and 5% to 6% with 0.25.

## Bug fixes relative to the original code

* A single continuous biomarker used in a direction (`<=` or `>=`) took its
  threshold from the pooled quantiles of all numeric candidates instead of its
  own quantiles.
* A single binary factor used in a subgroup produced the rule `x == ''`
  (empty subgroup), so such subgroups could never be selected.
* A single categorical factor (when several were candidates) stopped the
  search with an error; factors with exactly three levels could only define
  single-level subsets, with a wrong mapping of the genes.
* The binary prescreening used the deviance instead of the p-value and kept
  the least predictive biomarkers.
* The best-subgroup summary could produce vacuous conditions such as
  `x >= a | x <= b`.
* The final generation was ignored when choosing the best individual.
