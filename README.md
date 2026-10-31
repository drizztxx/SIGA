# SIGA

**Subgroup Identification based on a Genetic Algorithm**

SIGA searches two-arm randomised clinical trial data for a subgroup of
patients, defined by a combination of biomarker conditions, in which the
treatment effect is much larger than in the rest of the population. Candidate
subgroups are encoded as real-valued chromosomes, scored by the
treatment-by-subgroup interaction p-value and evolved with the genetic
algorithm operators of [gatbxr](https://github.com/drizztxx/gatbxr) in several
isolated sub-populations. A resampling scheme controls multiplicity (false
subgroup discovery, parameter `eta`) and complexity (noise biomarkers,
parameter `theta`).

* Endpoints: time-to-event (Cox), binary (logistic) and continuous (linear).
* Biomarkers: continuous, ordinal, binary and categorical.
* Prescreening of candidate biomarkers, cross-validation, prediction for new
  patients, plots, parallel computing with reproducible results, and the
  simulation models of the SIGA paper.

## Installation

```r
# install.packages("remotes")
remotes::install_github("drizztxx/gatbxr")   # genetic algorithm toolbox
remotes::install_github("drizztxx/SIGA")
```

## Example

```r
library(SIGA)

## simulated trial: model B3 of the SIGA paper, true subgroup x4 > 1 & x5 = 1
dat <- siga_simulate(300, model = "B3", seed = 1)

fit <- siga(Surv(time, status) ~ trt | ., data = dat,
            eta = 0.2, theta = 0.5, n_cores = 4, seed = 1)
fit                      # decision, selected biomarkers, subgroup and effects
summary(fit)
plot(fit, "effect")      # Kaplan-Meier curves in the subgroup and its complement

update(fit, eta = 0.25, theta = 0.7)     # other tuning parameters, no resampling
predict(fit, newdata = dat)               # subgroup membership
subgroup_effect(fit, data = siga_simulate(2000, model = "B3", seed = 2))
```

Binary and continuous endpoints use the same interface, e.g.
`siga(response ~ arm | ., data, control_arm = "placebo")`.

## Documentation

* `vignette("SIGA")`: getting started.
* `vignette("siga-method")`: the method, its parameters and the differences
  from the original code.
* `vignette("siga-endpoints")`: binary and continuous endpoints.
* `vignette("siga-simulation")`: reproducing the simulation study.

## Reference

Zhao L, Zhang W, Wu Y, Cao L, Wang L, Li K (2024). A promising subgroup
identification method based on a genetic algorithm for censored survival data.
*Journal of Biopharmaceutical Statistics* 34(1), 55–77.
[doi:10.1080/10543406.2023.2170397](https://doi.org/10.1080/10543406.2023.2170397)

`citation("SIGA")` gives the reference in BibTeX format.

## License

MIT. The genetic algorithm operators come from the
[gatbxr](https://github.com/drizztxx/gatbxr) package (GPL >= 2), an R port of
the Sheffield Genetic Algorithm Toolbox for MATLAB.
