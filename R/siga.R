#' Subgroup identification based on a genetic algorithm
#'
#' Runs the complete SIGA procedure: candidate prescreening, the resampling
#' scheme that controls multiplicity and complexity, and the final
#' multi-population genetic algorithm search on the predictors that were
#' selected frequently during resampling.
#'
#' @details
#' **Resampling.** `n_resample` training subsamples of `resample_prop` (63.2%)
#' of the patients are drawn without replacement (`resample_method = "exact"`)
#' or by independent Bernoulli inclusion (`"bernoulli"`, the behaviour of the
#' original code). In every training subsample the predictors are prescreened
#' (`screen`, `n_candidates`) and the genetic algorithm is run with
#' `resample_subpop` sub-populations. The best subgroup is evaluated by its
#' treatment effect in the held-out patients.
#'
#' **Decision rule.** Let \eqn{r} be the treatment effect ratio of the complete
#' data and \eqn{\bar r^*} the mean ratio of the resampled subgroups in the
#' held-out data, where lower values indicate more benefit. SIGA reports
#' "no subgroup" unless \eqn{\bar r^* < 1} and \eqn{\bar r^* \le (1-\eta) r}
#' (multiplicity control). Only predictors used by the best subgroups of at
#' least a proportion \eqn{\theta} of the resamples enter the final search
#' (complexity control); if there is none, SIGA also reports "no subgroup".
#' The ratio is the hazard ratio for time-to-event endpoints (as in the SIGA
#' paper), the inverse odds ratio for binary endpoints and
#' \eqn{\exp(-\pi d / \sqrt 3)} for continuous endpoints, where \eqn{d} is the
#' standardised mean difference (mean difference divided by the pooled
#' standard deviation of the complete data; the logistic approximation that
#' converts \eqn{d} to a log odds ratio). With \eqn{\eta = 0.2} a continuous
#' subgroup must improve the standardised effect by about 0.12 standard
#' deviations. In a held-out set the ratio is set to 1 when the subgroup or its
#' complement is empty and to missing when it is not estimable or extreme
#' (above 100 or below 0.001).
#'
#' **Final search.** The genetic algorithm with `control$n_subpop` isolated
#' sub-populations is applied to the complete data with the selected predictors.
#' The best individual of all generations and sub-populations defines the
#' subgroup. Thresholds are reported as the boundary values observed in the
#' subgroup and redundant conditions are removed.
#'
#' **Defaults.** The resampling defaults are those of the original SIGA code:
#' 10 resamples, each searched with as many sub-populations as the final search
#' (`control$n_subpop`, 5). With them the package reproduces the false
#' discovery rates of the simulation study of the SIGA paper, for which the
#' recommended tuning parameters `eta` 0.2 (0.15 to 0.25) and `theta` 0.5 (0.5
#' to 0.8) were derived (paper, section 2.3.6). In the same simulations with
#' binary endpoints (`siga_simulate(endpoint = "binary")`) the false discovery
#' rate was 8.5% with `eta = 0.2` and 5% to 6% with `eta = 0.25`, which is
#' therefore the default for binary endpoints. The genetic algorithm uses
#' population size 50, 50 generations and 5 sub-populations. More resamples make
#' the decision less dependent on the seed but stricter, because the mean
#' held-out ratio varies less: with 100 resamples false discoveries are rare at
#' `eta = 0.2`, and `eta` has to be lowered (to about 0.05) for error rates
#' similar to the default. One sub-population per resample, as described in the
#' paper, is faster but finds the predictive biomarkers less often; see
#' `vignette("siga-method")`.
#'
#' @inheritParams siga_search
#' @param eta multiplicity control parameter \eqn{\eta} in \[0, 1): the required
#'   relative improvement of the mean held-out effect ratio over the overall
#'   ratio. Default (`NULL`): 0.2, or 0.25 for binary endpoints (see Details).
#' @param theta complexity control parameter \eqn{\theta} in (0, 1]: minimum
#'   proportion of resamples in which a predictor must be selected. Default 0.5.
#' @param n_resample number of resamples \eqn{q}. Default 10, as in the
#'   original code (see Details).
#' @param screen prescreening method used in every resample: `"lrt"`
#'   (interaction likelihood-ratio test, default), `"rf"` (random forests) or
#'   `"none"`; see [siga_screen()].
#' @param n_candidates number of candidates kept by the prescreening (default:
#'   \eqn{\lceil p/2 \rceil} if \eqn{p \le 20}, else 10).
#' @param screen_args list of further arguments for the random forest
#'   screening.
#' @param resample_prop proportion of patients in every training subsample.
#' @param resample_method `"exact"` (default) or `"bernoulli"`.
#' @param resample_subpop number of sub-populations of the genetic algorithm in
#'   every resample. Default `control$n_subpop` (5), as in the original code.
#'   The SIGA paper describes one sub-population per resample.
#' @param control genetic algorithm settings, see [siga_control()].
#' @param n_cores number of cores. Resamples (and the sub-populations of the
#'   final search) are distributed over the cores.
#' @param seed integer seed; see [siga_search()].
#' @param verbose print progress messages.
#' @return An object of class `siga` with components
#'   \item{subgroup}{logical: was a subgroup identified?}
#'   \item{rule}{the identified subgroup ([siga_rule]) or `NULL`.}
#'   \item{decision}{the decision rule: mean held-out ratio, threshold,
#'     selected predictors and the reason for "no subgroup".}
#'   \item{resamples}{one row per resample: status, sizes, objective, effect
#'     ratios in the training and test data and the best subgroup.}
#'   \item{frequency}{selection frequency of every predictor.}
#'   \item{screening}{prescreening of the complete data (for information).}
#'   \item{final}{the final search (see [siga_search()]) or `NULL`.}
#'   \item{effects}{treatment effects in the subgroup, its complement and the
#'     interaction, see [subgroup_effect()].}
#'   \item{overall_ratio}{the effect ratio of the complete data.}
#' @references Zhao L, Zhang W, Wu Y, Cao L, Wang L, Li K (2024). A promising
#'   subgroup identification method based on a genetic algorithm for censored
#'   survival data. *Journal of Biopharmaceutical Statistics* 34(1), 55--77.
#'   \doi{10.1080/10543406.2023.2170397}
#' @seealso [update.siga()], [siga_cv()], [subgroup_effect()], [predict.siga()],
#'   [plot.siga()]
#' @examples
#' \donttest{
#' dat <- siga_simulate(300, model = "B1", seed = 10)
#' fit <- siga(Surv(time, status) ~ trt | ., data = dat,
#'             control = siga_control(pop_size = 30, max_gen = 30, n_subpop = 2),
#'             seed = 1)
#' fit
#' summary(fit)
#' ## stricter complexity control without repeating the resampling
#' update(fit, theta = 0.8)
#' }
#' @export
siga <- function(formula, data, control_arm = NULL, endpoint = NULL, direction = NULL,
                 eta = NULL, theta = 0.5, n_resample = 10,
                 screen = c("lrt", "rf", "none"), n_candidates = NULL, screen_args = list(),
                 resample_prop = 0.632, resample_method = c("exact", "bernoulli"),
                 resample_subpop = control$n_subpop, control = siga_control(),
                 predictor_types = NULL, n_cores = 1, seed = NULL, verbose = TRUE) {
  screen <- match.arg(screen)
  resample_method <- match.arg(resample_method)
  if (!is.null(eta)) {
    check_number(eta, "eta", lower = 0, upper = 1)
    if (eta >= 1) stop("'eta' must be smaller than 1.", call. = FALSE)
  }
  check_number(theta, "theta", lower = 0, upper = 1)
  if (theta <= 0) stop("'theta' must be positive.", call. = FALSE)
  check_number(n_resample, "n_resample", lower = 1, integer = TRUE)
  check_number(resample_prop, "resample_prop", lower = 0.1, upper = 0.95)
  ## checked before resample_subpop, whose default is control$n_subpop
  if (!inherits(control, "siga_control"))
    stop("'control' must be created by siga_control().", call. = FALSE)
  check_number(resample_subpop, "resample_subpop", lower = 1, integer = TRUE)
  if (!is.null(n_candidates)) check_number(n_candidates, "n_candidates", lower = 1, integer = TRUE)
  n_cores <- check_cores(n_cores)
  t0 <- proc.time()[["elapsed"]]

  prep <- prepare_siga_data(formula, data, control_arm, endpoint, direction, predictor_types,
                            quiet = !verbose)
  eta <- eta %||% default_eta(prep$ep$type)
  seed <- resolve_seed(seed)
  q <- as.integer(n_resample)
  streams <- make_streams(q + 1L + control$n_subpop, seed)
  settings <- list(eta = eta, theta = theta, n_resample = q, screen = screen,
                   n_candidates = n_candidates, screen_args = screen_args,
                   resample_prop = resample_prop, resample_method = resample_method,
                   resample_subpop = as.integer(resample_subpop), control = control,
                   predictors = names(prep$types))

  overall_ratio <- exp(-effect_theta(prep$ep, rep(TRUE, prep$ep$n)))
  screening <- with_rng_state(streams[[q + 1L]],
                              screen_predictors(prep, settings$predictors, screen, n_candidates,
                                                screen_args))
  vmessage(verbose, sprintf("SIGA: %s endpoint, %d patients, %d predictors; %d resamples on %d core(s)",
                            prep$info$endpoint, prep$ep$n, length(prep$types), q, n_cores))
  res <- siga_lapply(seq_len(q), one_resample, prep = prep, settings = settings,
                     .streams = streams[seq_len(q)], .n_cores = n_cores)
  rtab <- resample_table(res)
  freq <- variable_frequency(res, settings$predictors)

  fit <- structure(list(call = match.call(), formula = formula, data = prep, settings = settings,
                        seed = seed, overall_ratio = overall_ratio, screening = screening,
                        resamples = rtab, resample_rules = lapply(res, `[[`, "rule"),
                        frequency = freq),
                   class = "siga")
  fit <- siga_finalize(fit, eta, theta, n_cores, verbose)
  fit$elapsed <- proc.time()[["elapsed"]] - t0
  vmessage(verbose, sprintf("SIGA finished in %.1f s: %s", fit$elapsed,
                            if (fit$subgroup) format(fit$rule) else "no subgroup"))
  fit
}

## Default multiplicity control: 0.2 as in the SIGA paper; 0.25 for binary
## endpoints, whose false discovery rate is higher at 0.2 (validation, V6).
default_eta <- function(type) if (type == "binary") 0.25 else 0.2

## Apply the decision rule and run the final search if needed. The final
## search always uses the same random number streams, so the result for given
## (eta, theta) does not depend on whether it came from siga() or update().
siga_finalize <- function(fit, eta, theta, n_cores = 1L, verbose = FALSE, previous = NULL) {
  dec <- siga_decision(fit$resamples, fit$frequency, fit$overall_ratio, eta, theta)
  fit$settings$eta <- eta
  fit$settings$theta <- theta
  fit$decision <- dec
  fit$subgroup <- FALSE
  fit$rule <- NULL
  fit$final <- NULL
  fit$effects <- NULL
  if (!dec$subgroup) {
    vmessage(verbose, "SIGA: no subgroup (", dec$reason, ")")
    return(fit)
  }
  control <- fit$settings$control
  if (!is.null(previous) && identical(previous$vars, dec$selected)) {
    final <- previous
  } else {
    vmessage(verbose, sprintf("SIGA: final search with %d sub-population(s) on %s",
                              control$n_subpop, paste(dec$selected, collapse = ", ")))
    q <- fit$settings$n_resample
    streams <- make_streams(q + 1L + control$n_subpop, fit$seed)[q + 1L + seq_len(control$n_subpop)]
    final <- run_search(fit$data, dec$selected, control, streams, n_cores)
  }
  fit$final <- final
  if (final$valid) {
    fit$subgroup <- TRUE
    fit$rule <- final$rule
    fit$effects <- effect_summary(fit$data$ep, final$mask)
  } else {
    fit$decision$subgroup <- FALSE
    fit$decision$reason <- "the final search found no valid subgroup"
  }
  fit
}

#' Re-apply the SIGA decision rule with other tuning parameters
#'
#' Recomputes the multiplicity and complexity decision of a [siga()] fit for new
#' values of `eta` and `theta` using the stored resampling results, and reruns
#' the final search only when the set of selected predictors changes. The
#' result is identical to calling [siga()] again with the same seed.
#'
#' @param object a `siga` fit.
#' @param eta,theta new tuning parameters (default: those of `object`).
#' @param n_cores number of cores for the final search.
#' @param verbose print progress messages.
#' @param ... unused.
#' @return An updated `siga` object.
#' @seealso [siga()]
#' @examples
#' \donttest{
#' dat <- siga_simulate(300, model = "B3", seed = 11)
#' fit <- siga(Surv(time, status) ~ trt | ., data = dat,
#'             control = siga_control(pop_size = 30, max_gen = 30, n_subpop = 2), seed = 2)
#' update(fit, eta = 0.25, theta = 0.7)
#' }
#' @export
update.siga <- function(object, eta = object$settings$eta, theta = object$settings$theta,
                        n_cores = 1, verbose = FALSE, ...) {
  check_number(eta, "eta", lower = 0, upper = 1)
  if (eta >= 1) stop("'eta' must be smaller than 1.", call. = FALSE)
  check_number(theta, "theta", lower = 0, upper = 1)
  if (theta <= 0) stop("'theta' must be positive.", call. = FALSE)
  n_cores <- check_cores(n_cores)
  fit <- siga_finalize(object, eta, theta, n_cores, verbose, previous = object$final)
  fit$call$eta <- eta
  fit$call$theta <- theta
  fit
}
