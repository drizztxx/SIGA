## Multi-population genetic algorithm search ----------------------------------

prep_subset <- function(prep, idx) {
  list(ep = ep_subset(prep$ep, idx), X = prep$X[idx, , drop = FALSE], types = prep$types,
       info = prep$info)
}

## Combine the runs of several isolated sub-populations. The best subgroup has
## the smallest objective; ties are broken by the larger subgroup, then by the
## smaller number of conditions, then by the sub-population index.
combine_runs <- function(runs, ep, vars) {
  valid <- vapply(runs, `[[`, logical(1), "valid")
  obj <- vapply(runs, `[[`, numeric(1), "objective")
  size <- bal <- ncond <- rep(NA_real_, length(runs))
  rules <- rep(NA_character_, length(runs))
  for (i in which(valid)) {
    m <- runs[[i]]$mask
    n1 <- sum(m & ep$z == 1)
    size[i] <- sum(m)
    bal[i] <- 2 * min(n1, size[i] - n1)
    ncond[i] <- length(runs[[i]]$rule$conditions)
    rules[i] <- format(runs[[i]]$rule)
  }
  cand <- data.frame(subpop = seq_along(runs), objective = obj, size = size,
                     proportion = size / ep$n, n_conditions = ncond, rule = rules,
                     stringsAsFactors = FALSE)
  best <- if (any(valid)) {
    ok <- which(valid)
    ok[order(obj[ok], -bal[ok], ncond[ok], ok)[1L]]
  } else NA_integer_
  history <- do.call(cbind, lapply(runs, `[[`, "history"))
  if (!is.null(history)) colnames(history) <- paste0("subpop", seq_along(runs))
  list(valid = any(valid), best = best,
       rule = if (any(valid)) runs[[best]]$rule else NULL,
       mask = if (any(valid)) runs[[best]]$mask else NULL,
       objective = if (any(valid)) obj[best] else min(obj),
       candidates = cand, history = history, vars = vars,
       rules = lapply(runs, `[[`, "rule"),
       n_eval = sum(vapply(runs, `[[`, integer(1), "n_eval")))
}

## Search on the predictors `vars` of `prep` with one stream per sub-population.
run_search <- function(prep, vars, control, streams, n_cores = 1L) {
  enc <- encode_predictors(prep$X, vars, prep$types)
  ## sequential sub-populations share one objective cache (each parallel
  ## worker receives its own copy)
  cache <- new_objective_cache()
  runs <- siga_lapply(seq_len(control$n_subpop),
                      function(i, ep, X, enc, control, cache)
                        ga_run(ep, X, enc, control, cache = cache),
                      ep = prep$ep, X = prep$X, enc = enc, control = control, cache = cache,
                      .streams = streams, .n_cores = n_cores)
  combine_runs(runs, prep$ep, vars)
}

#' Genetic algorithm search for a subgroup
#'
#' Runs the SIGA genetic algorithm (one or several isolated sub-populations) on
#' a given set of candidate predictors, without prescreening or resampling.
#' This is the search step of [siga()]; it always returns the best subgroup
#' found, which is subject to selection bias, and should be used for
#' exploration or with predictors chosen in advance.
#'
#' @param formula a formula `outcome ~ treatment | predictors`. The outcome is
#'   `Surv(time, status)` for time-to-event endpoints, a binary variable
#'   (logical, two-level factor or 0/1) or a numeric variable for continuous
#'   endpoints. `treatment` is a variable with two values. `predictors` are
#'   column names separated by `+`, or `.` for all remaining columns
#'   (`. - id` excludes `id`).
#' @param data a data frame.
#' @param control_arm value of the treatment variable that identifies the
#'   control arm. By default 0 or `FALSE` for numeric or logical treatment
#'   variables and otherwise the first level (with a message).
#' @param endpoint `"tte"`, `"binary"` or `"continuous"`; detected from the
#'   outcome by default.
#' @param direction which direction of the treatment effect indicates benefit:
#'   `"lower"` hazard (default for time-to-event), `"higher"` probability of
#'   the event level (default for binary; the event level is `TRUE`, 1 or the
#'   second factor level) or `"higher"` mean (default for continuous).
#' @param predictor_types optional named character vector overriding the
#'   predictor types (`"numeric"`, `"ordinal"`, `"binary"`, `"categorical"`).
#'   By default numeric columns are `"numeric"`, ordered factors `"ordinal"`,
#'   and factors, logical and character columns `"binary"` or `"categorical"`
#'   depending on the number of levels.
#' @param control genetic algorithm settings, see [siga_control()].
#' @param n_cores number of cores used to run the sub-populations in parallel.
#' @param seed integer seed. With the same seed results are identical for any
#'   `n_cores`. By default a seed is drawn from the current random number
#'   generator, so [set.seed()] also makes the run reproducible.
#' @return An object of class `siga_search` with components
#'   \item{rule}{the best subgroup ([siga_rule] object, `NULL` if no valid
#'     subgroup was found).}
#'   \item{objective}{its objective value.}
#'   \item{candidates}{the best subgroup of every sub-population.}
#'   \item{history}{best objective value per generation and sub-population.}
#'   \item{effects}{treatment effects in the subgroup and its complement, see
#'     [subgroup_effect()].}
#' @seealso [siga()], [siga_control()]
#' @examples
#' dat <- siga_simulate(300, model = "B3", seed = 2)
#' fit <- siga_search(Surv(time, status) ~ trt | x2 + x4 + x5 + x9, data = dat,
#'                    control = siga_control(pop_size = 30, max_gen = 20, n_subpop = 2),
#'                    seed = 1)
#' fit
#' @export
siga_search <- function(formula, data, control_arm = NULL, endpoint = NULL, direction = NULL,
                        predictor_types = NULL, control = siga_control(), n_cores = 1,
                        seed = NULL) {
  if (!inherits(control, "siga_control")) stop("'control' must be created by siga_control().",
                                               call. = FALSE)
  n_cores <- check_cores(n_cores)
  prep <- prepare_siga_data(formula, data, control_arm, endpoint, direction, predictor_types)
  seed <- resolve_seed(seed)
  t0 <- proc.time()[["elapsed"]]
  res <- run_search(prep, names(prep$types), control, make_streams(control$n_subpop, seed), n_cores)
  res$effects <- if (res$valid) effect_summary(prep$ep, res$mask) else NULL
  res$formula <- formula
  res$data <- prep
  res$control <- control
  res$seed <- seed
  res$elapsed <- proc.time()[["elapsed"]] - t0
  res$call <- match.call()
  class(res) <- "siga_search"
  res
}
