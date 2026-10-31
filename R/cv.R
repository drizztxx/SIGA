#' Cross-validation of the SIGA procedure
#'
#' Estimates how the subgroups identified by SIGA perform in new patients by
#' (repeated) k-fold cross-validation of the complete procedure: [siga()] is run
#' on every training part and the subgroup found is evaluated on the held-out
#' fold. This is the procedure used for the case study of the SIGA paper
#' (section 4.2), where 5-fold cross-validation was repeated 100 times, and can
#' be used to compare tuning parameters.
#'
#' @inheritParams siga_search
#' @param k number of folds.
#' @param repeats number of repetitions with different fold assignments.
#' @param ... arguments passed to [siga()] (e.g. `eta`, `theta`, `n_resample`,
#'   `control`, `control_arm`, `endpoint`, `direction`).
#' @param n_cores number of cores; the folds are distributed over the cores and
#'   every [siga()] call runs on one core.
#' @param verbose print progress messages.
#' @return An object of class `siga_cv`; its `results` data frame has one row
#'   per fold with the identified subgroup, its proportion in the held-out fold,
#'   the treatment effect estimate (hazard ratio, odds ratio or mean
#'   difference) in the held-out subgroup, the one-sided p-value of the
#'   treatment effect in the subgroup and the one-sided p-value of the
#'   treatment-by-subgroup interaction. `summary()` gives the proportion of
#'   folds with a subgroup and the medians of these statistics over the folds
#'   with a subgroup.
#' @seealso [siga()]
#' @examples
#' \donttest{
#' dat <- siga_simulate(300, model = "B1", seed = 12)
#' cv <- siga_cv(Surv(time, status) ~ trt | ., data = dat, k = 3, n_resample = 5,
#'               control = siga_control(pop_size = 20, max_gen = 20, n_subpop = 2), seed = 3)
#' summary(cv)
#' }
#' @export
siga_cv <- function(formula, data, k = 5, repeats = 1, ..., n_cores = 1, seed = NULL,
                    verbose = TRUE) {
  dots <- list(...)
  if (any(c("n_cores", "seed", "verbose") %in% names(dots)))
    stop("'n_cores', 'seed' and 'verbose' are arguments of siga_cv() itself.", call. = FALSE)
  check_number(k, "k", lower = 2, integer = TRUE)
  check_number(repeats, "repeats", lower = 1, integer = TRUE)
  n_cores <- check_cores(n_cores)
  pe <- prepare_endpoint(formula, data, dots$control_arm, dots$endpoint, dots$direction,
                         quiet = !verbose)
  data <- as.data.frame(data)[pe$rows, , drop = FALSE]
  n <- nrow(data)
  if (n < 2 * k) stop("Too few observations for ", k, " folds.", call. = FALSE)
  seed <- resolve_seed(seed)
  streams <- make_streams(repeats * (k + 1L), seed)
  folds <- lapply(seq_len(repeats), function(r)
    with_rng_state(streams[[r]], sample(rep_len(seq_len(k), n))))
  tasks <- lapply(seq_len(repeats * k), function(i)
    list(rep = (i - 1L) %/% k + 1L, fold = (i - 1L) %% k + 1L))
  vmessage(verbose, sprintf("SIGA cross-validation: %d x %d-fold on %d core(s)", repeats, k, n_cores))
  t0 <- proc.time()[["elapsed"]]
  res <- siga_lapply(tasks, cv_task, formula = formula, data = data, folds = folds, dots = dots,
                     .streams = streams[repeats + seq_len(repeats * k)], .n_cores = n_cores)
  structure(list(results = do.call(rbind, res), k = k, repeats = repeats, seed = seed,
                 info = pe$info, call = match.call(),
                 elapsed = proc.time()[["elapsed"]] - t0),
            class = "siga_cv")
}

cv_task <- function(task, formula, data, folds, dots) {
  f <- folds[[task$rep]]
  train <- data[f != task$fold, , drop = FALSE]
  test <- data[f == task$fold, , drop = FALSE]
  row <- data.frame(rep = task$rep, fold = task$fold, n_train = nrow(train), n_test = nrow(test),
                    subgroup = FALSE, rule = NA_character_, proportion = NA_real_,
                    estimate = NA_real_, lower = NA_real_, upper = NA_real_,
                    p_treatment = NA_real_, p_interaction = NA_real_, mean_test_ratio = NA_real_,
                    error = NA_character_, stringsAsFactors = FALSE)
  args <- c(list(formula = formula, data = train, n_cores = 1L,
                 seed = sample.int(.Machine$integer.max, 1L), verbose = FALSE), dots)
  fit <- tryCatch(suppressMessages(do.call(siga, args)), error = function(e) e)
  if (inherits(fit, "error")) {
    row$error <- conditionMessage(fit)
    return(row)
  }
  row$mean_test_ratio <- fit$decision$mean_ratio
  if (!fit$subgroup) return(row)
  row$subgroup <- TRUE
  row$rule <- format(fit$rule)
  eff <- tryCatch(subgroup_effect(fit, data = test), error = function(e) NULL)
  if (is.null(eff)) return(row)
  sub <- eff$table[eff$table$group == "subgroup", ]
  row$proportion <- sub$n / eff$table$n[1L]
  row$estimate <- sub$estimate
  row$lower <- sub$lower
  row$upper <- sub$upper
  row$p_treatment <- sub$p_benefit
  row$p_interaction <- eff$interaction$p_benefit
  row
}

#' @export
print.siga_cv <- function(x, ...) {
  r <- x$results
  cat(sprintf("SIGA %d x %d-fold cross-validation (%s endpoint)\n", x$repeats, x$k, x$info$endpoint))
  cat(sprintf("  subgroup identified in %d of %d folds\n", sum(r$subgroup), nrow(r)))
  if (any(!is.na(r$error))) cat(sprintf("  %d fold(s) failed\n", sum(!is.na(r$error))))
  s <- summary(x)
  if (any(r$subgroup)) {
    cat("  medians over folds with a subgroup (held-out data):\n")
    print(s$medians, row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.siga_cv <- function(object, ...) {
  r <- object$results[object$results$subgroup, , drop = FALSE]
  med <- function(v) if (length(v) && any(!is.na(v))) stats::median(v, na.rm = TRUE) else NA_real_
  medians <- data.frame(proportion = med(r$proportion), estimate = med(r$estimate),
                        p_treatment = med(r$p_treatment), p_interaction = med(r$p_interaction))
  rules <- if (nrow(r)) sort(table(r$rule), decreasing = TRUE) else table(character())
  structure(list(discovery = mean(object$results$subgroup), medians = medians, rules = rules,
                 results = object$results), class = "summary.siga_cv")
}

#' @export
print.summary.siga_cv <- function(x, ...) {
  cat(sprintf("Proportion of folds with a subgroup: %.2f\n", x$discovery))
  cat("Medians over folds with a subgroup (held-out data):\n")
  print(x$medians, row.names = FALSE)
  if (length(x$rules)) {
    cat("Most frequent subgroups:\n")
    print(utils::head(x$rules, 5L))
  }
  invisible(x)
}
