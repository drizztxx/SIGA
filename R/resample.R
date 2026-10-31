## Resampling scheme (SIGA paper, section 2.3.5) ------------------------------
##
## Step 1  draw a training subsample of 63.2% of the patients without
##         replacement; the remaining 36.8% form the test set;
## Step 2  screen the predictors and run SIGA on the training subsample;
##         evaluate the best subgroup by its treatment effect in the test set;
## Step 3  stop with "no subgroup" unless the mean test effect ratio is below 1
##         and at least a factor (1 - eta) below the overall ratio, and keep the
##         predictors selected in at least a proportion theta of the resamples;
## Step 4  final multi-population search on all data with these predictors.

one_resample <- function(h, prep, settings) {
  n <- prep$ep$n
  train <- if (settings$resample_method == "exact") {
    sort(sample.int(n, round(settings$resample_prop * n)))
  } else {
    which(stats::rbinom(n, 1L, settings$resample_prop) == 1L)
  }
  test <- setdiff(seq_len(n), train)
  out <- list(id = h, status = "failed", n_train = length(train), n_test = length(test),
              objective = NA_real_, rule = NULL, vars = character(), candidates = character(),
              train_ratio = NA_real_, test_ratio = NA_real_, train_prop = NA_real_,
              test_prop = NA_real_)
  if (length(test) < 2L) return(out)
  ptrain <- prep_subset(prep, train)
  ptest <- prep_subset(prep, test)
  if (min(sum(ptrain$ep$z == 1), sum(ptrain$ep$z == 0)) < 2L) return(out)
  vars <- tryCatch(
    screen_predictors(ptrain, settings$predictors, settings$screen, settings$n_candidates,
                      settings$screen_args)$selected,
    error = function(e) NULL)
  if (!length(vars)) return(out)
  out$candidates <- vars
  control <- settings$control
  control$n_subpop <- settings$resample_subpop
  control$keep_history <- FALSE
  enc <- encode_predictors(ptrain$X, vars, ptrain$types)
  objective <- make_objective(ptrain$ep, control$objective, control$min_size)
  cache <- new_objective_cache()
  runs <- lapply(seq_len(control$n_subpop), function(i)
    ga_run(ptrain$ep, ptrain$X, enc, control, objective = objective, cache = cache))
  res <- combine_runs(runs, ptrain$ep, vars)
  out$objective <- res$objective
  if (!res$valid) return(out)
  test_mask <- rule_members(res$rule, ptest$X)
  out$status <- "ok"
  out$rule <- res$rule
  out$vars <- rule_variables(res$rule)
  out$train_ratio <- exp(-effect_theta(ptrain$ep, res$mask))
  out$test_ratio <- test_ratio(ptest$ep, test_mask)
  out$train_prop <- mean(res$mask)
  out$test_prop <- mean(test_mask)
  out
}

resample_table <- function(res) {
  data.frame(
    id = vapply(res, `[[`, numeric(1), "id"),
    status = vapply(res, `[[`, character(1), "status"),
    n_train = vapply(res, `[[`, numeric(1), "n_train"),
    n_test = vapply(res, `[[`, numeric(1), "n_test"),
    objective = vapply(res, `[[`, numeric(1), "objective"),
    train_ratio = vapply(res, `[[`, numeric(1), "train_ratio"),
    test_ratio = vapply(res, `[[`, numeric(1), "test_ratio"),
    train_prop = vapply(res, `[[`, numeric(1), "train_prop"),
    test_prop = vapply(res, `[[`, numeric(1), "test_prop"),
    rule = vapply(res, function(r) if (is.null(r$rule)) NA_character_ else format(r$rule),
                  character(1)),
    stringsAsFactors = FALSE)
}

## Proportion of successful resamples whose best subgroup uses each predictor,
## in decreasing order (ties in alphabetical order, as in the legacy code).
variable_frequency <- function(res, predictors) {
  ok <- vapply(res, function(r) identical(r$status, "ok"), logical(1))
  if (!any(ok)) return(setNames(numeric(length(predictors)), predictors)[0L])
  used <- unlist(lapply(res[ok], `[[`, "vars"))
  cnt <- table(factor(used, levels = sort(predictors)))
  freq <- as.numeric(cnt) / sum(ok)
  names(freq) <- names(cnt)
  freq[order(-freq)]
}

## The eta / theta decision rule. `overall_ratio` is the effect ratio
## r = exp(-theta) of the complete data (the hazard ratio for time-to-event
## endpoints).
siga_decision <- function(rtab, freq, overall_ratio, eta, theta) {
  ok <- rtab$status == "ok"
  mean_ratio <- if (any(ok)) mean(rtab$test_ratio[ok], na.rm = TRUE) else NA_real_
  if (is.nan(mean_ratio)) mean_ratio <- NA_real_
  threshold <- (1 - eta) * overall_ratio
  selected <- names(freq)[freq >= theta & freq > 0]
  reason <- if (!any(ok) || is.na(mean_ratio)) {
    "no resample produced an evaluable subgroup"
  } else if (mean_ratio >= 1) {
    sprintf("the mean test-set effect ratio (%.3f) is not below 1", mean_ratio)
  } else if (mean_ratio > threshold) {
    sprintf("the mean test-set effect ratio (%.3f) exceeds (1 - eta) x overall ratio (%.3f)",
            mean_ratio, threshold)
  } else if (!length(selected)) {
    sprintf("no predictor was selected in at least %.0f%% of the resamples", 100 * theta)
  } else NULL
  list(subgroup = is.null(reason), reason = reason, eta = eta, theta = theta,
       mean_ratio = mean_ratio, overall_ratio = overall_ratio, threshold = threshold,
       improvement = 1 - mean_ratio / overall_ratio, selected = selected,
       n_ok = sum(ok), n_failed = sum(!ok))
}
