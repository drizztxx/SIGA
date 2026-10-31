## Candidate biomarker prescreening -------------------------------------------

## Number of candidates recommended in the SIGA paper (section 2.2).
default_n_candidates <- function(p) if (p <= 20) as.integer(ceiling(p / 2)) else 10L

## Design columns of one predictor for the screening models: numeric values,
## level codes for ordinal predictors, dummy variables (first level as
## reference) for binary and categorical predictors.
screen_design <- function(x, type) {
  switch(type,
    numeric = matrix(as.numeric(x), ncol = 1L),
    ordinal = matrix(as.numeric(as.integer(x)), ncol = 1L),
    {
      lev <- levels(x)
      d <- vapply(lev[-1L], function(l) as.numeric(x == l), numeric(length(x)))
      matrix(d, nrow = length(x))
    })
}

## Likelihood-ratio test of the treatment-by-predictor interaction for every
## predictor (F test for continuous endpoints), on the complete cases of that
## predictor. Failed fits get p = 1, as in the original implementation.
screen_lrt <- function(ep, X, vars, types) {
  rows <- lapply(vars, function(v) {
    xd <- screen_design(X[[v]], types[[v]])
    ok <- which(stats::complete.cases(xd))
    res <- data.frame(variable = v, type = types[[v]], statistic = NA_real_,
                      df = ncol(xd), p_value = 1, stringsAsFactors = FALSE)
    if (length(ok) < 10L) return(res)
    epv <- ep_subset(ep, ok)
    xd <- xd[ok, , drop = FALSE]
    zx <- epv$z * xd
    f0 <- model_fit_stat(epv, cbind(epv$z, xd))
    f1 <- model_fit_stat(epv, cbind(epv$z, xd, zx))
    if (is.null(f0) || is.null(f1)) return(res)
    df <- ncol(xd)
    if (ep$type == "continuous") {
      df2 <- f1$n - f1$df
      if (df2 < 1L || f1$rss <= 0) return(res)
      stat <- ((f0$rss - f1$rss) / df) / (f1$rss / df2)
      p <- stats::pf(stat, df, df2, lower.tail = FALSE)
    } else {
      stat <- 2 * (f1$loglik - f0$loglik)
      p <- stats::pchisq(stat, df, lower.tail = FALSE)
    }
    res$statistic <- stat
    res$p_value <- if (is.na(p)) 1 else p
    res
  })
  do.call(rbind, rows)
}

## Random forest screening (legacy method 1): in each arm grow a random
## forest and compute the permutation importance; rank predictors by the
## absolute difference of importance between the arms.
screen_rf <- function(ep, X, vars, types, ntree = 2000, ...) {
  if (!requireNamespace("randomForestSRC", quietly = TRUE))
    stop("screen = \"rf\" requires the 'randomForestSRC' package.", call. = FALSE)
  old <- options(rf.cores = 1L, mc.cores = 1L)
  on.exit(options(old), add = TRUE)
  d <- X[vars]
  if (ep$type == "tte") {
    d$.time <- ep$time
    d$.status <- ep$status
    form <- stats::as.formula("Surv(.time, .status) ~ .")
  } else {
    d$.y <- if (ep$type == "binary") factor(ep$y) else ep$y
    form <- stats::as.formula(".y ~ .")
  }
  vim <- function(arm) {
    fit <- randomForestSRC::rfsrc(form, data = d[ep$z == arm, , drop = FALSE], ntree = ntree,
                                  importance = TRUE, na.action = "na.omit", block.size = 1,
                                  seed = -sample.int(1e6, 1L), ...)
    imp <- fit$importance
    if (is.matrix(imp)) imp <- imp[, 1L]
    imp[vars]
  }
  score <- vim(1) - vim(0)
  score[is.na(score)] <- 0
  data.frame(variable = vars, type = unname(types[vars]), statistic = unname(score),
             df = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
}

screen_predictors <- function(prep, vars = names(prep$types), method = "lrt",
                              n_candidates = NULL, screen_args = list()) {
  p <- length(vars)
  k <- min(n_candidates %||% default_n_candidates(p), p)
  tab <- switch(method,
    lrt = screen_lrt(prep$ep, prep$X, vars, prep$types),
    rf = do.call(screen_rf, c(list(prep$ep, prep$X, vars, prep$types), screen_args)),
    none = data.frame(variable = vars, type = unname(prep$types[vars]), statistic = NA_real_,
                      df = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
  ord <- switch(method, lrt = order(tab$p_value), rf = order(-abs(tab$statistic)),
                none = seq_len(p))
  tab <- tab[ord, , drop = FALSE]
  rownames(tab) <- NULL
  tab$rank <- seq_len(p)
  if (method == "none") k <- p
  tab$selected <- tab$rank <= k
  list(table = tab, selected = tab$variable[tab$selected], method = method, n_candidates = k)
}

#' Prescreening of candidate biomarkers
#'
#' Ranks the predictors by their potential to interact with the treatment and
#' selects the top candidates for the genetic algorithm search (SIGA paper,
#' section 2.2). With `method = "lrt"` (default) a model with the treatment,
#' the predictor and their interaction is compared with the model without the
#' interaction by a likelihood-ratio test (Cox regression for time-to-event,
#' logistic regression for binary endpoints) or an F test (linear regression
#' for continuous endpoints). With `method = "rf"` random forests (package
#' \pkg{randomForestSRC}) are grown in each arm and predictors are ranked by the
#' absolute difference of their permutation importance between the arms.
#'
#' @inheritParams siga_search
#' @param method `"lrt"`, `"rf"` or `"none"` (keep all predictors).
#' @param n_candidates number of candidates \eqn{k} to keep. The default
#'   follows the recommendation of the SIGA paper: \eqn{k = \lceil p/2 \rceil}
#'   for \eqn{p \le 20} predictors and \eqn{k = 10} otherwise.
#' @param ... further arguments passed to [randomForestSRC::rfsrc()] for
#'   `method = "rf"` (e.g. `ntree`, default 2000).
#' @return An object of class `siga_screen` with the ranking `table` and the
#'   `selected` predictors.
#' @seealso [siga()]
#' @examples
#' dat <- siga_simulate(300, model = "B3", seed = 5)
#' siga_screen(Surv(time, status) ~ trt | ., data = dat)
#' @export
siga_screen <- function(formula, data, method = c("lrt", "rf", "none"), n_candidates = NULL,
                        control_arm = NULL, endpoint = NULL, direction = NULL,
                        predictor_types = NULL, seed = NULL, ...) {
  method <- match.arg(method)
  if (!is.null(n_candidates)) check_number(n_candidates, "n_candidates", lower = 1, integer = TRUE)
  prep <- prepare_siga_data(formula, data, control_arm, endpoint, direction, predictor_types)
  seed <- resolve_seed(seed)
  res <- with_rng_state(make_streams(1L, seed)[[1L]],
                        screen_predictors(prep, method = method, n_candidates = n_candidates,
                                          screen_args = list(...)))
  res$info <- prep$info
  res$call <- match.call()
  class(res) <- "siga_screen"
  res
}

#' @export
print.siga_screen <- function(x, ...) {
  cat(sprintf("Candidate prescreening (%s), %s endpoint: %d of %d predictors selected\n",
              x$method, x$info$endpoint, length(x$selected), nrow(x$table)))
  tab <- x$table
  tab$statistic <- formatC(tab$statistic, digits = 3, format = "fg")
  tab$p_value <- fmt_p(tab$p_value)
  print(tab, row.names = FALSE)
  invisible(x)
}
