## Fast model fitting ------------------------------------------------------
##
## The genetic algorithm evaluates thousands of candidate subgroups, so the
## models are fitted with the low-level fitters of the survival and stats
## packages instead of the formula interface. Results are identical to
## coxph(), glm() and lm().

cox_fit <- function(x, y, ctl) {
  storage.mode(x) <- "double"
  tryCatch(suppressWarnings(
    survival::coxph.fit(x, y, strata = NULL, offset = NULL, init = NULL, control = ctl,
                        weights = NULL, method = "efron", rownames = NULL, resid = FALSE)),
    error = function(e) NULL)
}

## Estimate, standard error and residual degrees of freedom (Inf for
## z-tests) of coefficient `j` of the design `x` (without intercept).
## Returns NULL when the coefficient is not estimable.
coef_stat <- function(ep, x, j) {
  if (ep$type == "tte") {
    fit <- cox_fit(x, ep$surv, ep$cox_ctl)
    if (is.null(fit)) return(NULL)
    b <- unname(fit$coefficients[j])
    v <- fit$var[j, j]
    if (!is.finite(b) || !is.finite(v) || v <= 0) return(NULL)
    return(c(b, sqrt(v), Inf))
  }
  x1 <- cbind(1, x)
  storage.mode(x1) <- "double"
  j1 <- j + 1L
  p <- ncol(x1)
  if (ep$type == "binary") {
    fit <- tryCatch(suppressWarnings(stats::glm.fit(x1, ep$y, family = ep$family)),
                    error = function(e) NULL)
    if (is.null(fit) || fit$rank < p) return(NULL)
    covu <- chol2inv(fit$qr$qr[seq_len(p), seq_len(p), drop = FALSE])
    pos <- match(j1, fit$qr$pivot)
    b <- fit$coefficients[j1]
    v <- covu[pos, pos]
    df <- Inf
  } else {
    fit <- tryCatch(stats::.lm.fit(x1, ep$y), error = function(e) NULL)
    if (is.null(fit) || fit$rank < p || any(fit$pivot != seq_len(p))) return(NULL)
    df <- nrow(x1) - p
    if (df < 1L) return(NULL)
    s2 <- sum(fit$residuals^2) / df
    covu <- chol2inv(fit$qr[seq_len(p), seq_len(p), drop = FALSE])
    b <- fit$coefficients[j1]
    v <- s2 * covu[j1, j1]
  }
  if (!is.finite(b) || !is.finite(v) || v <= 0) return(NULL)
  c(unname(b), sqrt(v), df)
}

## One-sided p-value for a statistic oriented so that large values indicate
## benefit.
p_upper <- function(stat, df = Inf) {
  p <- if (is.finite(df)) stats::pt(stat, df, lower.tail = FALSE) else
    stats::pnorm(stat, lower.tail = FALSE)
  if (is.na(p)) 1 else p
}

## Log-likelihood (or minus half the deviance / residual sum of squares for
## the F test) used by the screening likelihood-ratio tests.
model_fit_stat <- function(ep, x) {
  if (ep$type == "tte") {
    fit <- cox_fit(x, ep$surv, ep$cox_ctl)
    if (is.null(fit) || any(!is.finite(fit$coefficients))) return(NULL)
    return(list(loglik = fit$loglik[2L], df = ncol(x)))
  }
  x1 <- cbind(1, x)
  storage.mode(x1) <- "double"
  if (ep$type == "binary") {
    fit <- tryCatch(suppressWarnings(stats::glm.fit(x1, ep$y, family = ep$family)),
                    error = function(e) NULL)
    if (is.null(fit) || fit$rank < ncol(x1)) return(NULL)
    return(list(loglik = -fit$deviance / 2, df = ncol(x1)))
  }
  fit <- tryCatch(stats::.lm.fit(x1, ep$y), error = function(e) NULL)
  if (is.null(fit) || fit$rank < ncol(x1)) return(NULL)
  list(rss = sum(fit$residuals^2), df = ncol(x1), n = nrow(x1))
}

## Treatment effect within the observations `mask`, on an oriented log-ratio
## scale theta (theta > 0 means benefit):
##   tte        theta = -log(HR)               (direction "lower")
##   binary     theta =  log(OR)               (direction "higher"; Haldane
##                                               correction if a cell is empty)
##   continuous theta = (pi / sqrt(3)) * d      d = mean difference / pooled SD
## The ratio used by the eta rule is r = exp(-theta) (lower is better), i.e.
## the hazard ratio for survival endpoints as in the SIGA paper.
effect_theta <- function(ep, mask) {
  z <- ep$z[mask]
  n1 <- sum(z == 1)
  n0 <- length(z) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  est <- switch(ep$type,
    tte = {
      fit <- cox_fit(matrix(z, ncol = 1L), ep$surv[mask], ep$cox_ctl)
      if (is.null(fit)) NA_real_ else unname(fit$coefficients[1L])
    },
    binary = {
      y <- ep$y[mask]
      a <- sum(y[z == 1]); b <- n1 - a
      c0 <- sum(y[z == 0]); d <- n0 - c0
      if (min(a, b, c0, d) == 0) {
        a <- a + 0.5; b <- b + 0.5; c0 <- c0 + 0.5; d <- d + 0.5
      }
      log((a * d) / (b * c0))
    },
    continuous = {
      y <- ep$y[mask]
      (pi / sqrt(3)) * (mean(y[z == 1]) - mean(y[z == 0])) / ep$sd_y
    })
  if (!is.finite(est)) return(NA_real_)
  ep$sign * est
}

## Ratio r = exp(-theta) of a subgroup in a test set, with the conventions of
## the legacy implementation: r = 1 if the subgroup or its complement is empty
## and NA if the effect is not estimable or extreme (r > 100 or r < 0.001).
test_ratio <- function(ep, mask) {
  if (!any(mask) || all(mask)) return(1)
  r <- exp(-effect_theta(ep, mask))
  if (!is.finite(r) || r > 100 || r < 0.001) NA_real_ else r
}
