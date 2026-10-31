## Treatment effects in a subgroup, its complement and the interaction ------

effect_measure <- function(type) {
  switch(type, tte = "hazard ratio", binary = "odds ratio", continuous = "mean difference")
}

km_median <- function(time, status) {
  if (!length(time)) return(NA_real_)
  fit <- survival::survfit(survival::Surv(time, status) ~ 1)
  unname(summary(fit)$table["median"])
}

## Estimate (HR, OR or mean difference, treatment vs control), 95% CI,
## two-sided p-value and one-sided p-value in the direction of benefit from a
## coefficient statistic c(estimate, se, df).
stat_ci <- function(s, ep) {
  if (is.null(s)) return(c(estimate = NA, lower = NA, upper = NA, p_value = NA, p_benefit = NA))
  q <- if (is.finite(s[3L])) stats::qt(0.975, s[3L]) else stats::qnorm(0.975)
  ci <- s[1L] + c(0, -1, 1) * q * s[2L]
  if (ep$type != "continuous") ci <- exp(ci)
  stat <- s[1L] / s[2L]
  p2 <- if (is.finite(s[3L])) 2 * stats::pt(-abs(stat), s[3L]) else 2 * stats::pnorm(-abs(stat))
  c(estimate = ci[1L], lower = ci[2L], upper = ci[3L], p_value = p2,
    p_benefit = p_upper(ep$sign * stat, s[3L]))
}

group_effect <- function(ep, mask) {
  idx <- which(mask)
  z <- ep$z[idx]
  n1 <- sum(z == 1)
  n0 <- length(z) - n1
  arm <- function(f, sel) if (any(sel)) f(sel) else NA_real_
  if (ep$type == "tte") {
    st <- ep$status[idx]; tm <- ep$time[idx]
    ev <- c(sum(st[z == 1]), sum(st[z == 0]))
    smry <- c(arm(function(s) km_median(tm[s], st[s]), z == 1),
              arm(function(s) km_median(tm[s], st[s]), z == 0))
  } else {
    y <- ep$y[idx]
    ev <- if (ep$type == "binary") c(sum(y[z == 1]), sum(y[z == 0])) else c(NA_real_, NA_real_)
    smry <- c(arm(function(s) mean(y[s]), z == 1), arm(function(s) mean(y[s]), z == 0))
  }
  s <- if (n1 > 0 && n0 > 0) coef_stat(ep_subset(ep, idx), matrix(z, ncol = 1L), 1L) else NULL
  est <- stat_ci(s, ep)
  data.frame(n = length(idx), n_trt = n1, n_ctrl = n0, events_trt = ev[1L], events_ctrl = ev[2L],
             summary_trt = smry[1L], summary_ctrl = smry[2L], estimate = est[["estimate"]],
             lower = est[["lower"]], upper = est[["upper"]], p_value = est[["p_value"]],
             p_benefit = est[["p_benefit"]])
}

## Effects in the overall population, the subgroup and its complement, and
## the treatment-by-subgroup interaction (ratio of hazard ratios / odds ratios,
## or difference of mean differences).
effect_summary <- function(ep, mask) {
  groups <- list(overall = rep(TRUE, ep$n), subgroup = mask, complement = !mask)
  tab <- do.call(rbind, lapply(groups, group_effect, ep = ep))
  tab <- cbind(group = names(groups), tab)
  rownames(tab) <- NULL
  I <- as.numeric(mask)
  s <- if (any(mask) && !all(mask)) coef_stat(ep, cbind(ep$z, I, ep$z * I), 3L) else NULL
  inter <- stat_ci(s, ep)
  structure(list(table = tab, interaction = as.list(inter), endpoint = ep$type,
                 direction = ep$direction, measure = effect_measure(ep$type)),
            class = "siga_effect")
}

#' Treatment effects in a subgroup and its complement
#'
#' Estimates the treatment effect (hazard ratio, odds ratio or mean difference
#' of the treated versus the control arm) in the overall population, in a
#' subgroup and in its complement, together with the treatment-by-subgroup
#' interaction. Use it to evaluate a subgroup identified by [siga()] or
#' [siga_search()] on the training data or on independent validation data.
#'
#' Hazard ratios and odds ratios come from Cox and logistic regression with the
#' treatment indicator as the only covariate, mean differences from linear
#' regression; confidence intervals are Wald (z or t) intervals. For time-to-event
#' endpoints `summary_trt` and `summary_ctrl` are Kaplan-Meier median survival
#' times, for binary endpoints event rates and for continuous endpoints means.
#' The interaction is the ratio of hazard ratios (odds ratios) or the difference
#' of mean differences between subgroup and complement; `p_benefit` is the
#' one-sided p-value in the direction of a larger benefit in the subgroup.
#'
#' @param object a `siga` or `siga_search` fit, or a [siga_rule].
#' @param data data frame with the outcome, treatment and predictors. By default
#'   the data used to fit `object`.
#' @param formula,control_arm,endpoint,direction only needed when `object` is a
#'   rule; see [siga()].
#' @return An object of class `siga_effect` with components `table` (one row per
#'   group) and `interaction`.
#' @seealso [siga()], [predict.siga()]
#' @examples
#' dat <- siga_simulate(300, model = "B1", seed = 3)
#' fit <- siga_search(Surv(time, status) ~ trt | x1 + x2 + x4, data = dat,
#'                    control = siga_control(pop_size = 20, max_gen = 15, n_subpop = 1),
#'                    seed = 1)
#' subgroup_effect(fit)
#' ## validation in independent data
#' new <- siga_simulate(1000, model = "B1", seed = 4)
#' subgroup_effect(fit, data = new)
#' @export
subgroup_effect <- function(object, data = NULL, formula = NULL, control_arm = NULL,
                            endpoint = NULL, direction = NULL) {
  if (inherits(object, c("siga", "siga_search"))) {
    rule <- object$rule
    if (is.null(rule)) stop("No subgroup was identified.", call. = FALSE)
    if (is.null(data)) {
      ep <- object$data$ep
      mask <- rule_members(rule, object$data$X)
    } else {
      info <- object$data$info
      pe <- prepare_endpoint(object$formula, data, info$control_arm, info$endpoint,
                             info$direction, quiet = TRUE)
      ep <- pe$ep
      mask <- rule_members(rule, as.data.frame(data)[pe$rows, , drop = FALSE])
    }
  } else if (inherits(object, "siga_rule")) {
    if (is.null(formula) || is.null(data))
      stop("'formula' and 'data' are required when 'object' is a rule.", call. = FALSE)
    rule <- object
    pe <- prepare_endpoint(formula, data, control_arm, endpoint, direction, quiet = TRUE)
    ep <- pe$ep
    mask <- rule_members(rule, as.data.frame(data)[pe$rows, , drop = FALSE])
  } else {
    stop("'object' must be a 'siga' or 'siga_search' fit or a 'siga_rule'.", call. = FALSE)
  }
  res <- effect_summary(ep, mask)
  res$rule <- format(rule)
  res
}

#' @export
print.siga_effect <- function(x, digits = 3, ...) {
  tab <- x$table
  lab <- switch(x$endpoint, tte = c("events", "median"), binary = c("events", "rate"),
                continuous = c(NA, "mean"))
  if (!is.null(x$rule)) cat("Subgroup:", x$rule, "\n")
  cat(sprintf("Treatment effect (%s, treated vs control; benefit = %s):\n", x$measure,
              if (x$direction == "lower") "lower" else "higher"))
  f <- function(v) formatC(v, digits = digits, format = "fg")
  out <- data.frame(group = tab$group, n = tab$n,
                    `n trt/ctrl` = paste0(tab$n_trt, "/", tab$n_ctrl), check.names = FALSE)
  if (!is.na(lab[1L])) out[[paste(lab[1L], "trt/ctrl")]] <-
    paste0(tab$events_trt, "/", tab$events_ctrl)
  out[[paste(lab[2L], "trt/ctrl")]] <- paste0(f(tab$summary_trt), "/", f(tab$summary_ctrl))
  out$estimate <- f(tab$estimate)
  out$`95% CI` <- paste0("(", f(tab$lower), ", ", f(tab$upper), ")")
  out$p <- fmt_p(tab$p_value)
  print(out, row.names = FALSE, right = FALSE)
  it <- x$interaction
  cat(sprintf("Interaction (%s): %s (95%% CI %s, %s), one-sided p = %s\n",
              if (x$endpoint == "continuous") "difference of mean differences" else
                paste("ratio of", paste0(x$measure, "s")),
              f(it$estimate), f(it$lower), f(it$upper), fmt_p(it$p_benefit)))
  invisible(x)
}
