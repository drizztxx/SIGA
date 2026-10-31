## Methods for siga and siga_search objects -----------------------------------

endpoint_label <- function(info) {
  switch(info$endpoint, tte = "time-to-event", binary = "binary", continuous = "continuous")
}

ratio_label <- function(type) {
  switch(type, tte = "hazard ratio", binary = "inverse odds ratio",
         continuous = "effect ratio exp(-pi d / sqrt(3))")
}

print_effect_line <- function(eff, digits = 3) {
  if (is.null(eff)) return(invisible())
  tab <- eff$table
  f <- function(v) formatC(v, digits = digits, format = "fg")
  sub <- tab[tab$group == "subgroup", ]
  comp <- tab[tab$group == "complement", ]
  cat(sprintf("  %s in subgroup %s (95%% CI %s, %s), in complement %s; interaction one-sided p = %s\n",
              eff$measure, f(sub$estimate), f(sub$lower), f(sub$upper), f(comp$estimate),
              fmt_p(eff$interaction$p_benefit)))
}

#' @export
print.siga <- function(x, digits = 3, ...) {
  info <- x$data$info
  ep <- x$data$ep
  d <- x$decision
  f <- function(v) formatC(v, digits = digits, format = "fg")
  cat("SIGA: subgroup identification based on a genetic algorithm\n")
  cat(sprintf("Endpoint: %s, %s; %d patients (%d %s = %s, %d control = %s)\n",
              endpoint_label(info), info$outcome, ep$n, sum(ep$z == 1), info$treatment,
              info$treatment_arm, sum(ep$z == 0), info$control_arm))
  cat(sprintf("Overall %s: %s\n", ratio_label(info$endpoint), f(x$overall_ratio)))
  cat(sprintf("Resampling: %d subsamples (%d evaluable); mean test-set ratio %s vs threshold (1 - %.2f) x %s = %s\n",
              x$settings$n_resample, d$n_ok, f(d$mean_ratio), d$eta, f(d$overall_ratio),
              f(d$threshold)))
  fr <- x$frequency[x$frequency > 0]
  if (length(fr)) {
    top <- utils::head(fr, 6L)
    cat("Selection frequency:", paste0(names(top), " (", formatC(top, digits = 2, format = "f"), ")",
                                       collapse = ", "),
        if (length(fr) > 6L) "..." else "", "\n")
  }
  if (x$subgroup) {
    cat(sprintf("Predictors selected (theta = %.2f): %s\n", d$theta, paste(d$selected, collapse = ", ")))
    sub <- x$effects$table[x$effects$table$group == "subgroup", ]
    cat(sprintf("Subgroup: %s  [%d patients, %.1f%%]\n", format(x$rule), sub$n, 100 * sub$n / ep$n))
    print_effect_line(x$effects, digits)
  } else {
    cat("No subgroup identified:", d$reason, "\n")
  }
  invisible(x)
}

#' Summarise a SIGA fit
#'
#' @param object a `siga` fit.
#' @param ... unused.
#' @return An object of class `summary.siga` with the decision, the predictor
#'   frequencies, the prescreening of the complete data, a summary of the
#'   resamples, the candidates of the final sub-populations and the effect table.
#' @export
summary.siga <- function(object, ...) {
  r <- object$resamples
  ok <- r$status == "ok"
  structure(list(
    fit = object,
    frequency = data.frame(variable = names(object$frequency),
                           frequency = unname(object$frequency),
                           selected = names(object$frequency) %in% object$decision$selected),
    resampling = c(n = nrow(r), evaluable = sum(ok),
                   stats::quantile(r$test_ratio[ok], c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                   mean = object$decision$mean_ratio),
    screening = object$screening$table,
    candidates = if (!is.null(object$final)) object$final$candidates else NULL,
    effects = object$effects), class = "summary.siga")
}

#' @export
print.summary.siga <- function(x, digits = 3, ...) {
  print(x$fit, digits = digits)
  cat("\nTest-set effect ratios of the resampled subgroups:\n")
  print(round(x$resampling, digits))
  cat("\nPredictor selection frequencies:\n")
  fr <- x$frequency[x$frequency$frequency > 0, , drop = FALSE]
  fr$frequency <- round(fr$frequency, 3)
  print(fr, row.names = FALSE)
  cat(sprintf("\nPrescreening of the complete data (%s):\n", x$fit$settings$screen))
  sc <- x$screening
  sc$statistic <- signif(sc$statistic, digits)
  sc$p_value <- fmt_p(sc$p_value)
  print(utils::head(sc, 20L), row.names = FALSE)
  if (!is.null(x$candidates)) {
    cat("\nBest subgroup of every sub-population of the final search:\n")
    cand <- x$candidates
    cand$objective <- signif(cand$objective, digits)
    cand$proportion <- round(cand$proportion, 3)
    print(cand, row.names = FALSE)
  }
  if (!is.null(x$effects)) {
    cat("\n")
    print(x$effects, digits = digits)
  }
  invisible(x)
}

#' Subgroup membership of new patients
#'
#' @param object a `siga` or `siga_search` fit.
#' @param newdata a data frame with the predictors of the subgroup rule. By
#'   default the training data.
#' @param ... unused.
#' @return A logical vector: `TRUE` for patients in the identified subgroup.
#'   If no subgroup was identified all values are `FALSE` (with a warning).
#' @examples
#' dat <- siga_simulate(200, model = "B1", seed = 6)
#' fit <- siga_search(Surv(time, status) ~ trt | x1 + x4, data = dat,
#'                    control = siga_control(pop_size = 20, max_gen = 10, n_subpop = 1), seed = 1)
#' mean(predict(fit, siga_simulate(500, model = "B1", seed = 7)))
#' @export
predict.siga <- function(object, newdata = NULL, ...) {
  X <- if (is.null(newdata)) object$data$X else as.data.frame(newdata)
  if (is.null(object$rule)) {
    warning("No subgroup was identified; returning FALSE for all patients.", call. = FALSE)
    return(rep(FALSE, nrow(X)))
  }
  rule_members(object$rule, X)
}

#' @rdname predict.siga
#' @export
predict.siga_search <- predict.siga

#' Plot a SIGA fit
#'
#' @param x a `siga` or `siga_search` fit.
#' @param type `"frequency"` (selection frequency of the predictors during
#'   resampling, with \eqn{\theta}), `"resamples"` (held-out effect ratios of the
#'   resampled subgroups, with their mean, the overall ratio and the threshold),
#'   `"convergence"` (best objective value per generation of every
#'   sub-population of the final search) or `"effect"` (Kaplan-Meier curves,
#'   event rates or outcome distributions by arm in the subgroup and its
#'   complement).
#' @param ... further graphical arguments.
#' @return `x`, invisibly.
#' @export
plot.siga <- function(x, type = c("frequency", "resamples", "convergence", "effect"), ...) {
  type <- match.arg(type)
  switch(type,
    frequency = {
      fr <- utils::head(x$frequency[x$frequency > 0], 20L)
      if (!length(fr)) stop("No predictor was selected in any resample.", call. = FALSE)
      op <- par(mar = c(4, max(4, max(nchar(names(fr))) * 0.6), 2, 1))
      on.exit(par(op))
      barplot(rev(fr), horiz = TRUE, las = 1, xlim = c(0, 1), xlab = "Selection frequency",
              main = "Predictor selection in the resamples", ...)
      abline(v = x$decision$theta, lty = 2)
    },
    resamples = {
      r <- x$resamples$test_ratio[x$resamples$status == "ok"]
      r <- r[!is.na(r)]
      if (!length(r)) stop("No evaluable resample.", call. = FALSE)
      d <- x$decision
      hist(r, breaks = 20, main = "Held-out effect ratios of the resampled subgroups",
           xlab = ratio_label(x$data$info$endpoint),
           xlim = range(c(r, d$overall_ratio, d$threshold, 1)), ...)
      abline(v = c(d$mean_ratio, d$overall_ratio, d$threshold), lty = c(1, 2, 3), lwd = 2)
      legend("topright", c("mean", "overall", "threshold"), lty = 1:3, lwd = 2, bty = "n")
    },
    convergence = {
      if (is.null(x$final)) stop("No final search was run.", call. = FALSE)
      plot_convergence(x$final$history, ...)
    },
    effect = plot_effect(x, ...))
  invisible(x)
}

#' @rdname plot.siga
#' @export
plot.siga_search <- function(x, type = c("convergence", "effect"), ...) {
  type <- match.arg(type)
  if (type == "convergence") plot_convergence(x$history, ...) else plot_effect(x, ...)
  invisible(x)
}

plot_convergence <- function(history, ...) {
  if (is.null(history)) stop("No convergence history stored (keep_history = FALSE).", call. = FALSE)
  h <- log10(pmax(history, 1e-300))
  gen <- seq_len(nrow(h)) - 1L
  plot(range(gen), range(h), type = "n", xlab = "Generation",
       ylab = "log10(best objective value)", main = "Convergence of the sub-populations", ...)
  for (j in seq_len(ncol(h))) lines(gen, h[, j], col = j, lty = j)
  if (ncol(h) > 1L) legend("topright", colnames(h), col = seq_len(ncol(h)), lty = seq_len(ncol(h)),
                           bty = "n")
}

plot_effect <- function(x, ...) {
  if (is.null(x$rule)) stop("No subgroup was identified.", call. = FALSE)
  ep <- x$data$ep
  mask <- rule_members(x$rule, x$data$X)
  info <- x$data$info
  arms <- c(paste0(info$treatment, " = ", info$treatment_arm), paste0(info$treatment, " = ", info$control_arm))
  op <- par(mfrow = c(1, 2))
  on.exit(par(op))
  for (g in c("subgroup", "complement")) {
    m <- if (g == "subgroup") mask else !mask
    main <- if (g == "subgroup") paste("Subgroup:", format(x$rule)) else "Complement"
    if (ep$type == "tte") {
      fit <- survival::survfit(survival::Surv(ep$time[m], ep$status[m]) ~ factor(ep$z[m], levels = c(1, 0)))
      plot(fit, lty = 1:2, xlab = "Time", ylab = "Survival probability", main = main,
           cex.main = 0.8, ...)
      legend("bottomleft", arms, lty = 1:2, bty = "n")
    } else if (ep$type == "binary") {
      rates <- c(mean(ep$y[m & ep$z == 1]), mean(ep$y[m & ep$z == 0]))
      barplot(rates, names.arg = arms, ylim = c(0, 1), ylab = "Event rate", main = main,
              cex.main = 0.8, ...)
    } else {
      graphics::boxplot(ep$y[m] ~ factor(ep$z[m], levels = c(1, 0), labels = arms),
                        xlab = "", ylab = "Outcome", main = main, cex.main = 0.8, ...)
    }
  }
}

#' @export
print.siga_search <- function(x, digits = 3, ...) {
  info <- x$data$info
  cat(sprintf("SIGA genetic algorithm search (%s endpoint, %d patients, %d predictors, %d sub-population(s))\n",
              endpoint_label(info), x$data$ep$n, length(x$vars), nrow(x$candidates)))
  if (!x$valid) {
    cat("No valid subgroup found.\n")
    return(invisible(x))
  }
  sub <- x$effects$table[x$effects$table$group == "subgroup", ]
  cat(sprintf("Best subgroup: %s  [%d patients, %.1f%%]\n", format(x$rule), sub$n,
              100 * sub$n / x$data$ep$n))
  cat(sprintf("  objective (%s) = %s\n",
              if (is.function(x$control$objective)) "user-defined" else x$control$objective,
              formatC(x$objective, digits = digits, format = "g")))
  print_effect_line(x$effects, digits)
  if (nrow(x$candidates) > 1L) {
    cat("Best subgroup of every sub-population:\n")
    cand <- x$candidates
    cand$objective <- signif(cand$objective, digits)
    cand$proportion <- round(cand$proportion, 3)
    print(cand, row.names = FALSE)
  }
  invisible(x)
}
