#' Settings of the SIGA genetic algorithm
#'
#' Parameters of the genetic algorithm search used by [siga()],
#' [siga_search()] and [siga_cv()]. The defaults are the values recommended in
#' the SIGA paper (section 2.3.6).
#'
#' @param pop_size number of individuals (candidate subgroups) per
#'   sub-population, \eqn{m}. Default 50.
#' @param max_gen number of generations, \eqn{MaxG}. Default 50.
#' @param n_subpop number of isolated sub-populations, \eqn{s}, of the final
#'   search. Default 5.
#' @param ggap generation gap: fraction of the population replaced by
#'   offspring in every generation. Default 0.9.
#' @param min_size minimum subgroup size, \eqn{MinS}, as a proportion of the
#'   sample. As in the original implementation, the size of a subgroup is
#'   measured by \eqn{2 \min(n_1, n_0)}, where \eqn{n_1} and \eqn{n_0} are the
#'   numbers of treated and control patients in the subgroup, and subgroups
#'   smaller than `min_size` times the sample size receive the worst objective
#'   value. Default 0.1.
#' @param objective objective function to be minimised: `"interaction"`
#'   (default, the one-sided p-value of the treatment-by-subgroup interaction),
#'   `"subgroup"` (the one-sided p-value of the treatment effect in the
#'   subgroup, set to 1 when the complement has a smaller p-value) or
#'   `"differential"` (the one-sided p-value of the difference between the
#'   Wald statistics of the treatment effect in the subgroup and in its
#'   complement). A user-defined function `function(y, trt, subgroup, endpoint)`
#'   can be supplied instead; it receives the outcome (a `Surv` object for
#'   time-to-event endpoints, otherwise a numeric vector), the 0/1 treatment
#'   indicator, the logical subgroup membership and the endpoint type, and must
#'   return a single number (smaller is better).
#' @param mutation_rate probability that a gene is mutated; `NULL` (default)
#'   uses \eqn{1/k}, where \eqn{k} is the number of candidate predictors.
#' @param selection,recombination,mutation names of the \pkg{gatbxr} operators:
#'   selection (`"sus"` stochastic universal sampling or `"rws"` roulette wheel),
#'   recombination (`"reclin"` line, `"recint"` intermediate or `"recdis"`
#'   discrete recombination) and mutation (`"mutbga"` breeder genetic algorithm
#'   mutation).
#' @param keep_history keep the best objective value of every generation.
#' @return A list of class `siga_control`.
#' @seealso [siga()], [siga_search()]
#' @examples
#' siga_control()
#' siga_control(pop_size = 30, max_gen = 30, n_subpop = 3)
#' @export
siga_control <- function(pop_size = 50, max_gen = 50, n_subpop = 5, ggap = 0.9,
                         min_size = 0.1, objective = "interaction", mutation_rate = NULL,
                         selection = "sus", recombination = "reclin", mutation = "mutbga",
                         keep_history = TRUE) {
  check_number(pop_size, "pop_size", lower = 4, integer = TRUE)
  check_number(max_gen, "max_gen", lower = 1, integer = TRUE)
  check_number(n_subpop, "n_subpop", lower = 1, integer = TRUE)
  check_number(ggap, "ggap", lower = 0.05, upper = 1)
  check_number(min_size, "min_size", lower = 0, upper = 1)
  if (!is.function(objective)) {
    if (!is.character(objective) || length(objective) != 1L || !objective %in% objective_names)
      stop("'objective' must be one of ", paste(dQuote(objective_names, FALSE), collapse = ", "),
           " or a function.", call. = FALSE)
  }
  if (!is.null(mutation_rate)) check_number(mutation_rate, "mutation_rate", lower = 0, upper = 1)
  selection <- match.arg(selection, c("sus", "rws"))
  recombination <- match.arg(recombination, c("reclin", "recint", "recdis"))
  mutation <- match.arg(mutation, "mutbga")
  structure(list(pop_size = as.integer(pop_size), max_gen = as.integer(max_gen),
                 n_subpop = as.integer(n_subpop), ggap = ggap, min_size = min_size,
                 objective = objective, mutation_rate = mutation_rate,
                 selection = selection, recombination = recombination,
                 mutation = mutation, keep_history = isTRUE(keep_history)),
            class = "siga_control")
}

#' @export
print.siga_control <- function(x, ...) {
  obj <- if (is.function(x$objective)) "user-defined function" else x$objective
  cat("SIGA genetic algorithm settings\n")
  cat(sprintf("  population size %d, generations %d, sub-populations %d\n",
              x$pop_size, x$max_gen, x$n_subpop))
  cat(sprintf("  generation gap %.2f, minimum subgroup size %.0f%%\n", x$ggap, 100 * x$min_size))
  cat(sprintf("  objective: %s\n", obj))
  cat(sprintf("  operators: %s / %s / %s (mutation rate %s)\n", x$selection, x$recombination,
              x$mutation, if (is.null(x$mutation_rate)) "1/k" else format(x$mutation_rate)))
  invisible(x)
}
