## Subgroup rules ------------------------------------------------------------
##
## A rule is a list of conditions combined with AND. Every condition has
##   var    predictor name
##   type   "numeric", "ordinal", "binary" or "categorical"
##   op     "<=", ">=" (numeric / ordinal), "==" (binary) or "in" (categorical)
##   value  threshold (numeric; level code for ordinal) or level(s)
##   levels level labels (ordinal, binary, categorical)

new_siga_rule <- function(conditions = list()) {
  structure(list(conditions = conditions), class = "siga_rule")
}

#' Subgroup rules
#'
#' Subgroups found by SIGA are stored as `siga_rule` objects: a set of
#' conditions on predictors combined with AND. These functions print a rule,
#' convert it to text or to an R expression, list its variables and evaluate
#' subgroup membership on new data.
#'
#' Numeric thresholds are stored at full precision and printed with `digits`
#' significant digits. Patients with a missing value of a predictor used in the
#' rule are not members of the subgroup.
#'
#' @param x a `siga_rule` object.
#' @param digits number of significant digits for numeric thresholds; `NULL`
#'   prints full precision.
#' @param style `"text"` for a readable rule (e.g. `x4 >= 1.02 & x5 = 1`) or
#'   `"R"` for an R expression that can be used with [subset()].
#' @param data a data frame containing the predictors used in the rule.
#' @param ... unused.
#' @return `format()` and `as.character()` return a character string,
#'   `rule_variables()` a character vector and `rule_members()` a logical
#'   vector with one element per row of `data`.
#' @name siga_rule
#' @examples
#' dat <- siga_simulate(200, model = "B1", seed = 1)
#' fit <- siga_search(Surv(time, status) ~ trt | x1 + x4, data = dat,
#'                    control = siga_control(pop_size = 20, max_gen = 10, n_subpop = 1),
#'                    seed = 1)
#' fit$rule
#' format(fit$rule, style = "R")
#' table(rule_members(fit$rule, dat))
NULL

#' @rdname siga_rule
#' @export
format.siga_rule <- function(x, digits = 4, style = c("text", "R"), ...) {
  style <- match.arg(style)
  conds <- x$conditions
  if (!length(conds)) return(if (style == "R") "TRUE" else "<all patients>")
  q <- function(v) if (style == "R") paste0("\"", gsub("\"", "\\\\\"", v), "\"") else v
  vn <- function(v) if (style == "R" && make.names(v) != v) paste0("`", v, "`") else v
  parts <- vapply(conds, function(cd) {
    switch(cd$type,
      numeric = paste(vn(cd$var), cd$op, fmt_num(cd$value, digits)),
      ordinal = {
        code <- if (cd$op == "<=") floor(cd$value + 1e-8) else ceiling(cd$value - 1e-8)
        code <- min(max(code, 1L), length(cd$levels))
        paste(vn(cd$var), cd$op, q(cd$levels[code]))
      },
      binary = if (style == "R") paste(vn(cd$var), "==", q(cd$value)) else
        paste(cd$var, "=", cd$value),
      categorical = if (style == "R") {
        paste0(vn(cd$var), " %in% c(", paste(q(cd$value), collapse = ", "), ")")
      } else {
        paste0(cd$var, " in {", paste(cd$value, collapse = ", "), "}")
      })
  }, character(1))
  paste(parts, collapse = " & ")
}

#' @rdname siga_rule
#' @export
as.character.siga_rule <- function(x, ...) format(x, ...)

#' @rdname siga_rule
#' @export
print.siga_rule <- function(x, digits = 4, ...) {
  cat(format(x, digits = digits), "\n")
  invisible(x)
}

#' @rdname siga_rule
#' @export
rule_variables <- function(x) {
  if (is.null(x)) return(character())
  unique(vapply(x$conditions, `[[`, character(1), "var"))
}

## Membership for one condition.
condition_members <- function(cd, data) {
  x <- data[[cd$var]]
  if (is.null(x)) stop(sprintf("Variable '%s' not found in 'data'.", cd$var), call. = FALSE)
  hit <- switch(cd$type,
    numeric = {
      if (is.factor(x)) x <- if (is.ordered(x)) as.integer(x) else as.character(x)
      x <- suppressWarnings(as.numeric(x))
      if (cd$op == "<=") x <= cd$value else x >= cd$value
    },
    ordinal = {
      code <- match(as.character(x), cd$levels)
      if (cd$op == "<=") code <= cd$value else code >= cd$value
    },
    binary = , categorical = as.character(x) %in% as.character(cd$value) & !is.na(x))
  hit & !is.na(hit)
}

#' @rdname siga_rule
#' @export
rule_members <- function(x, data) {
  if (!inherits(x, "siga_rule")) stop("'x' must be a 'siga_rule' object.", call. = FALSE)
  m <- rep(TRUE, nrow(data))
  for (cd in x$conditions) m <- m & condition_members(cd, data)
  m
}

## Replace thresholds by the boundary values observed in the subgroup (legacy
## SIGA convention: e.g. "x <= 60" instead of the interpolated quantile 62.5)
## and level sets by the levels present in the subgroup. Membership on the
## training data is unchanged.
canonicalize_conditions <- function(conds, X, mask) {
  lapply(conds, function(cd) {
    x <- X[[cd$var]][mask]
    if (cd$type == "numeric" || cd$type == "ordinal") {
      xv <- if (cd$type == "ordinal") as.integer(x) else as.numeric(x)
      xv <- xv[!is.na(xv)]
      if (length(xv)) cd$value <- if (cd$op == ">=") min(xv) else max(xv)
    } else if (cd$type == "categorical") {
      present <- intersect(cd$levels, unique(as.character(x)))
      cd$value <- intersect(cd$value, present)
    }
    cd
  })
}

## Smallest subset of conditions (fewest conditions first, in combn order)
## reproducing the membership `target` on X.
simplify_conditions <- function(conds, X, target) {
  k <- length(conds)
  if (k <= 1L) return(conds)
  masks <- lapply(conds, condition_members, data = X)
  for (size in seq_len(k - 1L)) {
    for (cb in utils::combn(k, size, simplify = FALSE)) {
      if (identical(Reduce(`&`, masks[cb]), target)) return(conds[cb])
    }
  }
  conds
}

## Canonical rule for a chromosome: decode, canonicalise thresholds, drop
## redundant conditions, order conditions by variable name.
chrom_rule <- function(chrom, enc, X) {
  conds <- decode_conditions(chrom, enc)
  rule0 <- new_siga_rule(conds)
  mask <- rule_members(rule0, X)
  canon <- canonicalize_conditions(conds, X, mask)
  if (!identical(rule_members(new_siga_rule(canon), X), mask)) canon <- conds
  simp <- simplify_conditions(canon, X, mask)
  if (length(simp)) simp <- simp[order(vapply(simp, `[[`, character(1), "var"))]
  list(rule = new_siga_rule(simp), mask = mask)
}
