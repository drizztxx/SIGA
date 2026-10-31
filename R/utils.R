## Internal helpers --------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

## Sample quantile of type 7 from pre-sorted, NA-free values. The arithmetic
## is that of stats::quantile(x, probs, type = 7), so that the thresholds match
## the legacy SIGA code exactly, but sorting is done only once per variable.
q7 <- function(sorted, probs) {
  n <- length(sorted)
  index <- 1 + max(n - 1, 0) * probs
  lo <- floor(index)
  hi <- ceiling(index)
  qs <- sorted[lo]
  i <- which(index > lo & sorted[hi] != qs)
  h <- (index - lo)[i]
  qs[i] <- (1 - h) * qs[i] + h * sorted[hi[i]]
  qs
}

## Cache of objective values: membership keys and values in two vectors (looked
## up with match()). The keys must not be used as names in an environment:
## every name passed to assign() or env[[name]] becomes an R symbol, symbols
## are never garbage collected, and long simulations ran out of memory.
new_objective_cache <- function() {
  cache <- new.env(parent = emptyenv())
  cache$keys <- character()
  cache$vals <- numeric()
  cache
}

## Compact, collision-free keys for the columns of a logical matrix
## (one subgroup membership vector per column). Used to cache objective values.
mask_keys <- function(M) {
  n <- nrow(M)
  pad <- (8L - n %% 8L) %% 8L
  if (pad) M <- rbind(M, matrix(FALSE, pad, ncol(M)))
  nb <- (n + pad) %/% 8L
  hex <- as.character(packBits(as.vector(M)))
  apply(matrix(hex, nrow = nb), 2L, paste, collapse = "")
}

## Readable numbers for subgroup rules: integers are printed in full, other
## values with `digits` significant digits in fixed notation.
fmt_num <- function(x, digits = 4) {
  vapply(x, function(v) {
    if (!is.finite(v)) return(as.character(v))
    if (is.null(digits) || (v == round(v) && abs(v) < 1e15))
      return(format(v, scientific = FALSE, trim = TRUE, digits = 15))
    formatC(signif(v, digits), digits = digits, format = "fg")
  }, character(1))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 1e-4, formatC(p, digits = 2, format = "e"),
                                formatC(p, digits = 4, format = "f")))
}

check_number <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lower || x > upper ||
      (integer && x != round(x)))
    stop(sprintf("'%s' must be a single %s in [%s, %s].", name,
                 if (integer) "integer" else "number", lower, upper), call. = FALSE)
  invisible(x)
}

check_cores <- function(n_cores) {
  check_number(n_cores, "n_cores", lower = 1, integer = TRUE)
  avail <- parallel::detectCores()
  if (!is.na(avail) && n_cores > avail)
    stop(sprintf("'n_cores' (%d) exceeds the number of available cores (%d).",
                 n_cores, avail), call. = FALSE)
  as.integer(n_cores)
}

vmessage <- function(verbose, ...) if (isTRUE(verbose)) message(...)
