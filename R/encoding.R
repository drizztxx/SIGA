## Genetic representation of subgroups --------------------------------------
##
## Every individual is a vector of genes B_j in [0, 4], one per candidate
## predictor X_j (SIGA paper, section 2.3.1):
##   B_j <= 2                 X_j is not used;
##   numeric / ordinal X_j:   2 < B_j <= 3  ->  X_j <= P(B_j - 2)
##                            3 < B_j <= 4  ->  X_j >= P(B_j - 3)
##                            (P(u) = type-7 sample quantile of X_j at u);
##   binary X_j:              2 < B_j <= 3  ->  X_j == first level
##                            3 < B_j <= 4  ->  X_j == second level;
##   categorical X_j (c levels): (2, 4] is cut into 2^c - 2 blocks, block s
##                            ->  X_j in C_s, the s-th non-empty proper subset
##                            of the levels (subsets ordered by size, then
##                            lexicographically by level order).
## The conditions of the used genes are combined with AND. Patients with a
## missing value of a used predictor are not in the subgroup.

max_categorical_levels <- 10L

## Pre-compute everything needed to decode genes of `vars` on the data `X`.
encode_predictors <- function(X, vars, types) {
  info <- lapply(vars, function(v) {
    x <- X[[v]]
    type <- types[[v]]
    switch(type,
      numeric = , ordinal = {
        xv <- if (type == "ordinal") as.numeric(as.integer(x)) else as.numeric(x)
        nas <- is.na(xv)
        list(var = v, type = type,
             levels = if (type == "ordinal") levels(x) else NULL,
             x_le = ifelse(nas, Inf, xv), x_ge = ifelse(nas, -Inf, xv),
             sorted = sort(xv[!nas]))
      },
      binary = {
        lev <- levels(x)
        list(var = v, type = type, levels = lev,
             is1 = !is.na(x) & x == lev[1L], is2 = !is.na(x) & x == lev[2L])
      },
      categorical = {
        lev <- levels(x)
        nl <- length(lev)
        if (nl > max_categorical_levels)
          stop(sprintf(paste0("Categorical predictor '%s' has %d levels; at most %d are supported ",
                              "(2^c - 2 subsets). Collapse levels or set predictor_types."),
                       v, nl, max_categorical_levels), call. = FALSE)
        combos <- unlist(lapply(seq_len(nl - 1L), function(s) utils::combn(nl, s, simplify = FALSE)),
                         recursive = FALSE)
        list(var = v, type = type, levels = lev,
             codes = ifelse(is.na(x), 0L, as.integer(x)),
             combos = combos, ncomb = length(combos))
      },
      stop("Unknown predictor type: ", type))
  })
  for (i in seq_along(info))
    if (info[[i]]$type %in% c("numeric", "ordinal") && !length(info[[i]]$sorted))
      stop(sprintf("Predictor '%s' has no observed values.", vars[i]), call. = FALSE)
  structure(list(vars = vars, types = unname(types[vars]), info = info, n = nrow(X)),
            class = "siga_encoding")
}

## Subgroup membership of every individual (column) of the population `Chrom`.
decode_masks <- function(Chrom, enc) {
  m <- nrow(Chrom)
  M <- matrix(TRUE, enc$n, m)
  for (j in seq_along(enc$info)) {
    g <- Chrom[, j]
    if (!any(g > 2)) next
    v <- enc$info[[j]]
    if (v$type == "numeric" || v$type == "ordinal") {
      le <- which(g > 2 & g <= 3)
      if (length(le)) M[, le] <- M[, le] & outer(v$x_le, q7(v$sorted, g[le] - 2), "<=")
      ge <- which(g > 3)
      if (length(ge)) M[, ge] <- M[, ge] & outer(v$x_ge, q7(v$sorted, g[ge] - 3), ">=")
    } else if (v$type == "binary") {
      i1 <- which(g > 2 & g <= 3)
      if (length(i1)) M[, i1] <- M[, i1] & v$is1
      i2 <- which(g > 3)
      if (length(i2)) M[, i2] <- M[, i2] & v$is2
    } else {
      used <- which(g > 2)
      blocks <- ceiling((g[used] - 2) * (v$ncomb / 2))
      for (t in seq_along(used))
        M[, used[t]] <- M[, used[t]] & (v$codes %in% v$combos[[blocks[t]]])
    }
  }
  M
}

## Conditions encoded by one chromosome (numeric vector of genes).
decode_conditions <- function(chrom, enc) {
  conds <- list()
  for (j in seq_along(enc$info)) {
    g <- chrom[j]
    if (g <= 2) next
    v <- enc$info[[j]]
    cond <- switch(v$type,
      numeric = , ordinal = if (g <= 3) {
        list(op = "<=", value = q7(v$sorted, g - 2))
      } else {
        list(op = ">=", value = q7(v$sorted, g - 3))
      },
      binary = list(op = "==", value = if (g <= 3) v$levels[1L] else v$levels[2L]),
      categorical = {
        block <- ceiling((g - 2) * (v$ncomb / 2))
        list(op = "in", value = v$levels[v$combos[[block]]])
      })
    conds[[length(conds) + 1L]] <- c(list(var = v$var, type = v$type), cond,
                                    list(levels = v$levels))
  }
  conds
}
