## Single-population genetic algorithm ----------------------------------------
##
## One run of the SIGA genetic algorithm on the predictors encoded in `enc`,
## using the current random number stream. The sequence of gatbxr operators is
## that of the original implementation (SIGA.tte / SIGA.binary):
##   crtrp -> [ranking -> select -> recombin -> mutate -> evaluate -> reins] x max_gen
## with elitist reinsertion (reins(InsOpt = c(1, 1))): the offspring replace
## the worst parents, so the best individual always survives.

## `cache` (see new_objective_cache()) may be shared by runs on the same data:
## objective values are a deterministic function of the subgroup membership,
## so sharing only saves time and never changes the result.
ga_run <- function(ep, X, enc, control, objective = NULL, cache = NULL) {
  k <- length(enc$vars)
  objective <- objective %||% make_objective(ep, control$objective, control$min_size)
  penalty <- if (is.function(control$objective)) Inf else 1
  FieldDR <- matrix(c(0, 4), 2L, k)
  if (is.null(cache)) cache <- new_objective_cache()
  n_eval <- 0L
  evaluate <- function(Chrom) {
    M <- decode_masks(Chrom, enc)
    keys <- mask_keys(M)
    idx <- match(keys, cache$keys)
    miss <- which(is.na(idx))
    if (length(miss)) {
      ## new subgroups, each evaluated once (also when it occurs repeatedly)
      first <- miss[!duplicated(keys[miss])]
      cache$keys <- c(cache$keys, keys[first])
      cache$vals <- c(cache$vals, vapply(first, function(j) objective(M[, j]), numeric(1)))
      n_eval <<- n_eval + length(first)
      idx[miss] <- match(keys[miss], cache$keys)
    }
    list(objv = cache$vals[idx], keys = keys, M = M)
  }
  keep <- function(Ch, ev, sel) list(chrom = Ch[sel, , drop = FALSE], keys = ev$keys[sel],
                                    M = ev$M[, sel, drop = FALSE])

  Chrom <- gatbxr::crtrp(control$pop_size, FieldDR)
  ev <- evaluate(Chrom)
  ObjV <- ev$objv
  best <- min(ObjV)
  elite <- keep(Chrom, ev, ObjV == best)
  history <- numeric(control$max_gen + 1L)

  for (gen in seq_len(control$max_gen)) {
    history[gen] <- min(ObjV)
    FitnV <- gatbxr::ranking(ObjV)
    SelCh <- gatbxr::select(control$selection, Chrom, FitnV, control$ggap)
    SelCh <- gatbxr::recombin(control$recombination, SelCh)
    SelCh <- gatbxr::mutate(control$mutation, SelCh, FieldDR = FieldDR, MutOpt = control$mutation_rate)
    evs <- evaluate(SelCh)
    rs <- gatbxr::reins(Chrom, SelCh, 1, c(1, 1), ObjV, evs$objv)
    Chrom <- rs$Chrom
    ObjV <- rs$ObjVCh
    mo <- min(evs$objv)
    if (mo < best) {
      best <- mo
      elite <- keep(SelCh, evs, evs$objv == best)
    } else if (mo == best) {
      add <- keep(SelCh, evs, evs$objv == best)
      elite <- list(chrom = rbind(elite$chrom, add$chrom), keys = c(elite$keys, add$keys),
                    M = cbind(elite$M, add$M))
    }
  }
  history[control$max_gen + 1L] <- min(ObjV)

  out <- list(objective = best, valid = best < penalty, rule = NULL, mask = NULL, chrom = NULL,
              history = if (control$keep_history) history else NULL, n_eval = n_eval,
              vars = enc$vars)
  if (out$valid) out[c("rule", "mask", "chrom")] <- pick_best(elite, ep, X, enc)
  out
}

## Among the individuals with the best objective value choose the largest
## subgroup (2 * min(n_treated, n_control), then total size) and, among its
## representations, the one with the fewest conditions after simplification.
pick_best <- function(elite, ep, X, enc) {
  ukeys <- unique(elite$keys)
  first <- match(ukeys, elite$keys)
  bal <- size <- numeric(length(ukeys))
  for (i in seq_along(ukeys)) {
    m <- elite$M[, first[i]]
    n1 <- sum(m & ep$z == 1)
    size[i] <- sum(m)
    bal[i] <- 2 * min(n1, size[i] - n1)
  }
  key <- ukeys[order(-bal, -size)[1L]]
  rows <- utils::head(which(elite$keys == key), 25L)
  ## identical chromosomes can be stored repeatedly; decode each only once
  rows <- rows[!duplicated(elite$chrom[rows, , drop = FALSE])]
  cands <- lapply(rows, function(r) chrom_rule(elite$chrom[r, ], enc, X))
  ncond <- vapply(cands, function(cr) length(cr$rule$conditions), integer(1))
  i <- which.min(ncond)
  list(rule = cands[[i]]$rule, mask = cands[[i]]$mask, chrom = elite$chrom[rows[i], ])
}
