## Reproducible (parallel) computing -----------------------------------------
##
## Every independent task (a resample, a sub-population, a cross-validation
## fold) gets its own L'Ecuyer-CMRG random number stream derived from one seed.
## Results are therefore identical whether tasks run sequentially or on any
## number of cores. The random number generator state of the user is restored.

save_rng <- function() {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  list(had = had, seed = if (had) get(".Random.seed", envir = genv) else NULL,
       kind = RNGkind())
}

restore_rng <- function(state) {
  genv <- globalenv()
  if (state$had) {
    assign(".Random.seed", state$seed, envir = genv)
  } else {
    RNGkind(state$kind[1L], state$kind[2L], state$kind[3L])
    if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
  }
}

## Evaluate `expr` with the global RNG state `seed` (a .Random.seed vector).
with_rng_state <- function(seed, expr) {
  state <- save_rng()
  on.exit(restore_rng(state), add = TRUE)
  assign(".Random.seed", seed, envir = globalenv())
  force(expr)
}

## Draw an integer seed from the current RNG (so that set.seed() by the user
## makes a run reproducible) unless one is given.
resolve_seed <- function(seed) {
  if (is.null(seed)) return(sample.int(.Machine$integer.max, 1L))
  check_number(seed, "seed", integer = TRUE)
  as.integer(seed)
}

## `n` independent L'Ecuyer-CMRG streams from an integer seed.
make_streams <- function(n, seed) {
  if (n < 1L) return(list())
  state <- save_rng()
  on.exit(restore_rng(state), add = TRUE)
  set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  s <- get(".Random.seed", envir = globalenv())
  out <- vector("list", n)
  for (i in seq_len(n)) {
    out[[i]] <- s
    s <- parallel::nextRNGStream(s)
  }
  out
}

## lapply over `.x` where task i runs with RNG stream `.streams[[i]]`,
## optionally on a PSOCK cluster with `.n_cores` workers (each loads SIGA).
## Arguments in `...` are passed to `.f`.
siga_lapply <- function(.x, .f, ..., .streams, .n_cores = 1L) {
  stopifnot(length(.streams) == length(.x))
  ## extra arguments travel as one list so that names such as `X` cannot
  ## collide with the arguments of lapply() / parLapplyLB()
  .args <- list(...)
  task <- function(i, .args) with_rng_state(.streams[[i]], do.call(.f, c(list(.x[[i]]), .args)))
  .n_cores <- min(as.integer(.n_cores), length(.x))
  if (.n_cores <= 1L) return(lapply(seq_along(.x), task, .args = .args))
  cl <- parallel::makePSOCKcluster(.n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  ## the workers unserialize SIGA functions, which needs the installed package
  ok <- unlist(parallel::clusterCall(cl, function() {
    suppressPackageStartupMessages(requireNamespace("SIGA", quietly = TRUE))
  }))
  if (!all(ok))
    stop("Parallel computing requires the SIGA package to be installed in the library ",
         "used by the worker processes.", call. = FALSE)
  parallel::parLapplyLB(cl, seq_along(.x), task, .args = .args)
}

## Is SIGA installed in a library (and not only loaded with pkgload)?
siga_installed <- function() {
  length(find.package("SIGA", lib.loc = .libPaths(), quiet = TRUE)) > 0L
}
