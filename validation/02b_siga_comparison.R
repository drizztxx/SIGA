## V3b. The complete SIGA procedure: original code versus package.
##
## The original SIGA() function (../legacy) is run with the settings of the
## published simulations (5 sub-populations and 10 resamples with Bernoulli
## sampling, objF.cox3, varScreen = c(2, 10), eta = 0.2, theta = 0.5) on the
## datasets B1 (n = 300) and B4 (n = 600) of 04_simulation.R (original
## simulation program): trials 1 to 8 of B4 and 1 to <b1> (default 24) of B1. The package results for the same datasets and
## settings come from
##   Rscript 04_simulation.R endpoint=tte R=24 s=5 q=10 method=bernoulli sim=legacy
## (label tte_bernoulli_legacydata), which must be run first.
##
## Usage (from the validation/ directory): Rscript 02b_siga_comparison.R [workers] [b1]
## The parallel calls of the original code are replaced by sequential ones and the
## datasets are distributed over the workers. One dataset takes about 25 minutes;
## the results are cached in %LOCALAPPDATA%/SIGA-validation-original. The
## package is also run with three further resampling seeds on the B1 datasets
## (about 10 minutes on 4 workers).

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args)) as.integer(args[1]) else 3L
n_b1 <- if (length(args) > 1) as.integer(args[2]) else 24L
vroot <- normalizePath(".", winslash = "/")
out_dir <- file.path(Sys.getenv("LOCALAPPDATA", tempdir()), "SIGA-validation-original")
pkg_dir <- file.path(Sys.getenv("LOCALAPPDATA", tempdir()), "SIGA-validation", "tte_bernoulli_legacydata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seed_of <- function(model, r, n) {
  mi <- match(model, c("A1", "A2", "B1", "B2", "B3", "B4"))
  1e5 * (10 * mi + if (n == 300) 1 else 2) + 10 * r
}
tasks <- c(unlist(lapply(1:8, function(r) list(list(model = "B4", n = 600L, r = r),
                                                list(model = "B1", n = 300L, r = r))), recursive = FALSE),
           lapply(seq_len(max(0L, n_b1 - 8L)) + 8L, function(r) list(model = "B1", n = 300L, r = r)))
tasks <- lapply(tasks, function(t) {
  t$seed <- seed_of(t$model, t$r, t$n)
  t$file <- file.path(out_dir, sprintf("legacy_%s_n%d_r%03d.rds", t$model, t$n, t$r))
  t
})
todo <- tasks[!vapply(tasks, function(t) file.exists(t$file), logical(1))]
cat(sprintf("[%s] %d of %d datasets to run with the original code on %d workers\n",
            format(Sys.time()), length(todo), length(tasks), workers))

worker_init <- function(vroot) {
  suppressPackageStartupMessages({ library(SIGA); library(survival) })
  Sys.setenv(SIGA_VALIDATION_DIR = vroot)
  source(file.path(vroot, "R", "legacy.R"))
  L <- load_legacy()
  L$detectCores <- function(...) 64L
  L$makeCluster <- function(...) NULL
  L$stopCluster <- function(...) invisible(NULL)
  L$clusterExport <- function(...) invisible(NULL)
  L$clusterEvalQ <- function(...) invisible(NULL)
  L$parLapply <- L$parLapplyLB <- function(cl, X, fun, ...) lapply(X, fun, ...)
  for (f in c("parSIGA.R", "boot.SIGA.R", "predict.SIGA.R", "freqvars.R", "SIGA.R"))
    sys.source(file.path(vroot, "..", "legacy", f), envir = L)
  assign("L", L, envir = globalenv())
  invisible(TRUE)
}

run_task <- function(t) {
  dat <- SIGA::siga_simulate(t$n, t$model, legacy = TRUE, seed = t$seed)
  ldat <- data.frame(Ts = dat$time, event = dat$status, z = factor(dat$trt, levels = c(1, 0)),
                     dat[paste0("x", 1:20)])
  t0 <- proc.time()[["elapsed"]]
  set.seed(t$seed)
  res <- tryCatch({
    junk <- utils::capture.output(r <- suppressWarnings(
      L$SIGA(ds = ldat, trt = "z=1", ctrl = "z=0", time = "Ts", status = "event",
             objFunction = "objF.cox3", cutoff = 10, varScreen = c(2, 10), resptyp = "tte",
             clnum = 5, popSize = 50, MAXGEN = 50, n.boot = 10, freq = 0.5,
             multiplicity = 0.2, CPUs = 1)))
    list(overall = r$overall.hr, meanhr = r$meanhr, bootres = r$bootres, varused = r$varused,
         bestsub = r$bestsub, error = NA_character_)
  }, error = function(e) list(error = conditionMessage(e)))
  res$seconds <- proc.time()[["elapsed"]] - t0
  res$task <- t
  saveRDS(res, t$file)
  res$seconds
}

## The package with three further resampling seeds on the B1 datasets: how much
## does the final rule depend on the resamples?
seed_runs <- expand.grid(r = seq_len(n_b1), k = 1:3)
seed_runs$file <- file.path(out_dir, sprintf("package_B1_n300_r%03d_seed%d.rds", seed_runs$r, seed_runs$k))
seed_todo <- which(!file.exists(seed_runs$file))
run_seed <- function(s) {
  seed <- 1e5 * (10 * 3 + 1) + 10 * s$r
  dat <- SIGA::siga_simulate(300, "B1", legacy = TRUE, seed = seed)
  fit <- SIGA::siga(Surv(time, status) ~ trt | ., dat, eta = 0.2, theta = 0.5,
                    n_resample = 10, resample_subpop = 5, resample_method = "bernoulli",
                    seed = seed + s$k, verbose = FALSE)
  v <- if (fit$subgroup) sort(SIGA::rule_variables(fit$rule)) else character()
  res <- data.frame(r = s$r, k = s$k, found = fit$subgroup, vars = paste(v, collapse = " "),
                    exact = fit$subgroup && identical(v, "x4"), all = "x4" %in% v,
                    max_noise = max(c(0, fit$frequency[names(fit$frequency) != "x4"])))
  saveRDS(res, s$file)
  invisible(TRUE)
}

if (length(todo) || length(seed_todo)) {
  cl <- parallel::makePSOCKcluster(workers)
  invisible(parallel::clusterCall(cl, worker_init, vroot))
  if (length(todo)) {
    secs <- parallel::parLapplyLB(cl, todo, run_task, chunk.size = 1L)
    cat(sprintf("[%s] done: mean %.0f s per dataset\n", format(Sys.time()), mean(unlist(secs))))
  }
  if (length(seed_todo))
    invisible(parallel::parLapplyLB(cl, lapply(seed_todo, function(i) seed_runs[i, ]), run_seed,
                                    chunk.size = 1L))
  parallel::stopCluster(cl)
}

## Comparison, per dataset and method: mean held-out / overall hazard ratio of
## the resampled subgroups, number of variables per resampled rule, selection
## frequencies of the predictive and of the most frequent noise variable, and
## the final subgroup (eta = 0.2, theta = 0.5).
pv <- list(B1 = "x4", B4 = c("x1", "x5", "x6"))
vars_of <- function(s, sep) {
  if (is.na(s) || s == "no subgroup") return(character())
  unique(sub("^\\s*\\(?\\s*([A-Za-z0-9_.]+).*$", "\\1", trimws(strsplit(s, sep)[[1]])))
}
summ <- function(method, model, r, rules, test_ratio, overall, selected, final_vars, found) {
  V <- lapply(rules, vars_of, sep = if (method == "original") "[&|]" else "&")
  freq <- table(unlist(V)) / length(V)
  p <- pv[[model]]
  fp <- vapply(p, function(v) if (v %in% names(freq)) freq[[v]] else 0, numeric(1))
  noise <- freq[setdiff(names(freq), p)]
  data.frame(method = method, model = model, rep = r, resamples = length(V),
             test_over_overall = mean(test_ratio, na.rm = TRUE) / overall,
             vars_per_rule = mean(lengths(V)), min_freq_predictive = min(fp),
             max_freq_noise = if (length(noise)) max(noise) else 0,
             selected = paste(sort(selected), collapse = " "), found = found,
             final = paste(sort(final_vars), collapse = " "),
             all = found && all(p %in% final_vars), exact = found && setequal(final_vars, p),
             stringsAsFactors = FALSE)
}
rows <- list()
for (t in tasks) {
  if (!file.exists(t$file)) next
  x <- readRDS(t$file)
  if (!is.na(x$error)) { cat("error:", basename(t$file), x$error, "\n"); next }
  pf <- file.path(pkg_dir, sprintf("%s_n%d_r%03d.rds", t$model, t$n, t$r))
  if (!file.exists(pf)) { cat("package result missing:", pf, "\n"); next }
  found <- !identical(as.character(x$bestsub), "no subgroup")
  rows[[length(rows) + 1L]] <- summ("original", t$model, t$r, as.character(x$bootres$subs),
                                    x$bootres$testd, x$overall, x$varused,
                                    if (found) vars_of(as.character(x$bestsub), "[&|]") else character(),
                                    found)
  y <- readRDS(pf)
  rs <- y$resamples[y$resamples$status == "ok", ]
  rr <- y$rows[y$rows$eta == 0.2 & y$rows$theta == 0.5, ]
  rows[[length(rows) + 1L]] <- summ("package", t$model, t$r, rs$rule, rs$test_ratio, y$overall_ratio,
                                    names(y$frequency)[y$frequency >= 0.5],
                                    if (rr$found) strsplit(rr$vars, " ")[[1]] else character(), rr$found)
}
d <- do.call(rbind, rows)
agg <- aggregate(cbind(test_over_overall, vars_per_rule, min_freq_predictive, max_freq_noise, found,
                       all, exact) ~ model + method, d, mean)
agg$datasets <- aggregate(found ~ model + method, d, length)$found
agg$same_decision <- vapply(seq_len(nrow(agg)), function(i) {
  w <- d[d$model == agg$model[i], ]
  o <- w[w$method == "original", ]; p <- w[w$method == "package", ]
  mean(o$found[match(p$rep, o$rep)] == p$found)
}, numeric(1))
## paired comparison of "final rule exactly the predictive biomarkers"
o <- d[d$method == "original", ]; p <- d[d$method == "package", ]
p <- p[match(paste(o$model, o$rep), paste(p$model, p$rep)), ]
paired <- do.call(rbind, lapply(split(seq_len(nrow(o)), o$model), function(i) {
  b <- sum(o$exact[i] & !p$exact[i]); c <- sum(!o$exact[i] & p$exact[i])
  data.frame(model = o$model[i[1]], original_only = b, package_only = c,
             p_value = if (b + c) stats::binom.test(b, b + c)$p.value else 1)
}))
sr <- do.call(rbind, lapply(seed_runs$file[file.exists(seed_runs$file)], readRDS))
seeds <- if (!is.null(sr)) aggregate(cbind(found, all, exact) ~ k, sr, mean)
dir.create("results", showWarnings = FALSE)
saveRDS(list(datasets = d, summary = agg, paired = paired, seed_runs = sr, seeds = seeds),
        "results/02b_siga_comparison.rds")
options(width = 200)
print(d, digits = 3, row.names = FALSE)
cat("\nMeans over the datasets (all / exact: proportion of datasets whose final rule uses all /",
    "exactly the predictive biomarkers):\n")
print(agg, digits = 3, row.names = FALSE)
cat("\nDatasets whose final rule is exactly the predictive biomarkers with one code only:\n")
print(paired, digits = 2, row.names = FALSE)
if (!is.null(seeds)) {
  cat("\nPackage on the B1 datasets with three further resampling seeds (proportions):\n")
  print(seeds, digits = 3, row.names = FALSE)
}
