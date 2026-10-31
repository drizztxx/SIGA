## V5 / V6. Simulation study of the SIGA paper (section 3, Tables 3 to 6).
##
## Usage (from the validation/ directory):
##   Rscript 04_simulation.R endpoint=tte R=100 cores=8 [s=5] [q=10]
##                           [method=exact|bernoulli] [sim=proper|legacy]
##                           [sizes=300,600] [models=A1,A2,...] [cache=<dir>]
## The defaults are the package defaults: s = 5 sub-populations of the genetic
## algorithm per resample (control$n_subpop), q = 10 resamples and exactly
## 63.2% of the patients per resample (method = "exact"; "bernoulli" is the
## sampling of the original code). Time-to-event data are generated with the
## observed time min(T, C) as described in the paper (sim = "proper",
## siga_simulate(legacy = FALSE)); sim = "legacy" uses the original simulation
## program (siga_simulate(legacy = TRUE)), which keeps the failure time of
## censored patients. Runs with other settings are cached separately (label
## suffixes "_s<s>", "_q<q>", "_bernoulli", "_legacydata").
##
## Scenarios: models A1, A2, B1-B4 x training size 300 / 600, 20 covariates.
## Censoring is about 30%; binary and continuous data come from the analogous
## models. Every trial generates a training set L1 and a test set L2 of 5000
## patients, runs SIGA with the settings of the paper (k = 10 candidates,
## m = 50, MaxG = 50, s = 5, GGAP = 0.9, MinS = 10%) and evaluates the 9
## combinations eta in {0.15, 0.2, 0.25} x theta in {0.5, 0.7, 0.9} from the
## same resampling with update(). One result file per trial is written to the
## cache directory, so the run can be interrupted and resumed.

suppressPackageStartupMessages(library(SIGA))
kv <- commandArgs(trailingOnly = TRUE)
opt <- setNames(sub("^[^=]*=", "", kv), sub("=.*$", "", kv))
endpoint <- if ("endpoint" %in% names(opt)) opt[["endpoint"]] else "tte"
n_rep <- if ("R" %in% names(opt)) as.integer(opt[["R"]]) else 100L
n_cores <- if ("cores" %in% names(opt)) as.integer(opt[["cores"]]) else
  max(1L, parallel::detectCores() - 1L)
sizes <- if ("sizes" %in% names(opt)) as.integer(strsplit(opt[["sizes"]], ",")[[1]]) else c(300L, 600L)
n_subpop <- if ("s" %in% names(opt)) as.integer(opt[["s"]]) else 5L
n_resample <- if ("q" %in% names(opt)) as.integer(opt[["q"]]) else 10L
method <- if ("method" %in% names(opt)) opt[["method"]] else "exact"
stopifnot(method %in% c("exact", "bernoulli"))
sim <- if ("sim" %in% names(opt)) opt[["sim"]] else "proper"
stopifnot(sim %in% c("legacy", "proper"))
label <- paste0(endpoint, if (n_subpop != 5L) paste0("_s", n_subpop),
                if (n_resample != 10L) paste0("_q", n_resample),
                if (method == "bernoulli") "_bernoulli",
                if (endpoint == "tte" && sim == "legacy") "_legacydata")
cache <- if ("cache" %in% names(opt)) opt[["cache"]] else
  file.path(Sys.getenv("LOCALAPPDATA", tempdir()), "SIGA-validation", label)
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
resampling <- list(resample_subpop = n_subpop, resample_method = method)

models <- c("A1", "A2", "B1", "B2", "B3", "B4")
run_models <- if ("models" %in% names(opt)) strsplit(opt[["models"]], ",")[[1]] else models
stopifnot(all(run_models %in% models))
size_index <- c(`300` = 1L, `600` = 2L)
grid_eta <- c(0.15, 0.2, 0.25)
grid_theta <- c(0.5, 0.7, 0.9)

tasks <- list()
for (mi in which(models %in% run_models)) for (n in sizes) for (r in seq_len(n_rep)) {
  ni <- size_index[[as.character(n)]]
  file <- file.path(cache, sprintf("%s_n%d_r%03d.rds", models[mi], n, r))
  if (!file.exists(file))
    tasks[[length(tasks) + 1L]] <- list(model = models[mi], n = n, rep = r, file = file,
                                        seed = 1e5 * (10 * mi + ni) + 10 * r)
}
cat(sprintf("[%s] %s (s = %d, q = %d, %s resampling%s): %d trials to run on %d cores (cache %s)\n",
            format(Sys.time()), endpoint, n_subpop, n_resample, method,
            if (endpoint == "tte") paste0(", ", sim, " data") else "", length(tasks), n_cores,
            if (nzchar(Sys.getenv("LOCALAPPDATA")))
              sub(Sys.getenv("LOCALAPPDATA"), "%LOCALAPPDATA%", cache, fixed = TRUE) else cache))

run_trial <- function(task, endpoint, grid_eta, grid_theta, resampling, n_resample, sim) {
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch({
    legacy <- endpoint == "tte" && sim == "legacy"
    L1 <- siga_simulate(task$n, task$model, endpoint, legacy = legacy, seed = task$seed)
    L2 <- siga_simulate(5000, task$model, endpoint, legacy = legacy, seed = task$seed + 5)
    f <- if (endpoint == "tte") Surv(time, status) ~ trt | . else y ~ trt | .
    fit <- siga(f, L1, eta = grid_eta[1], theta = grid_theta[1], n_resample = n_resample,
                seed = task$seed, verbose = FALSE, resample_subpop = resampling$resample_subpop,
                resample_method = resampling$resample_method)
    truth <- attr(L2, "true_subgroup")
    pv <- attr(L1, "predictive")
    ep2 <- SIGA:::prepare_endpoint(f, L2, quiet = TRUE)$ep
    rows <- list()
    cur <- fit
    for (th in grid_theta) for (et in grid_eta) {
      cur <- update(cur, eta = et, theta = th)
      row <- data.frame(model = task$model, n = task$n, rep = task$rep, eta = et, theta = th,
                        found = cur$subgroup, vars = NA_character_, one = NA, all = NA, exact = NA,
                        se = NA_real_, sp = NA_real_, auc = NA_real_, prop_L2 = NA_real_,
                        est_L1 = NA_real_, est_L2 = NA_real_, rule = NA_character_,
                        stringsAsFactors = FALSE)
      if (cur$subgroup) {
        v <- sort(rule_variables(cur$rule))
        m2 <- rule_members(cur$rule, L2)
        row$vars <- paste(v, collapse = " ")
        row$rule <- format(cur$rule)
        if (length(pv)) {
          row$one <- any(pv %in% v); row$all <- all(pv %in% v); row$exact <- setequal(v, pv)
          row$se <- mean(m2[truth]); row$sp <- mean(!m2[!truth]); row$auc <- (row$se + row$sp) / 2
        }
        row$prop_L2 <- mean(m2)
        row$est_L1 <- cur$effects$table$estimate[2]
        row$est_L2 <- SIGA:::effect_summary(ep2, m2)$table$estimate[2]
      }
      rows[[length(rows) + 1L]] <- row
    }
    list(rows = do.call(rbind, rows), mean_ratio = fit$decision$mean_ratio,
         overall_ratio = fit$overall_ratio, n_ok = fit$decision$n_ok, frequency = fit$frequency,
         resamples = fit$resamples, error = NA_character_)
  }, error = function(e) list(rows = NULL, error = conditionMessage(e)))
  res$task <- task[c("model", "n", "rep", "seed")]
  res$endpoint <- endpoint
  res$q <- n_resample
  res$sim <- sim
  res$resampling <- resampling
  res$seconds <- proc.time()[["elapsed"]] - t0
  saveRDS(res, task$file)
  res$seconds
}

if (length(tasks)) {
  cl <- parallel::makePSOCKcluster(n_cores)
  invisible(parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(SIGA))))
  ## one trial per chunk: trials are processed in order and the load is balanced
  secs <- parallel::parLapplyLB(cl, tasks, run_trial, endpoint = endpoint, grid_eta = grid_eta,
                                grid_theta = grid_theta, resampling = resampling,
                                n_resample = n_resample, sim = sim, chunk.size = 1L)
  parallel::stopCluster(cl)
  cat(sprintf("[%s] done: mean %.1f s per trial\n", format(Sys.time()), mean(unlist(secs))))
}
