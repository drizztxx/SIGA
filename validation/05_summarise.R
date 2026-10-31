## Summarise the simulation results (04_simulation.R) in the format of Tables 3
## to 6 of the SIGA paper and compare them with the published values.
##
## Usage (from the validation/ directory): Rscript 05_summarise.R [labels...]
## Without arguments every run in the cache is summarised (labels as in
## 04_simulation.R: "tte", "binary" and "continuous" for the package defaults).

args <- commandArgs(trailingOnly = TRUE)
cache_root <- file.path(Sys.getenv("LOCALAPPDATA", tempdir()), "SIGA-validation")
endpoints <- if (length(args)) args else sort(basename(list.dirs(cache_root, recursive = FALSE)))
paper <- read.csv("paper_tables.csv")

load_trials <- function(endpoint) {
  files <- list.files(file.path(cache_root, endpoint), pattern = "\\.rds$", full.names = TRUE)
  trials <- lapply(files, readRDS)
  err <- vapply(trials, function(t) !is.na(t$error), logical(1))
  list(rows = do.call(rbind, lapply(trials[!err], `[[`, "rows")),
       n_error = sum(err), seconds = vapply(trials, `[[`, numeric(1), "seconds"),
       errors = unique(vapply(trials[err], `[[`, character(1), "error")),
       trials = trials[!err])
}

## Number of resamples. For trials with 100 stored resamples: the decisions
## that q < 100 resamples would give. Resample h uses random number stream h for
## every q, so the first q resamples are exactly those of a run with
## n_resample = q and the same seed ("first"); "random" averages the decision
## over random subsets of q resamples.
rule_vars <- function(rule) unique(sub("^\\s*\\(?\\s*([A-Za-z0-9_.]+).*$", "\\1",
                                       trimws(strsplit(rule, "&", fixed = TRUE)[[1]])))
resample_number <- function(trials, qs = c(10, 20, 50, 100), n_random = 200) {
  grid <- expand.grid(eta = c(0.05, 0.1, 0.15, 0.2, 0.25), theta = c(0.5, 0.7, 0.9))
  decide <- function(rs, M, i, overall) {
    mr <- if (length(i)) mean(rs$test_ratio[i], na.rm = TRUE) else NA_real_
    mf <- if (length(i)) max(colMeans(M[i, , drop = FALSE])) else 0
    !is.na(mr) & mr < 1 & mr <= (1 - grid$eta) * overall & mf >= grid$theta
  }
  set.seed(1)
  out <- lapply(trials, function(t) {
    if (is.null(t$resamples) || nrow(t$resamples) < max(qs)) return(NULL)
    ids <- t$resamples$id
    rs <- t$resamples[t$resamples$status == "ok", ]
    V <- lapply(rs$rule, rule_vars)
    M <- do.call(rbind, lapply(V, function(v) sort(unique(unlist(V))) %in% v))
    do.call(rbind, lapply(qs, function(q) {
      random <- numeric(nrow(grid))
      for (b in seq_len(if (q < length(ids)) n_random else 1L))
        random <- random + decide(rs, M, which(rs$id %in% sample(ids, q)), t$overall_ratio)
      data.frame(model = t$task$model, n = t$task$n, rep = t$task$rep, grid, q = q,
                 first = as.numeric(decide(rs, M, which(rs$id <= q), t$overall_ratio)),
                 random = random / if (q < length(ids)) n_random else 1)
    }))
  })
  d <- do.call(rbind, out)
  if (is.null(d)) return(NULL)
  agg <- aggregate(cbind(first, random) ~ model + n + eta + theta + q, d, mean)
  agg$n_trials <- aggregate(first ~ model + n + eta + theta + q, d, length)$first
  agg[order(agg$model, agg$n, agg$eta, agg$theta, agg$q), ]
}

summarise_rows <- function(rows) {
  key <- unique(rows[c("model", "n", "eta", "theta")])
  out <- lapply(seq_len(nrow(key)), function(i) {
    k <- key[i, ]
    r <- rows[rows$model == k$model & rows$n == k$n & rows$eta == k$eta & rows$theta == k$theta, ]
    f <- r[r$found, ]
    nt <- nrow(r); nf <- nrow(f)
    mk <- function(metric, value, sd = NA_real_, mcse = NA_real_)
      data.frame(k, metric = metric, value = value, sd = sd, mc_se = mcse, n_trials = nt, n_found = nf)
    p <- mean(r$found)
    res <- list(mk(if (k$model %in% c("A1", "A2")) "fdr" else "discovery", p,
                   mcse = sqrt(p * (1 - p) / nt)))
    if (!k$model %in% c("A1", "A2") && nf > 0) {
      for (m in c("one", "all", "exact")) {
        q <- mean(f[[m]])
        res[[length(res) + 1L]] <- mk(paste0("p_", m), q, mcse = sqrt(q * (1 - q) / nf))
      }
      rel <- (f$est_L1 - f$est_L2) / f$est_L2
      for (m in c("se", "sp", "auc", "prop_L2", "est_L1", "est_L2")) {
        v <- f[[m]]
        res[[length(res) + 1L]] <- mk(m, mean(v, na.rm = TRUE), stats::sd(v, na.rm = TRUE),
                                      stats::sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))
      }
      res[[length(res) + 1L]] <- mk("rel_bias", mean(rel, na.rm = TRUE), stats::sd(rel, na.rm = TRUE),
                                    stats::sd(rel, na.rm = TRUE) / sqrt(sum(!is.na(rel))))
    }
    do.call(rbind, res)
  })
  do.call(rbind, out)
}

dir.create("results", showWarnings = FALSE)
for (ep in endpoints) {
  if (!dir.exists(file.path(cache_root, ep))) next
  tr <- load_trials(ep)
  if (is.null(tr$rows)) next
  s <- summarise_rows(tr$rows)
  s$endpoint <- ep
  if (startsWith(ep, "tte")) {
    s <- merge(s, paper[c("eta", "theta", "n", "model", "metric", "value", "se")],
               by = c("eta", "theta", "n", "model", "metric"), all.x = TRUE,
               suffixes = c("", "_paper"))
    names(s)[names(s) == "se"] <- "sd_paper"
    s$diff <- s$value - s$value_paper
  }
  s <- s[order(s$metric, s$model, s$n, s$eta, s$theta), ]
  saveRDS(list(summary = s, rows = tr$rows, n_error = tr$n_error, errors = tr$errors,
               seconds = tr$seconds), sprintf("results/04_simulation_%s.rds", ep))
  write.csv(s, sprintf("results/04_simulation_%s_summary.csv", ep), row.names = FALSE)
  cat(sprintf("\n=== %s: %d trials (%d errors), mean %.0f s per trial ===\n", ep,
              length(tr$seconds), tr$n_error, mean(tr$seconds)))
  fd <- s[s$metric %in% c("fdr", "discovery") & s$eta == 0.2 & s$theta == 0.5, ]
  print(fd[intersect(c("model", "n", "metric", "value", "mc_se", "value_paper"), names(fd))],
        row.names = FALSE, digits = 3)
  rn <- resample_number(tr$trials)
  if (!is.null(rn)) {
    saveRDS(rn, sprintf("results/04_resample_number_%s.rds", ep))
    cat("Rates with the first q of the 100 resamples (eta = 0.2, theta = 0.5):\n")
    x <- rn[rn$eta == 0.2 & rn$theta == 0.5, ]
    print(reshape(x[c("model", "n", "q", "first", "n_trials")], idvar = c("model", "n", "n_trials"),
                  timevar = "q", direction = "wide"), row.names = FALSE, digits = 3)
  }
}
