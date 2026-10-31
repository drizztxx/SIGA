## V3. Search quality of the SIGA package versus the original code.
##
## Run from the validation/ directory:  Rscript 02_ga_comparison.R [n_datasets]
## On the same simulated B3 datasets (n = 300, 20 covariates, original
## simulation program) and the same 10 prescreened candidates, one
## single-population genetic algorithm run (m = 50, MaxG = 50) of the original
## SIGA.tte() is compared with one run of the package. Both best subgroups are
## scored with the same criterion (one-sided interaction p-value) and their
## variables are compared with the true predictive biomarkers x4 and x5.

suppressPackageStartupMessages({ library(SIGA); library(survival) })
source("R/legacy.R")
L <- load_legacy()
args <- commandArgs(trailingOnly = TRUE)
n_data <- if (length(args)) as.integer(args[1]) else 50L
ctrl <- siga_control(pop_size = 50, max_gen = 50, n_subpop = 1)
vars_of <- function(rule_text) {
  if (is.na(rule_text) || rule_text == "no subgroup") return(character())
  sort(unique(sub("^\\(?\\s*([A-Za-z0-9_.]+).*$", "\\1", trimws(strsplit(rule_text, "&")[[1]]))))
}

rows <- lapply(seq_len(n_data), function(d) {
  dat <- siga_simulate(300, model = "B3", legacy = TRUE, seed = 5000 + d)
  cand <- siga_screen(Surv(time, status) ~ trt | ., data = dat)$selected
  prep <- SIGA:::prepare_siga_data(Surv(time, status) ~ trt | ., dat, quiet = TRUE)
  obj <- SIGA:::make_objective(prep$ep, "interaction", 0.1)
  ldat <- data.frame(Ts = dat$time, event = dat$status, z = dat$trt, dat[paste0("x", 1:20)])
  set.seed(d)
  t_old <- system.time(old <- L$SIGA.tte(ldat, trt = "z=1", ctrl = "z=0", status = "event",
                                         time = "Ts", varused = cand, popSize = 50,
                                         objFunction = "objF.cox3", MAXGEN = 50, cutoff = 10,
                                         verbose = FALSE))[["elapsed"]]
  ov <- unlist(old$Objv)
  best_old <- unlist(old$Subset)[which.min(ov)]
  m_old <- if (best_old == "no subgroup") rep(TRUE, nrow(dat)) else
    with(ldat, { m <- eval(parse(text = best_old)); m & !is.na(m) })
  f <- as.formula(paste("Surv(time, status) ~ trt |", paste(cand, collapse = " + ")))
  t_new <- system.time(new <- siga_search(f, data = dat, control = ctrl, seed = d))[["elapsed"]]
  m_new <- predict(new)
  v_old <- vars_of(best_old)
  v_new <- sort(rule_variables(new$rule))
  data.frame(dataset = d, x4_x5_candidates = all(c("x4", "x5") %in% cand),
             p_old = obj(m_old), p_new = obj(m_new),
             p_true = obj(attr(dat, "true_subgroup")),
             vars_old = paste(v_old, collapse = " "), vars_new = paste(v_new, collapse = " "),
             all_old = all(c("x4", "x5") %in% v_old), all_new = all(c("x4", "x5") %in% v_new),
             exact_old = identical(v_old, c("x4", "x5")), exact_new = identical(v_new, c("x4", "x5")),
             time_old = t_old, time_new = t_new)
})
res <- do.call(rbind, rows)
dir.create("results", showWarnings = FALSE)
saveRDS(res, "results/02_ga_comparison.rds")

cat(sprintf("Datasets: %d (x4 and x5 both among the 10 candidates in %d)\n", nrow(res),
            sum(res$x4_x5_candidates)))
cat(sprintf("Median log10 p of best subgroup: original %.2f, package %.2f, true subgroup %.2f\n",
            median(log10(res$p_old)), median(log10(res$p_new)), median(log10(res$p_true))))
cat(sprintf("Package subgroup at least as good (p_new <= p_old): %d of %d\n",
            sum(res$p_new <= res$p_old * (1 + 1e-8)), nrow(res)))
cat(sprintf("Best subgroup contains x4 and x5: original %d, package %d\n", sum(res$all_old), sum(res$all_new)))
cat(sprintf("Best subgroup uses exactly x4 and x5: original %d, package %d\n", sum(res$exact_old),
            sum(res$exact_new)))
cat(sprintf("Mean time per run: original %.1f s, package %.2f s (%.0f x faster)\n",
            mean(res$time_old), mean(res$time_new), mean(res$time_old) / mean(res$time_new)))
