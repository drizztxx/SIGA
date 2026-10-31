## V4. Regression check on the RADIANT case study (SIGA paper, section 4).
##
## The RADIANT data were obtained through ClinicalStudyDataRequest.com and are
## NOT part of the package or the repository. This script reads them from the
## author's local workspace and writes aggregate results only.
##
## Run from the validation/ directory:
##   Rscript 03_radiant_case_study.R <path to 20220212.RData> [n_cores]
## (or set the environment variable SIGA_RADIANT_DATA to the path).
##
## Published / original-code result (eta = 0.25, theta = 0.7, k = 10, q = 100,
## m = 50, MaxG = 50, s = 5): overall HR 0.945, mean test-set HR 0.70,
## predictors EGFR_IHC_PROP and EGFR_IHC_INTENSE, subgroup
## EGFR_IHC_INTENSE <= 1 & EGFR_IHC_PROP <= 60.

suppressPackageStartupMessages(library(SIGA))
args <- commandArgs(trailingOnly = TRUE)
rdata <- if (length(args)) args[1] else Sys.getenv("SIGA_RADIANT_DATA")
if (!nzchar(rdata) || !file.exists(rdata))
  stop("Give the path of the RADIANT workspace (20220212.RData) as the first argument ",
       "or in the environment variable SIGA_RADIANT_DATA.", call. = FALSE)
n_cores <- if (length(args) > 1) as.integer(args[2]) else max(1L, parallel::detectCores() - 1L)

e <- new.env()
load(rdata, envir = e)
d <- e$erlotinib.dfs
d$pmutate <- NULL                    # derived variable, not one of the 21 predictors
legacy_fit <- e$erlotinib.dfs.siga2  # original-code result stored in the workspace
rm(e)
stopifnot(ncol(d) == 24L)
f <- Surv(dfs, dfs.e) ~ trt | .

## The paper used q = 100 resamples for the case study. "default": package
## defaults otherwise (5 sub-populations per resample, exactly 63.2%);
## "bernoulli": the resampling of the original code; "one_subpop": one
## sub-population per resample, as described in the paper; "default_q10": the
## package defaults including q = 10, with more seeds.
runs <- list()
configs <- list(
  default = list(n_resample = 100),
  bernoulli = list(n_resample = 100, resample_method = "bernoulli"),
  one_subpop = list(n_resample = 100, resample_subpop = 1),
  default_q10 = list(n_resample = 10))
seeds <- list(default = 1:3, bernoulli = 1:3, one_subpop = 1:3, default_q10 = 1:10)
for (cfg in names(configs)) for (seed in seeds[[cfg]]) {
  t <- system.time(fit <- do.call(siga, c(list(formula = f, data = d, control_arm = "2",
                                                eta = 0.25, theta = 0.7, n_cores = n_cores,
                                                seed = seed, verbose = FALSE),
                                           configs[[cfg]])))[["elapsed"]]
  fit2 <- update(fit, eta = 0.2, theta = 0.5)
  sub <- if (fit$subgroup) fit$effects$table[fit$effects$table$group == "subgroup", ] else NULL
  runs[[length(runs) + 1L]] <- data.frame(
    config = cfg, seed = seed, n_resample = fit$settings$n_resample,
    resample_subpop = fit$settings$resample_subpop, overall_hr = fit$overall_ratio,
    mean_test_hr = fit$decision$mean_ratio, improvement = fit$decision$improvement,
    freq_ihc_prop = fit$frequency[["egfr_ihc_prop"]], freq_ihc_intense = fit$frequency[["egfr_ihc_intense"]],
    selected = paste(fit$decision$selected, collapse = " "),
    subgroup = if (fit$subgroup) format(fit$rule) else "no subgroup",
    n_subgroup = if (is.null(sub)) NA else sub$n, hr_subgroup = if (is.null(sub)) NA else sub$estimate,
    subgroup_eta0.2_theta0.5 = if (fit2$subgroup) format(fit2$rule) else "no subgroup",
    minutes = t / 60, stringsAsFactors = FALSE)
  cat(sprintf("[%s, seed %d] mean test HR %.3f; selected: %s; subgroup: %s (%.1f min)\n", cfg, seed,
              fit$decision$mean_ratio, paste(fit$decision$selected, collapse = ", "),
              runs[[length(runs)]]$subgroup, t / 60))
}
res <- do.call(rbind, runs)
ref <- data.frame(overall_hr = legacy_fit$overall.hr, mean_test_hr = legacy_fit$meanhr,
                  selected = paste(sort(legacy_fit$varused), collapse = " "),
                  subgroup = legacy_fit$bestsub)
dir.create("results", showWarnings = FALSE)
saveRDS(list(runs = res, original = ref), "results/03_radiant_case_study.rds")
cat("\nOriginal code (stored result):\n")
print(ref, row.names = FALSE)
cat("\nSIGA package:\n")
print(res[c("config", "seed", "overall_hr", "mean_test_hr", "selected", "subgroup", "hr_subgroup")],
      row.names = FALSE, digits = 3)
