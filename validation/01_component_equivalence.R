## V2. Component-level equivalence of the SIGA package with the original code.
##
## Run from the validation/ directory:  Rscript 01_component_equivalence.R
## Checks, on identical inputs and on code paths where the original code has no
## known bug:
##   (a) subgroup decoding (legacy decodsubset vs SIGA decoder);
##   (b) objective values (legacy objF.cox3 vs SIGA "interaction" objective);
##   (c) TTE prescreening p-values and ranking (legacy VariableScreen method 2);
##   (d) the complete genetic algorithm trajectory (legacy SIGA.tte vs ga_run)
##       with the same random seed;
##   (e) the simulation program (legacy simul.Cox vs siga_simulate(legacy = TRUE)).

suppressPackageStartupMessages({ library(SIGA); library(survival) })
source("R/legacy.R")
L <- load_legacy()
out <- list()
report <- function(name, ok, detail = "") {
  cat(sprintf("%-60s %s %s\n", name, if (ok) "PASS" else "FAIL", detail))
  out[[name]] <<- list(ok = ok, detail = detail)
}

dat <- siga_simulate(300, model = "B3", legacy = TRUE, seed = 101)
prep <- SIGA:::prepare_siga_data(Surv(time, status) ~ trt | ., dat, quiet = TRUE)
legacy_dat <- data.frame(Ts = dat$time, event = dat$status, z = dat$trt, dat[paste0("x", 1:20)])

## (a) decoding ---------------------------------------------------------------
## Non-buggy paths of decodsubset: >= 2 numeric genes per direction or none,
## >= 2 binary genes per level or none.
vars <- c("x1", "x2", "x4", "x9", "x5", "x6")
enc <- SIGA:::encode_predictors(prep$X, vars, prep$types)
set.seed(1)
Chrom <- matrix(runif(4000 * length(vars), 0, 4), ncol = length(vars))
num <- Chrom[, 1:4]; bin <- Chrom[, 5:6]
okrow <- function(g, lo, hi) { k <- rowSums(g > lo & g <= hi); k == 0 | k >= 2 }
keep <- okrow(num, 2, 3) & okrow(num, 3, 4) & okrow(bin, 2, 3) & okrow(bin, 3, 4)
Chrom <- Chrom[keep, ][1:300, ]
M_new <- SIGA:::decode_masks(Chrom, enc)
leg <- L$decodsubset(data = legacy_dat[vars], Chrom = Chrom)
M_old <- vapply(seq_len(nrow(Chrom)), function(i) {
  s <- leg$subgroup.desc[i]
  if (s == "no subgroup") rep(TRUE, nrow(dat)) else {
    m <- with(legacy_dat, eval(parse(text = s))); m & !is.na(m)
  }
}, logical(nrow(dat)))
n_diff <- sum(colSums(M_new != M_old) > 0)
report("(a) decoding: identical membership (300 chromosomes)", n_diff == 0,
       sprintf("[%d differing]", n_diff))

## (b) objective values ------------------------------------------------------
obj_new <- SIGA:::make_objective(prep$ep, "interaction", 0.1)
legacy_rules <- leg$subgroup.desc
p_old <- L$objF.cox3(legacy_rules, legacy_dat, cutoff = 10, trt = "z=1", ctrl = "z=0",
                     status = "event", time = "Ts")$objv
## legacy transform of the two-sided p-value, from the SIGA z statistic
z_new <- vapply(seq_len(ncol(M_new)), function(i) {
  m <- M_new[, i]
  if (obj_new(m) == 1) return(NA_real_)
  I <- as.numeric(m)
  s <- SIGA:::coef_stat(prep$ep, cbind(prep$ep$z, I, prep$ep$z * I), 3L)
  s[1] / s[2]
}, numeric(1))
p_legacy_form <- ifelse(is.na(z_new), 1, {
  p2 <- 2 * pnorm(-abs(z_new)); ifelse(z_new > 0 & p2 < 0.5, 1 - p2, p2)
})
maxdiff <- max(abs(p_legacy_form - p_old))
report("(b) objective: legacy p-values reproduced from SIGA z", maxdiff < 1e-10,
       sprintf("[max |diff| %.2e]", maxdiff))
p_new <- vapply(seq_len(ncol(M_new)), function(i) obj_new(M_new[, i]), numeric(1))
benef <- !is.na(z_new) & z_new < 0
same_rank <- identical(rank(p_new[benef]), rank(p_old[benef]))
report("(b) objective: identical ranking of beneficial subgroups", same_rank,
       sprintf("[%d subgroups with z < 0]", sum(benef)))
report("(b) objective: same penalised (p = 1) subgroups", identical(p_new == 1, p_old == 1))

## (c) screening ----------------------------------------------------------
vs_old <- L$VariableScreen(legacy_dat, trt = "z=1", ctrl = "z=0", status = "event", time = "Ts",
                           varScreen = c(2, 10), resptyp = "tte")
legacy_dat_f <- legacy_dat
legacy_dat_f$x5 <- factor(legacy_dat_f$x5); legacy_dat_f$x6 <- factor(legacy_dat_f$x6)
vs_new <- siga_screen(Surv(time, status) ~ trt | ., data = dat)$table
pd <- max(abs(vs_new$p_value[match(names(vs_old$vim), vs_new$variable)] - vs_old$vim))
report("(c) screening: identical LRT p-values (20 predictors)", pd < 1e-8, sprintf("[max |diff| %.2e]", pd))
report("(c) screening: identical top-10 ranking", identical(vs_old$col.vim[1:10], vs_new$variable[1:10]))

## (d) GA trajectory -------------------------------------------------------
## One predictor of each type, so that every chromosome decodes without a
## legacy bug (the categorical predictor has 4 levels: the legacy decoder is
## wrong for exactly 3 levels, see the validation report); legacy objective
## transform reproduced in the SIGA engine.
dg <- dat[c("time", "status", "trt", "x4", "x5")]
dg$m4 <- factor(cut(dat$x1, c(-Inf, -0.7, 0, 0.7, Inf), labels = c("a", "b", "c", "d")))
prep_g <- SIGA:::prepare_siga_data(Surv(time, status) ~ trt | x4 + x5 + m4, dg, quiet = TRUE)
enc_g <- SIGA:::encode_predictors(prep_g$X, c("x4", "x5", "m4"), prep_g$types)
obj_g <- SIGA:::make_objective(prep_g$ep, "interaction", 0.1)
## Both implementations round the objective to 10 significant digits: values
## that are equal in exact arithmetic (e.g. a subgroup and its complement under
## the legacy transform) then tie exactly in both, so that gatbxr::ranking()
## (which averages the fitness of exactly tied values) behaves identically.
legacy_obj <- function(mask) {
  if (obj_g(mask) == 1) return(1)
  I <- as.numeric(mask)
  s <- SIGA:::coef_stat(prep_g$ep, cbind(prep_g$ep$z, I, prep_g$ep$z * I), 3L)
  if (is.null(s)) return(1)
  z <- s[1] / s[2]; p2 <- 2 * pnorm(-abs(z))
  signif(if (z > 0 && p2 < 0.5) 1 - p2 else p2, 10)
}
L$objF.cox3_r10 <- function(...) {
  r <- L$objF.cox3(...)
  r$objv <- signif(r$objv, 10)
  r
}
ctrl <- siga_control(pop_size = 50, max_gen = 50, n_subpop = 1)
legacy_g <- data.frame(Ts = dg$time, event = dg$status, z = dg$trt, x4 = dg$x4, x5 = dg$x5, m4 = dg$m4)
traj_ok <- logical(5)
for (s in 1:5) {
  set.seed(s)
  old <- L$SIGA.tte(legacy_g, trt = "z=1", ctrl = "z=0", status = "event", time = "Ts",
                    varused = c("x4", "x5", "m4"), popSize = 50, objFunction = "objF.cox3_r10",
                    MAXGEN = 50, cutoff = 10, verbose = FALSE)
  set.seed(s)
  new <- SIGA:::ga_run(prep_g$ep, prep_g$X, enc_g, ctrl, objective = legacy_obj)
  traj_ok[s] <- isTRUE(all.equal(old$Best, new$history[1:50], tolerance = 1e-10))
}
report("(d) GA: identical best-objective trajectory (5 seeds x 50 generations)", all(traj_ok),
       sprintf("[%d/5 identical]", sum(traj_ok)))

## (e) simulation program --------------------------------------------------
same_sim <- vapply(c("A1", "A2", "B1", "B2", "B3", "B4"), function(m) {
  ex <- c(A1 = "0", A2 = "2*(x1+x2)", B1 = "0.1*(x1+x2)+2*z*(x4>1)", B2 = "0.2*(x1+x2)+1.2*z*(x1>0)",
          B3 = "0.3*(x2+x9)+2*z*(x4>1 & x5==1)", B4 = "2*z*(x1<0.5 & x5==1 & x6==1)")[[m]]
  pi_m <- c(A1 = 3.65, A2 = 4.42, B1 = 4.05, B2 = 3.97, B3 = 3.86, B4 = 3.84)[[m]]
  o <- L$simul.Cox(n = 300, exp = ex, weib.shape = 4.34, weib.scale = 1.2, censor.par = pi_m,
                   extra.noise = 10, seed = 7)
  n <- siga_simulate(300, model = m, legacy = TRUE, seed = 7)
  identical(o$Ts, n$time) && identical(o$event, n$status) &&
    all(vapply(paste0("x", 1:20), function(v) identical(as.character(o[[v]]), as.character(n[[v]])), TRUE))
}, logical(1))
report("(e) simulation: siga_simulate(legacy = TRUE) == simul.Cox (6 models)", all(same_sim))

dir.create("results", showWarnings = FALSE)
saveRDS(out, "results/01_component_equivalence.rds")
cat(sprintf("\n%d of %d checks passed\n", sum(vapply(out, `[[`, TRUE, "ok")), length(out)))
