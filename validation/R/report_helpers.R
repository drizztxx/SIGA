## Tables for the validation report and the manual.

f3 <- function(x, d = 3) ifelse(is.na(x), "-", formatC(x, digits = d, format = "f"))

## "ours (paper)" cells
pair <- function(ours, paper, d = 3) {
  ifelse(is.na(paper), f3(ours, d), paste0(f3(ours, d), " (", f3(paper, d), ")"))
}

sim_summary <- function(endpoint, dir = "results") {
  f <- file.path(dir, sprintf("04_simulation_%s.rds", endpoint))
  if (!file.exists(f)) return(NULL)
  readRDS(f)
}

## Table 3: false discovery rates (A1, A2) by eta, theta and n.
table3 <- function(s, with_paper = TRUE) {
  x <- s[s$metric == "fdr", ]
  key <- unique(x[c("eta", "theta")])
  key <- key[order(key$eta, key$theta), ]
  out <- data.frame(eta = key$eta, theta = key$theta)
  for (n in c(300, 600)) for (m in c("A1", "A2")) {
    v <- vapply(seq_len(nrow(key)), function(i) {
      r <- x[x$eta == key$eta[i] & x$theta == key$theta[i] & x$n == n & x$model == m, ]
      if (!nrow(r)) return("-")
      if (with_paper && "value_paper" %in% names(r)) pair(r$value, r$value_paper) else f3(r$value)
    }, character(1))
    out[[sprintf("%s n=%d", m, n)]] <- v
  }
  out
}

## Discovery rates of the non-null models.
table_discovery <- function(s) {
  x <- s[s$metric == "discovery", ]
  key <- unique(x[c("eta", "theta")])
  key <- key[order(key$eta, key$theta), ]
  out <- data.frame(eta = key$eta, theta = key$theta)
  for (n in c(300, 600)) for (m in c("B1", "B2", "B3", "B4")) {
    out[[sprintf("%s n=%d", m, n)]] <- vapply(seq_len(nrow(key)), function(i) {
      r <- x[x$eta == key$eta[i] & x$theta == key$theta[i] & x$n == n & x$model == m, ]
      if (nrow(r)) f3(r$value) else "-"
    }, character(1))
  }
  out
}

## Tables 4 to 6: conditional metrics of the non-null models for one n.
table_metrics <- function(s, metrics, n, with_paper = TRUE, d = 3) {
  x <- s[s$metric %in% metrics & s$n == n, ]
  key <- unique(x[c("eta", "theta")])
  key <- key[order(key$eta, key$theta), ]
  out <- data.frame(eta = key$eta, theta = key$theta)
  for (m in c("B1", "B2", "B3", "B4")) for (mt in metrics) {
    out[[paste(m, mt)]] <- vapply(seq_len(nrow(key)), function(i) {
      r <- x[x$eta == key$eta[i] & x$theta == key$theta[i] & x$model == m & x$metric == mt, ]
      if (!nrow(r)) return("-")
      if (with_paper && "value_paper" %in% names(r)) pair(r$value, r$value_paper, d) else f3(r$value, d)
    }, character(1))
  }
  out
}

## Agreement of the false discovery rates with the paper: number of cells
## within two Monte Carlo standard errors.
fdr_agreement <- function(s) {
  x <- s[s$metric == "fdr" & !is.na(s$value_paper), ]
  se <- sqrt(pmax(x$value_paper, 1 / x$n_trials) * (1 - x$value_paper) / x$n_trials)
  c(cells = nrow(x), within_2se = sum(abs(x$value - x$value_paper) <= 2 * se),
    mean_ours = mean(x$value), mean_paper = mean(x$value_paper))
}

## Agreement with the paper by metric: mean over cells and number of cells
## within two Monte Carlo standard errors (binomial standard error with the
## larger of the two proportions, otherwise the Monte Carlo standard error of
## the package value).
agreement <- function(s, metrics = c("fdr", "p_one", "p_all", "p_exact", "se", "sp", "auc",
                                     "prop_L2", "est_L1", "est_L2", "rel_bias")) {
  x <- s[!is.na(s$value_paper) & s$metric %in% metrics, ]
  prop <- x$metric %in% c("fdr", "p_one", "p_all", "p_exact")
  nn <- ifelse(x$metric == "fdr", x$n_trials, x$n_found)
  se <- x$mc_se
  se[prop] <- sqrt(pmax(x$value * (1 - x$value), x$value_paper * (1 - x$value_paper))[prop] / nn[prop])
  ok <- ifelse(x$value == x$value_paper, TRUE, abs(x$value - x$value_paper) <= 2 * se)
  do.call(rbind, lapply(intersect(metrics, unique(x$metric)), function(m) {
    i <- x$metric == m & !is.na(ok)
    data.frame(metric = m, cells = sum(i), package = mean(x$value[i]), paper = mean(x$value_paper[i]),
               within_2se = sum(ok[i]))
  }))
}

## One row per simulation configuration: means over the cells of the paper's
## tables (in brackets the share of cells within two Monte Carlo standard
## errors) and the discovery rates of B1 (n = 300) and B4 (n = 600) at
## eta = 0.2, theta = 0.5 / 0.9.
config_table <- function(configs, dir = "results") {
  rows <- lapply(names(configs), function(lab) {
    v <- sim_summary(lab, dir)
    if (is.null(v)) return(NULL)
    s <- v$summary
    a <- agreement(s, c("fdr", "p_one", "p_all", "p_exact", "auc", "rel_bias"))
    cell <- function(m, pct = TRUE) {
      r <- a[a$metric == m, ]
      if (!nrow(r)) return("-")
      if (pct) sprintf("%.3f (%.0f%%)", r$package, 100 * r$within_2se / r$cells) else sprintf("%.2f", r$package)
    }
    disc <- function(m, n) {
      v <- vapply(c(0.5, 0.9), function(th) {
        r <- s[s$metric == "discovery" & s$eta == 0.2 & s$theta == th & s$model == m & s$n == n, ]
        if (nrow(r)) sprintf("%.2f", r$value) else "-"
      }, "")
      paste(v, collapse = " / ")
    }
    data.frame(Configuration = configs[[lab]],
               Trials = max(s$n_trials[s$metric %in% c("fdr", "discovery")]),
               FDR = cell("fdr"), p_one = cell("p_one"), p_all = cell("p_all"),
               p_exact = cell("p_exact"), AUC = cell("auc"), rel_bias = cell("rel_bias", FALSE),
               `B1 n=300` = disc("B1", 300), `B4 n=600` = disc("B4", 600), check.names = FALSE)
  })
  do.call(rbind, rows)
}
