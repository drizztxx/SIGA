test_that("q7() reproduces stats::quantile(type = 7)", {
  set.seed(1)
  for (x in list(rnorm(37), c(1, 1, 2, 2, 2, 5), round(runif(50) * 4), 3)) {
    p <- c(0, 1e-9, 0.1, 0.25, 0.5, runif(20), 1)
    expect_identical(SIGA:::q7(sort(x), p), unname(quantile(x, p, type = 7)))
  }
})

test_that("mask keys are distinct for distinct masks", {
  M <- cbind(c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), c(FALSE, FALSE, TRUE))
  k <- SIGA:::mask_keys(M)
  expect_identical(k[1], k[2])
  expect_false(k[1] == k[3])
})

make_data <- function(n = 200) {
  set.seed(2)
  data.frame(a = rnorm(n), b = runif(n, 0, 100),
             f1 = factor(rbinom(n, 1, 0.5)), f2 = factor(rbinom(n, 1, 0.5)),
             m1 = factor(sample(letters[1:3], n, TRUE)), m2 = factor(sample(LETTERS[1:4], n, TRUE)),
             o = factor(sample(c("lo", "mid", "hi"), n, TRUE), levels = c("lo", "mid", "hi"),
                        ordered = TRUE))
}

test_that("a single continuous condition uses the quantile of its own variable (legacy bug 1)", {
  X <- make_data()
  tp <- SIGA:::type_predictors(X)
  enc <- SIGA:::encode_predictors(tp$X, names(tp$types), tp$types)
  chrom <- c(0, 3.5, 0, 0, 0, 0, 0)                       # only b, ">=" at its median
  cd <- SIGA:::decode_conditions(chrom, enc)
  expect_length(cd, 1L)
  expect_identical(cd[[1]]$op, ">=")
  expect_equal(cd[[1]]$value, unname(quantile(X$b, 0.5)))
  m <- SIGA:::decode_masks(matrix(chrom, 1), enc)[, 1]
  expect_identical(m, X$b >= quantile(X$b, 0.5))
})

test_that("a single binary condition selects the coded level (legacy bug 2)", {
  X <- make_data()
  tp <- SIGA:::type_predictors(X)
  enc <- SIGA:::encode_predictors(tp$X, names(tp$types), tp$types)
  m_hi <- SIGA:::decode_masks(matrix(c(0, 0, 3.5, 0, 0, 0, 0), 1), enc)[, 1]
  m_lo <- SIGA:::decode_masks(matrix(c(0, 0, 2.5, 0, 0, 0, 0), 1), enc)[, 1]
  expect_identical(m_hi, X$f1 == "1")
  expect_identical(m_lo, X$f1 == "0")
})

test_that("a single categorical condition decodes to a level subset (legacy bug 3)", {
  X <- make_data()
  tp <- SIGA:::type_predictors(X)
  enc <- SIGA:::encode_predictors(tp$X, names(tp$types), tp$types)
  ## m1 has 3 levels: 6 subsets a, b, c, ab, ac, bc; gene in (2, 4] -> block ceiling((g-2)*3)
  expected <- list("a", "b", "c", c("a", "b"), c("a", "c"), c("b", "c"))
  for (s in 1:6) {
    g <- 2 + (s - 0.5) / 3
    cd <- SIGA:::decode_conditions(c(0, 0, 0, 0, g, 0, 0), enc)
    expect_identical(cd[[1]]$value, expected[[s]])
    m <- SIGA:::decode_masks(matrix(c(0, 0, 0, 0, g, 0, 0), 1), enc)[, 1]
    expect_identical(m, as.character(X$m1) %in% expected[[s]])
  }
})

test_that("decode_masks() agrees with the decoded rule for random populations", {
  X <- make_data()
  X$a[c(3, 10)] <- NA
  X$f2[5] <- NA
  X$m2[7] <- NA
  tp <- SIGA:::type_predictors(X)
  enc <- SIGA:::encode_predictors(tp$X, names(tp$types), tp$types)
  set.seed(3)
  Chrom <- matrix(runif(60 * 7, 0, 4), 60)
  M <- SIGA:::decode_masks(Chrom, enc)
  for (i in seq_len(nrow(Chrom))) {
    rule <- SIGA:::new_siga_rule(SIGA:::decode_conditions(Chrom[i, ], enc))
    expect_identical(rule_members(rule, tp$X), M[, i])
  }
  ## missing predictor values are never in the subgroup when the predictor is used
  used_a <- which(Chrom[, 1] > 2)
  expect_false(any(M[c(3, 10), used_a]))
})

test_that("ordinal predictors use thresholds on the level order", {
  X <- make_data()
  tp <- SIGA:::type_predictors(X)
  expect_identical(unname(tp$types["o"]), "ordinal")
  enc <- SIGA:::encode_predictors(tp$X, names(tp$types), tp$types)
  rule <- SIGA:::new_siga_rule(SIGA:::decode_conditions(c(0, 0, 0, 0, 0, 0, 2.5), enc))
  expect_identical(rule_members(rule, X), as.integer(X$o) <= quantile(as.integer(X$o), 0.5))
})

test_that("predictor types are detected and can be overridden", {
  X <- data.frame(n = 1:10, l = rep(c(TRUE, FALSE), 5), ch = rep(c("u", "v", "w", "u", "v"), 2),
                  b01 = rep(0:1, 5), const = 1)
  expect_warning(tp <- SIGA:::type_predictors(X), "without variation")
  expect_identical(unname(tp$types), c("numeric", "binary", "categorical", "numeric"))
  tp2 <- suppressWarnings(SIGA:::type_predictors(X, c(b01 = "binary", n = "ordinal")))
  expect_identical(unname(tp2$types[c("n", "b01")]), c("ordinal", "binary"))
  expect_error(SIGA:::type_predictors(X["ch"], c(ch = "binary")), "cannot be treated as binary")
})
