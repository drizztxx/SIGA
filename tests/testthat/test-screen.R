test_that("LRT screening reproduces anova() for every endpoint", {
  set.seed(4)
  n <- 200
  d <- data.frame(trt = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n),
                  f = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$x2[1:5] <- NA
  d$time <- rexp(n, exp(-0.8 * d$trt * (d$x1 > 0)))
  d$status <- rbinom(n, 1, 0.85)
  d$yb <- rbinom(n, 1, plogis(1.2 * d$trt * d$x1))
  d$yc <- d$trt * d$x1 + rnorm(n)
  sc <- siga_screen(Surv(time, status) ~ trt | x1 + x2 + f, data = d, n_candidates = 3)$table
  for (v in c("x1", "x2", "f")) {
    dd <- d[!is.na(d[[v]]), ]
    a <- anova(survival::coxph(as.formula(paste("Surv(time, status) ~ trt +", v)), data = dd),
               survival::coxph(as.formula(paste("Surv(time, status) ~ trt *", v)), data = dd))
    expect_equal(sc$p_value[sc$variable == v], a[2, 4], tolerance = 1e-8)
  }
  sb <- siga_screen(yb ~ trt | x1 + x2, data = d[c("trt", "x1", "x2", "yb")])$table
  a <- anova(glm(yb ~ trt + x1, binomial, d), glm(yb ~ trt * x1, binomial, d), test = "Chisq")
  expect_equal(sb$p_value[sb$variable == "x1"], a[2, "Pr(>Chi)"], tolerance = 1e-8)
  scn <- siga_screen(yc ~ trt | x1 + x2, data = d[c("trt", "x1", "x2", "yc")])$table
  a <- anova(lm(yc ~ trt + x1, d), lm(yc ~ trt * x1, d))
  expect_equal(scn$p_value[scn$variable == "x1"], a[2, "Pr(>F)"], tolerance = 1e-8)
})

test_that("binary screening ranks a true predictor first (legacy bug 4)", {
  set.seed(1)
  n <- 400
  d <- data.frame(trt = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-0.5 + 0.3 * d$trt + 1.5 * d$trt * d$x1))
  s <- siga_screen(y ~ trt | ., data = d)
  expect_identical(s$table$variable[1], "x1")
  expect_lt(s$table$p_value[1], 1e-3)
  expect_length(s$selected, 2L)                     # k = ceiling(3 / 2) = 2
  expect_true(all(diff(s$table$p_value) >= 0))
})

test_that("the default number of candidates follows the paper", {
  expect_identical(SIGA:::default_n_candidates(5), 3L)
  expect_identical(SIGA:::default_n_candidates(20), 10L)
  expect_identical(SIGA:::default_n_candidates(21), 10L)
  expect_identical(SIGA:::default_n_candidates(100), 10L)
})

test_that("random forest screening runs when randomForestSRC is available", {
  skip_if_not_installed("randomForestSRC")
  skip_on_cran()
  d <- siga_simulate(200, model = "B1", n_noise = 0, seed = 8)
  s <- siga_screen(Surv(time, status) ~ trt | x1 + x2 + x4 + x5, data = d, method = "rf",
                   ntree = 100, seed = 1)
  expect_s3_class(s, "siga_screen")
  expect_setequal(s$table$variable, c("x1", "x2", "x4", "x5"))
})
