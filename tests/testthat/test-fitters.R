sim_ep <- function(type, n = 150, seed = 1) {
  set.seed(seed)
  d <- data.frame(trt = rbinom(n, 1, 0.5), x = rnorm(n), g = rbinom(n, 1, 0.4))
  if (type == "tte") {
    d$time <- rexp(n, exp(-0.5 * d$trt * d$g + 0.2 * d$x))
    d$status <- rbinom(n, 1, 0.8)
    f <- Surv(time, status) ~ trt | x + g
  } else if (type == "binary") {
    d$y <- rbinom(n, 1, plogis(-0.3 + 0.8 * d$trt * d$g + 0.3 * d$x))
    f <- y ~ trt | x + g
  } else {
    d$y <- 1 + d$trt * d$g + 0.5 * d$x + rnorm(n)
    f <- y ~ trt | x + g
  }
  list(data = d, prep = SIGA:::prepare_siga_data(f, d))
}

test_that("fast fitters match coxph(), glm() and lm()", {
  for (type in c("tte", "binary", "continuous")) {
    s <- sim_ep(type)
    d <- s$data
    I <- as.numeric(d$g == 1)
    X <- cbind(d$trt, I, d$trt * I)
    st <- SIGA:::coef_stat(s$prep$ep, X, 3L)
    ref <- switch(type,
      tte = summary(survival::coxph(Surv(time, status) ~ trt * I, data = cbind(d, I = I)))$coefficients[3, c(1, 3)],
      binary = summary(glm(y ~ trt * I, family = binomial, data = cbind(d, I = I)))$coefficients[4, 1:2],
      continuous = summary(lm(y ~ trt * I, data = cbind(d, I = I)))$coefficients[4, 1:2])
    expect_equal(unname(st[1:2]), unname(ref), tolerance = 1e-8, info = type)
  }
})

test_that("non-estimable coefficients return NULL", {
  s <- sim_ep("tte")
  d <- s$data
  X <- cbind(d$trt, d$g, d$g)                 # aliased column
  expect_null(SIGA:::coef_stat(s$prep$ep, X, 3L))
})

test_that("the interaction objective is the one-sided interaction p-value", {
  for (type in c("tte", "binary", "continuous")) {
    s <- sim_ep(type)
    d <- s$data
    obj <- SIGA:::make_objective(s$prep$ep, "interaction", min_size = 0.1)
    mask <- d$g == 1
    I <- as.numeric(mask)
    fit <- switch(type,
      tte = summary(survival::coxph(Surv(time, status) ~ trt * I, data = cbind(d, I = I)))$coefficients[3, ],
      binary = summary(glm(y ~ trt * I, family = binomial, data = cbind(d, I = I)))$coefficients[4, ],
      continuous = summary(lm(y ~ trt * I, data = cbind(d, I = I)))$coefficients[4, ])
    stat <- if (type == "tte") -fit[["z"]] else fit[[3]]
    p <- if (type == "continuous") pt(stat, nrow(d) - 4, lower.tail = FALSE) else
      pnorm(stat, lower.tail = FALSE)
    expect_equal(obj(mask), p, tolerance = 1e-8, info = type)
  }
})

test_that("small subgroups and empty complements get the worst objective", {
  s <- sim_ep("tte", n = 100)
  obj <- SIGA:::make_objective(s$prep$ep, "interaction", min_size = 0.1)
  z <- s$prep$ep$z
  tiny <- rep(FALSE, 100)
  tiny[which(z == 1)[1:4]] <- TRUE            # only treated patients
  expect_identical(obj(tiny), 1)
  small <- rep(FALSE, 100)
  small[c(which(z == 1)[1:4], which(z == 0)[1:4])] <- TRUE   # 2*min = 8 < 10
  expect_identical(obj(small), 1)
  expect_identical(obj(rep(TRUE, 100)), 1)
  ok <- rep(FALSE, 100)
  ok[c(which(z == 1)[1:5], which(z == 0)[1:5])] <- TRUE      # 2*min = 10 = 10% of 100
  expect_lt(obj(ok), 1)
})

test_that("alternative and custom objectives work", {
  s <- sim_ep("continuous")
  mask <- s$data$g == 1
  for (o in c("subgroup", "differential")) {
    v <- SIGA:::make_objective(s$prep$ep, o)(mask)
    expect_true(v >= 0 && v <= 1)
  }
  custom <- function(y, trt, subgroup, endpoint) -mean(y[subgroup & trt == 1])
  expect_equal(SIGA:::make_objective(s$prep$ep, custom)(mask),
               -mean(s$data$y[mask & s$data$trt == 1]))
  bad <- function(y, trt, subgroup, endpoint) stop("boom")
  expect_identical(SIGA:::make_objective(s$prep$ep, bad)(mask), Inf)
})

test_that("effect ratios follow the documented conventions", {
  s <- sim_ep("tte")
  ep <- s$prep$ep
  hr <- exp(coef(survival::coxph(Surv(time, status) ~ trt, data = s$data)))
  expect_equal(exp(-SIGA:::effect_theta(ep, rep(TRUE, ep$n))), unname(hr), tolerance = 1e-10)
  expect_identical(SIGA:::test_ratio(ep, rep(TRUE, ep$n)), 1)
  expect_identical(SIGA:::test_ratio(ep, rep(FALSE, ep$n)), 1)
  sb <- sim_ep("binary")
  tab <- table(sb$data$trt, sb$data$y)
  or <- (tab["1", "1"] * tab["0", "0"]) / (tab["1", "0"] * tab["0", "1"])
  expect_equal(exp(SIGA:::effect_theta(sb$prep$ep, rep(TRUE, sb$prep$ep$n))), or)
  sc <- sim_ep("continuous")
  d <- sc$data
  delta <- mean(d$y[d$trt == 1]) - mean(d$y[d$trt == 0])
  expect_equal(SIGA:::effect_theta(sc$prep$ep, rep(TRUE, nrow(d))),
               pi / sqrt(3) * delta / sc$prep$ep$sd_y)
})
