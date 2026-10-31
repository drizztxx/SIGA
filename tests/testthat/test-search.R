small <- siga_control(pop_size = 20, max_gen = 15, n_subpop = 2)

test_that("siga_search() finds a strong planted subgroup", {
  d <- siga_simulate(400, model = "B1", n_noise = 0, seed = 21)
  fit <- siga_search(Surv(time, status) ~ trt | x1 + x2 + x4, data = d, control = small, seed = 1)
  expect_s3_class(fit, "siga_search")
  expect_true(fit$valid)
  expect_identical(rule_variables(fit$rule), "x4")
  expect_identical(fit$rule$conditions[[1]]$op, ">=")
  expect_equal(fit$rule$conditions[[1]]$value, 1, tolerance = 0.25)
  expect_identical(predict(fit), rule_members(fit$rule, d))
  expect_lt(fit$effects$table$estimate[2], 0.4)
  expect_identical(dim(fit$history), c(16L, 2L))
  expect_true(all(diff(fit$history[, 1]) <= 0))       # elitism: never gets worse
})

test_that("results are reproducible and leave the user's RNG state untouched", {
  d <- siga_simulate(200, model = "B2", n_noise = 0, seed = 22)
  set.seed(99)
  before <- .Random.seed
  f1 <- siga_search(Surv(time, status) ~ trt | x1 + x2 + x3, data = d, control = small, seed = 5)
  expect_identical(.Random.seed, before)
  f2 <- siga_search(Surv(time, status) ~ trt | x1 + x2 + x3, data = d, control = small, seed = 5)
  expect_identical(f1$history, f2$history)
  expect_identical(format(f1$rule, digits = NULL), format(f2$rule, digits = NULL))
  set.seed(7)
  f3 <- siga_search(Surv(time, status) ~ trt | x1 + x2 + x3, data = d, control = small)
  set.seed(7)
  f4 <- siga_search(Surv(time, status) ~ trt | x1 + x2 + x3, data = d, control = small)
  expect_identical(f3$seed, f4$seed)
  expect_identical(f3$history, f4$history)
})

test_that("parallel and sequential runs give identical results", {
  skip_on_cran()
  skip_if_not(SIGA:::siga_installed(), "SIGA not installed")
  skip_if(parallel::detectCores() < 2L)
  d <- siga_simulate(200, model = "B1", n_noise = 0, seed = 23)
  f1 <- siga_search(Surv(time, status) ~ trt | x1 + x4, data = d, control = small, seed = 3)
  f2 <- siga_search(Surv(time, status) ~ trt | x1 + x4, data = d, control = small, seed = 3,
                    n_cores = 2)
  expect_identical(f1$history, f2$history)
  expect_identical(f1$candidates, f2$candidates)
})

test_that("binary and continuous endpoints are searched in the direction of benefit", {
  db <- siga_simulate(600, model = "B1", endpoint = "binary", n_noise = 0, seed = 24)
  fb <- siga_search(y ~ trt | x1 + x4, data = db, control = small, seed = 1)
  expect_identical(rule_variables(fb$rule), "x4")
  expect_gt(fb$effects$table$estimate[2], 2)
  dc <- siga_simulate(600, model = "B1", endpoint = "continuous", n_noise = 0, seed = 24)
  fc <- siga_search(y ~ trt | x1 + x4, data = dc, control = small, seed = 1)
  expect_identical(rule_variables(fc$rule), "x4")
  ## reversing the direction looks for subgroups where the treatment is worse
  fr <- siga_search(y ~ trt | x1 + x4, data = dc, control = small, seed = 1, direction = "lower")
  expect_lt(fr$effects$interaction$estimate, 0.5)
})

test_that("the objective cache stores each subgroup once, in vectors", {
  d <- siga_simulate(200, model = "B1", n_noise = 0, seed = 26)
  prep <- SIGA:::prepare_siga_data(Surv(time, status) ~ trt | x1 + x4, d, quiet = TRUE)
  enc <- SIGA:::encode_predictors(prep$X, c("x1", "x4"), prep$types)
  cache <- SIGA:::new_objective_cache()
  set.seed(1)
  r1 <- SIGA:::ga_run(prep$ep, prep$X, enc, small, cache = cache)
  ## keys are not stored as names (names become R symbols, which are never freed)
  expect_setequal(ls(cache), c("keys", "vals"))
  expect_identical(length(cache$keys), r1$n_eval)
  expect_false(anyDuplicated(cache$keys) > 0)
  ## a second run on the same data reuses the cache and gives the same result
  ## as a run with its own cache
  set.seed(2)
  r2 <- SIGA:::ga_run(prep$ep, prep$X, enc, small, cache = cache)
  set.seed(2)
  r3 <- SIGA:::ga_run(prep$ep, prep$X, enc, small)
  expect_identical(r2$history, r3$history)
  expect_lt(r2$n_eval, r3$n_eval)
})

test_that("a user-defined objective can be used", {
  d <- siga_simulate(300, model = "B1", endpoint = "continuous", n_noise = 0, seed = 25)
  obj <- function(y, trt, subgroup, endpoint) {
    -(mean(y[subgroup & trt == 1]) - mean(y[subgroup & trt == 0]))
  }
  f <- siga_search(y ~ trt | x1 + x4, data = d, seed = 2,
                   control = siga_control(pop_size = 20, max_gen = 10, n_subpop = 1, objective = obj))
  expect_true(f$valid)
  expect_true(is.finite(f$objective))
})
