tiny <- siga_control(pop_size = 20, max_gen = 15, n_subpop = 2)

test_that("siga() identifies a strong subgroup and reports the decision", {
  d <- siga_simulate(400, model = "B1", n_noise = 0, seed = 31)
  fit <- siga(Surv(time, status) ~ trt | ., data = d, n_resample = 8, control = tiny, seed = 4,
              verbose = FALSE)
  expect_s3_class(fit, "siga")
  expect_true(fit$subgroup)
  expect_true("x4" %in% rule_variables(fit$rule))
  expect_identical(nrow(fit$resamples), 8L)
  expect_true(all(fit$resamples$n_train == round(0.632 * 400)))
  ## defaults of the original code: 10 resamples, each with the sub-populations of the final search
  expect_identical(formals(siga)$n_resample, 10)
  expect_identical(fit$settings$resample_subpop, tiny$n_subpop)
  d1 <- fit$decision
  expect_equal(d1$threshold, (1 - 0.2) * fit$overall_ratio)
  expect_lt(d1$mean_ratio, d1$threshold)
  expect_equal(unname(fit$frequency[["x4"]]),
               mean(vapply(fit$resample_rules, function(r) "x4" %in% rule_variables(r), logical(1))))
  expect_output(print(fit), "Subgroup:")
  expect_output(print(summary(fit)), "selection frequencies")
  expect_s3_class(subgroup_effect(fit), "siga_effect")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  for (type in c("frequency", "resamples", "convergence", "effect")) expect_silent(plot(fit, type))
})

test_that("update() equals a fresh fit with the new tuning parameters", {
  d <- siga_simulate(300, model = "B3", n_noise = 0, seed = 32)
  f1 <- siga(Surv(time, status) ~ trt | ., data = d, n_resample = 6, control = tiny, seed = 9,
             verbose = FALSE)
  f2 <- siga(Surv(time, status) ~ trt | ., data = d, n_resample = 6, control = tiny, seed = 9,
             theta = 0.8, eta = 0.1, verbose = FALSE)
  u <- update(f1, theta = 0.8, eta = 0.1)
  expect_identical(u$decision, f2$decision)
  expect_identical(u$subgroup, f2$subgroup)
  if (f2$subgroup) expect_identical(format(u$rule, digits = NULL), format(f2$rule, digits = NULL))
  expect_identical(update(u, theta = f1$settings$theta, eta = f1$settings$eta)$decision, f1$decision)
})

test_that("siga() stops when the held-out effect does not improve", {
  d <- siga_simulate(300, model = "A1", n_noise = 0, seed = 33)
  fit <- siga(Surv(time, status) ~ trt | ., data = d, n_resample = 6, control = tiny, seed = 1,
              eta = 0.9, verbose = FALSE)
  expect_false(fit$subgroup)
  expect_null(fit$rule)
  expect_match(fit$decision$reason, "mean test-set effect ratio")
  expect_output(print(fit), "No subgroup")
  expect_warning(p <- predict(fit, d), "No subgroup")
  expect_false(any(p))
  expect_error(subgroup_effect(fit), "No subgroup")
})

test_that("siga() works for binary and continuous endpoints and bernoulli resampling", {
  db <- siga_simulate(400, model = "B1", endpoint = "binary", n_noise = 0, seed = 34)
  fb <- siga(y ~ trt | ., data = db, n_resample = 5, control = tiny, seed = 2, verbose = FALSE,
             resample_method = "bernoulli")
  expect_s3_class(fb, "siga")
  expect_false(all(fb$resamples$n_train == fb$resamples$n_train[1]))
  ## default eta: 0.25 for binary endpoints, 0.2 otherwise; an explicit value is kept
  expect_identical(fb$settings$eta, 0.25)
  expect_equal(fb$decision$threshold, (1 - 0.25) * fb$overall_ratio)
  expect_identical(update(fb, eta = 0.2)$settings$eta, 0.2)
  dc <- siga_simulate(400, model = "B1", endpoint = "continuous", n_noise = 0, seed = 34)
  fc <- siga(y ~ trt | x1 + x2 + x4 + x5, data = dc, n_resample = 5, control = tiny, seed = 2,
             verbose = FALSE, screen = "none")
  expect_identical(fc$settings$eta, 0.2)
  expect_true(fc$subgroup)
  expect_true("x4" %in% rule_variables(fc$rule))
})

test_that("input errors are informative", {
  d <- siga_simulate(100, model = "A1", n_noise = 0, seed = 35)
  expect_error(siga(Surv(time, status) ~ trt, data = d), "treatment \\| predictors")
  expect_error(siga(Surv(time, status) ~ trt + x1 | x2, data = d), "single variable")
  expect_error(siga(Surv(time, status) ~ trt | x1:x2, data = d), "plain column names")
  expect_error(siga(Surv(time, status) ~ nope | x1, data = d), "not found")
  expect_error(siga(Surv(time, status) ~ trt | x1, data = d, eta = 1), "eta")
  expect_error(siga(Surv(time, status) ~ trt | x1, data = d, theta = 0), "theta")
  expect_error(siga(Surv(time, status) ~ trt | x1, data = d, control = list()), "siga_control")
  d$arm <- sample(c("a", "b", "c"), nrow(d), TRUE)
  expect_error(siga(Surv(time, status) ~ arm | x1, data = d), "exactly two")
  d$arm2 <- ifelse(d$trt == 1, "active", "placebo")
  expect_error(siga_search(Surv(time, status) ~ arm2 | x1, data = d, control_arm = "nope"),
               "control_arm")
  expect_message(siga_search(Surv(time, status) ~ arm2 | x1 + x4, data = d, seed = 1,
                             control = siga_control(pop_size = 10, max_gen = 2, n_subpop = 1)),
                 "control arm")
})

test_that("siga_cv() returns one row per fold", {
  d <- siga_simulate(240, model = "B1", n_noise = 0, seed = 36)
  cv <- siga_cv(Surv(time, status) ~ trt | ., data = d, k = 3, n_resample = 4, control = tiny,
                seed = 1, verbose = FALSE)
  expect_s3_class(cv, "siga_cv")
  expect_identical(nrow(cv$results), 3L)
  expect_true(all(is.na(cv$results$error)))
  expect_output(print(summary(cv)), "Proportion of folds")
})
