legacy_exprs <- c(A1 = "0", A2 = "2*(x1+x2)", B1 = "0.1*(x1+x2)+2*z*(x4>1)",
                  B2 = "0.2*(x1+x2)+1.2*z*(x1>0)", B3 = "0.3*(x2+x9)+2*z*(x4>1 & x5==1)",
                  B4 = "2*z*(x1<0.5 & x5==1 & x6==1)")
legacy_pi <- c(A1 = 3.65, A2 = 4.42, B1 = 4.05, B2 = 3.97, B3 = 3.86, B4 = 3.84)

test_that("legacy = TRUE reproduces the original simulation program exactly", {
  for (m in names(legacy_exprs)) {
    old <- simul.Cox(n = 150, exp = legacy_exprs[[m]], weib.shape = 4.34, weib.scale = 1.2,
                     censor.par = legacy_pi[[m]], extra.noise = 10, seed = 11)
    new <- siga_simulate(150, model = m, legacy = TRUE, seed = 11)
    expect_identical(old$Ts, new$time, info = m)
    expect_identical(old$event, new$status, info = m)
    expect_identical(as.numeric(as.character(old$z)), as.numeric(new$trt), info = m)
    for (v in paste0("x", 1:20))
      expect_identical(as.character(old[[v]]), as.character(new[[v]]), info = paste(m, v))
  }
})

test_that("censoring and the true subgroups follow the paper", {
  for (m in c("A1", "B1", "B4")) {
    d <- siga_simulate(20000, model = m, n_noise = 0, seed = 12)
    expect_equal(1 - mean(d$status), 0.3, tolerance = 0.02, info = m)
    expect_true(all(d$time <= pmax(d$time, 0)))
  }
  d <- siga_simulate(500, model = "B3", seed = 13)
  expect_identical(attr(d, "true_subgroup"), d$x4 > 1 & d$x5 == "1")
  expect_identical(attr(d, "predictive"), c("x4", "x5"))
  expect_identical(ncol(d), 23L)
  expect_identical(siga_true_subgroup(d, "A2"), rep(FALSE, 500))
})

test_that("simulation with a seed is reproducible and does not change the RNG state", {
  set.seed(1)
  before <- .Random.seed
  a <- siga_simulate(50, model = "B2", seed = 3)
  expect_identical(.Random.seed, before)
  expect_identical(a, siga_simulate(50, model = "B2", seed = 3))
  expect_error(siga_simulate(50, model = "B2", legacy = TRUE), "seed")
  expect_error(siga_simulate(50, model = "B2", endpoint = "binary", legacy = TRUE, seed = 1), "time-to-event")
  b <- siga_simulate(50, model = "B2", endpoint = "binary", seed = 3)
  expect_true(all(b$y %in% 0:1))
})
