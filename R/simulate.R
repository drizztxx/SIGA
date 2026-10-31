## Simulation models of the SIGA paper (section 3.1, Table 2) ------------------

siga_models <- list(
  A1 = list(lp = function(d) rep(0, nrow(d)), predictive = character(), censor = 3.65),
  A2 = list(lp = function(d) 2 * (d$x1 + d$x2), predictive = character(), censor = 4.42),
  B1 = list(lp = function(d) 0.1 * (d$x1 + d$x2) + 2 * d$trt * (d$x4 > 1),
            predictive = "x4", censor = 4.05),
  B2 = list(lp = function(d) 0.2 * (d$x1 + d$x2) + 1.2 * d$trt * (d$x1 > 0),
            predictive = "x1", censor = 3.97),
  B3 = list(lp = function(d) 0.3 * (d$x2 + d$x9) + 2 * d$trt * (d$x4 > 1 & d$x5 == 1),
            predictive = c("x4", "x5"), censor = 3.86),
  B4 = list(lp = function(d) 2 * d$trt * (d$x1 < 0.5 & d$x5 == 1 & d$x6 == 1),
            predictive = c("x1", "x5", "x6"), censor = 3.84))

#' Simulate data from the SIGA simulation models
#'
#' Generates two-arm trial data from the six models of the SIGA paper (section
#' 3.1, Table 2): the null models A1 (no effect) and A2 (prognostic effects
#' only) and the models B1 to B4 with a predictive subgroup. Binary and
#' continuous analogues with the same linear predictor are also available.
#'
#' @details
#' Covariates: \eqn{x_1, x_7, \dots} standard normal, \eqn{(x_2, x_3)}
#' bivariate normal with correlation 0.5, \eqn{x_4} exponential with rate 1,
#' \eqn{x_5, x_6} Bernoulli(0.5) (returned as factors with levels `"0"`, `"1"`),
#' \eqn{x_7, \dots, x_{10}} equicorrelated normal (correlation 0.5) and
#' `n_noise` further noise variables in blocks of 10 with within-block
#' correlation `noise_rho` (0.5 in the original simulation program; the paper
#' describes them as independent). The treatment `trt` is Bernoulli(0.5).
#'
#' Linear predictors \eqn{L}:
#' \tabular{ll}{
#'   A1 \tab \eqn{0} \cr
#'   A2 \tab \eqn{2(x_1 + x_2)} \cr
#'   B1 \tab \eqn{0.1(x_1 + x_2) + 2\,trt\,I(x_4 > 1)} \cr
#'   B2 \tab \eqn{0.2(x_1 + x_2) + 1.2\,trt\,I(x_1 > 0)} \cr
#'   B3 \tab \eqn{0.3(x_2 + x_9) + 2\,trt\,I(x_4 > 1, x_5 = 1)} \cr
#'   B4 \tab \eqn{2\,trt\,I(x_1 < 0.5, x_5 = 1, x_6 = 1)}
#' }
#' Positive values of \eqn{L} mean benefit. Time-to-event: survival times
#' \eqn{T = 1.2\,(-\log(U)\,e^{L})^{1/4.34}} (Weibull baseline, hazard
#' proportional to \eqn{e^{-L}}, so the hazard ratio in the B1 subgroup is
#' \eqn{e^{-2}}), censoring times \eqn{C \sim U(0, \pi)} with the values of
#' \eqn{\pi} of the paper that give 30% censoring (A1 3.65, A2 4.42, B1 4.05,
#' B2 3.97, B3 3.86, B4 3.84), observed time \eqn{\min(T, C)}. Binary:
#' \eqn{P(y = 1) = \mathrm{logit}^{-1}(\beta_0 + L)}. Continuous:
#' \eqn{y = L + \varepsilon}, \eqn{\varepsilon \sim N(0, \sigma^2)}.
#'
#' `legacy = TRUE` reproduces the original simulation program of the paper
#' (`simul.Cox`) exactly for a given seed, including the re-seeding of the
#' random number generator for the noise variables and its use of the event
#' time \eqn{T} instead of \eqn{\min(T, C)} as the observed time (the event
#' indicator is \eqn{I(T < C)} in both cases).
#'
#' @param n number of patients.
#' @param model one of `"A1"`, `"A2"`, `"B1"`, `"B2"`, `"B3"`, `"B4"`.
#' @param endpoint `"tte"`, `"binary"` or `"continuous"`.
#' @param n_noise number of additional noise covariates (10 in the paper, for
#'   20 covariates in total).
#' @param noise_rho correlation of the noise covariates within blocks of 10.
#' @param censor_max upper limit \eqn{\pi} of the uniform censoring
#'   distribution; by default the value of the paper for `model` (30%
#'   censoring). `Inf` gives no censoring.
#' @param intercept intercept \eqn{\beta_0} of the binary model.
#' @param sigma residual standard deviation of the continuous model.
#' @param legacy reproduce the original simulation program (time-to-event only;
#'   requires `seed` and `n_noise` a multiple of 10).
#' @param seed optional integer seed.
#' @return A data frame with the outcome (`time` and `status`, or `y`), the
#'   treatment `trt` (0/1) and the covariates `x1`, `x2`, .... The attributes
#'   `model`, `endpoint`, `predictive` (the true predictive covariates) and
#'   `true_subgroup` (logical membership of the true subgroup) are set.
#' @seealso [siga_true_subgroup()]
#' @examples
#' dat <- siga_simulate(300, model = "B3", seed = 1)
#' head(dat)
#' mean(1 - dat$status)              # about 30% censoring
#' table(attr(dat, "true_subgroup"))
#' @export
siga_simulate <- function(n, model = c("A1", "A2", "B1", "B2", "B3", "B4"),
                          endpoint = c("tte", "binary", "continuous"), n_noise = 10,
                          noise_rho = 0.5, censor_max = NULL, intercept = stats::qlogis(0.3),
                          sigma = 2, legacy = FALSE, seed = NULL) {
  check_number(n, "n", lower = 10, integer = TRUE)
  model <- match.arg(model)
  endpoint <- match.arg(endpoint)
  check_number(n_noise, "n_noise", lower = 0, integer = TRUE)
  check_number(noise_rho, "noise_rho", lower = 0, upper = 0.99)
  spec <- siga_models[[model]]
  if (legacy) {
    if (endpoint != "tte") stop("'legacy = TRUE' is only available for time-to-event data.", call. = FALSE)
    if (is.null(seed)) stop("'legacy = TRUE' requires 'seed'.", call. = FALSE)
    if (n_noise %% 10 != 0) stop("'legacy = TRUE' requires 'n_noise' to be a multiple of 10.", call. = FALSE)
  }
  if (!is.null(seed)) {
    check_number(seed, "seed", integer = TRUE)
    state <- save_rng()
    on.exit(restore_rng(state), add = TRUE)
  }
  gen <- function() {
    if (!is.null(seed))
      set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    x1 <- stats::rnorm(n)
    x23 <- mvtnorm::rmvnorm(n, sigma = matrix(c(1, 0.5, 0.5, 1), ncol = 2L))
    x4 <- stats::rexp(n)
    x5 <- stats::rbinom(n, 1, 0.5)
    x6 <- stats::rbinom(n, 1, 0.5)
    s4 <- matrix(0.5, 4L, 4L)
    diag(s4) <- 1
    x710 <- mvtnorm::rmvnorm(n, sigma = s4)
    trt <- stats::rbinom(n, 1, 0.5)
    u <- stats::runif(n)
    d <- data.frame(trt = trt, x1 = x1, x2 = x23[, 1L], x3 = x23[, 2L], x4 = x4, x5 = x5, x6 = x6,
                    x7 = x710[, 1L], x8 = x710[, 2L], x9 = x710[, 3L], x10 = x710[, 4L])
    if (n_noise > 0) {
      blocks <- ceiling(n_noise / 10)
      noise <- do.call(cbind, lapply(seq_len(blocks), function(b) {
        if (legacy) set.seed(seed + b)
        size <- if (b < blocks || n_noise %% 10 == 0) 10L else n_noise %% 10
        sb <- matrix(noise_rho, size, size)
        diag(sb) <- 1
        xx <- mvtnorm::rmvnorm(n, sigma = sb)
        colnames(xx) <- paste0("x", 10L * b + seq_len(size))
        xx
      }))
      d <- cbind(d, noise)
    }
    L <- spec$lp(d)
    if (endpoint == "tte") {
      t_event <- 1.2 * (-log(u) * exp(L))^(1 / 4.34)
      cmax <- censor_max %||% spec$censor
      t_cens <- if (is.finite(cmax)) stats::runif(n, 0, cmax) else rep(Inf, n)
      status <- as.numeric(t_cens > t_event)
      time <- if (legacy) t_event else pmin(t_event, t_cens)
      out <- data.frame(time = time, status = status)
    } else if (endpoint == "binary") {
      out <- data.frame(y = stats::rbinom(n, 1, stats::plogis(intercept + L)))
    } else {
      out <- data.frame(y = L + stats::rnorm(n, sd = sigma))
    }
    d <- cbind(out, d)
    d$x5 <- factor(d$x5, levels = c(0, 1))
    d$x6 <- factor(d$x6, levels = c(0, 1))
    d
  }
  d <- gen()
  attr(d, "model") <- model
  attr(d, "endpoint") <- endpoint
  attr(d, "predictive") <- spec$predictive
  attr(d, "true_subgroup") <- siga_true_subgroup(d, model)
  d
}

#' True subgroup of the SIGA simulation models
#'
#' @param data a data frame with the covariates `x1`, `x4`, `x5` and `x6`, for
#'   example generated by [siga_simulate()].
#' @param model one of `"A1"`, `"A2"`, `"B1"`, `"B2"`, `"B3"`, `"B4"`.
#' @return A logical vector: membership of the true predictive subgroup
#'   (`x4 > 1` for B1, `x1 > 0` for B2, `x4 > 1 & x5 = 1` for B3,
#'   `x1 < 0.5 & x5 = 1 & x6 = 1` for B4, none for A1 and A2).
#' @examples
#' dat <- siga_simulate(100, model = "B4", seed = 1)
#' mean(siga_true_subgroup(dat, "B4"))
#' @export
siga_true_subgroup <- function(data, model = c("A1", "A2", "B1", "B2", "B3", "B4")) {
  model <- match.arg(model)
  is1 <- function(x) !is.na(x) & as.character(x) == "1"
  switch(model,
    A1 = , A2 = rep(FALSE, nrow(data)),
    B1 = data$x4 > 1,
    B2 = data$x1 > 0,
    B3 = data$x4 > 1 & is1(data$x5),
    B4 = data$x1 < 0.5 & is1(data$x5) & is1(data$x6))
}
