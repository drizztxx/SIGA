## Objective functions -------------------------------------------------------
##
## An objective maps a subgroup membership vector to a value that the genetic
## algorithm minimises. Built-in objectives are one-sided p-values oriented in
## the direction of treatment benefit:
##
##  "interaction"  (default; SIGA paper, section 2.3.2) p-value of the
##                 treatment-by-subgroup interaction gamma in
##                 Cox:      h(t) = h0(t) exp(alpha I + beta Z + gamma Z I)
##                 logistic: logit P(Y = 1) = mu + alpha I + beta Z + gamma Z I
##                 linear:   E(Y) = mu + alpha I + beta Z + gamma Z I (t test)
##  "subgroup"     p-value of the treatment effect within the subgroup; set to
##                 1 if the complement has a smaller p-value (legacy objF.cox,
##                 objF.glm2).
##  "differential" p-value of (Z_sub - Z_comp) / sqrt(2), where Z are the
##                 oriented Wald statistics of the treatment effect in the
##                 subgroup and its complement (legacy objF.cox2).
##
## Size rule (legacy implementation of MinS): the objective is 1 when
## 2 * min(n_treated, n_control) in the subgroup is less than
## min_size * N, or when the subgroup or its complement is empty.

objective_names <- c("interaction", "subgroup", "differential")

make_objective <- function(ep, objective = "interaction", min_size = 0.1) {
  z <- ep$z
  N <- ep$n
  custom <- is.function(objective)
  if (!custom) objective <- match.arg(objective, objective_names)
  penalty <- if (custom) Inf else 1
  arm_stat <- function(idx) {
    s <- coef_stat(ep_subset(ep, idx), matrix(z[idx], ncol = 1L), 1L)
    if (is.null(s)) NULL else c(ep$sign * s[1L] / s[2L], s[3L])
  }
  function(mask) {
    n1 <- sum(mask & z == 1)
    nsub <- sum(mask)
    n0 <- nsub - n1
    bal <- 2 * min(n1, n0)
    if (bal == 0 || bal * 100 / N < min_size * 100 || nsub == N) return(penalty)
    if (custom) {
      val <- tryCatch(objective(ep_response(ep), z, mask, ep$type), error = function(e) NA_real_)
      return(if (length(val) == 1L && is.finite(val)) as.numeric(val) else penalty)
    }
    switch(objective,
      interaction = {
        I <- as.numeric(mask)
        s <- coef_stat(ep, cbind(z, I, z * I), 3L)
        if (is.null(s)) 1 else p_upper(ep$sign * s[1L] / s[2L], s[3L])
      },
      subgroup = {
        s1 <- arm_stat(which(mask))
        s2 <- arm_stat(which(!mask))
        p1 <- if (is.null(s1)) 1 else p_upper(s1[1L], s1[2L])
        p2 <- if (is.null(s2)) 1 else p_upper(s2[1L], s2[2L])
        if (p2 < p1) 1 else p1
      },
      differential = {
        s1 <- arm_stat(which(mask))
        s2 <- arm_stat(which(!mask))
        if (is.null(s1) || is.null(s2)) 1 else
          stats::pnorm((s1[1L] - s2[1L]) / sqrt(2), lower.tail = FALSE)
      })
  }
}
