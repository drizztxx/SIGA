rule_data <- function() {
  set.seed(5)
  n <- 100
  data.frame(age = round(rnorm(n, 60, 10)), sex = factor(sample(c("F", "M"), n, TRUE)),
             race = factor(sample(c("A", "B", "C"), n, TRUE)),
             stage = factor(sample(c("I", "II", "III"), n, TRUE), levels = c("I", "II", "III"),
                            ordered = TRUE),
             flag = sample(c(TRUE, FALSE), n, TRUE))
}

test_that("rules format as text and as evaluable R code", {
  X <- rule_data()
  rule <- SIGA:::new_siga_rule(list(
    list(var = "age", type = "numeric", op = ">=", value = 61.23456),
    list(var = "sex", type = "binary", op = "==", value = "M", levels = c("F", "M")),
    list(var = "race", type = "categorical", op = "in", value = c("A", "C"), levels = c("A", "B", "C")),
    list(var = "stage", type = "ordinal", op = "<=", value = 2, levels = c("I", "II", "III")),
    list(var = "flag", type = "binary", op = "==", value = "TRUE", levels = c("FALSE", "TRUE"))))
  expect_identical(format(rule),
                   "age >= 61.23 & sex = M & race in {A, C} & stage <= II & flag = TRUE")
  r_code <- format(rule, style = "R", digits = NULL)
  m_eval <- with(X, eval(parse(text = r_code)))
  expect_identical(rule_members(rule, X), m_eval & !is.na(m_eval))
  expect_identical(rule_variables(rule), c("age", "sex", "race", "stage", "flag"))
  expect_identical(format(SIGA:::new_siga_rule()), "<all patients>")
})

test_that("canonical thresholds and simplification keep the membership", {
  X <- rule_data()
  conds <- list(
    list(var = "age", type = "numeric", op = "<=", value = 63.7),
    list(var = "age", type = "numeric", op = ">=", value = -Inf),  # redundant
    list(var = "race", type = "categorical", op = "in", value = c("A", "B", "C"),
         levels = c("A", "B", "C")))                              # redundant
  rule <- SIGA:::new_siga_rule(conds)
  mask <- rule_members(rule, X)
  canon <- SIGA:::canonicalize_conditions(conds, X, mask)
  expect_identical(canon[[1]]$value, max(X$age[mask]))
  expect_identical(rule_members(SIGA:::new_siga_rule(canon), X), mask)
  simp <- SIGA:::simplify_conditions(canon, X, mask)
  expect_length(simp, 1L)
  expect_identical(rule_members(SIGA:::new_siga_rule(simp), X), mask)
})

test_that("rule evaluation handles missing values and unknown variables", {
  X <- rule_data()
  X$age[1] <- NA
  rule <- SIGA:::new_siga_rule(list(list(var = "age", type = "numeric", op = "<=", value = Inf)))
  m <- rule_members(rule, X)
  expect_false(m[1])
  expect_true(all(m[-1]))
  bad <- SIGA:::new_siga_rule(list(list(var = "nope", type = "numeric", op = "<=", value = 1)))
  expect_error(rule_members(bad, X), "not found")
})
