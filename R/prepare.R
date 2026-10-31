## Formula parsing and data preparation -----------------------------------

## Split `outcome ~ treatment | predictors` into its parts. The predictor part
## may be `.` (all remaining columns), `. - x` or a sum of column names.
parse_siga_formula <- function(formula, data) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("'formula' must be a two-sided formula: outcome ~ treatment | predictors.",
         call. = FALSE)
  rhs <- formula[[3L]]
  if (!is.call(rhs) || !identical(rhs[[1L]], as.name("|")))
    stop("The right-hand side of 'formula' must have the form 'treatment | predictors', ",
         "e.g. Surv(time, status) ~ trt | x1 + x2 or y ~ trt | .", call. = FALSE)
  trt <- all.vars(rhs[[2L]])
  if (length(trt) != 1L || !is.name(rhs[[2L]]))
    stop("The treatment part of 'formula' must be a single variable name.", call. = FALSE)
  if (!trt %in% names(data))
    stop(sprintf("Treatment variable '%s' not found in 'data'.", trt), call. = FALSE)
  outcome_vars <- all.vars(formula[[2L]])
  miss <- setdiff(outcome_vars, names(data))
  if (length(miss))
    stop("Outcome variable(s) not found in 'data': ", paste(miss, collapse = ", "),
         call. = FALSE)
  pool <- setdiff(names(data), c(outcome_vars, trt))
  tt <- tryCatch(stats::terms(stats::as.formula(call("~", rhs[[3L]]), env = globalenv()),
                              data = data[pool]),
                 error = function(e) stop("Cannot interpret the predictor part of 'formula': ",
                                          conditionMessage(e), call. = FALSE))
  preds <- gsub("^`|`$", "", attr(tt, "term.labels"))
  bad <- setdiff(preds, pool)
  if (length(bad))
    stop("Predictors must be plain column names of 'data' other than the outcome and ",
         "treatment variables; problematic term(s): ", paste(bad, collapse = ", "),
         call. = FALSE)
  if (!length(preds)) stop("No predictors specified in 'formula'.", call. = FALSE)
  list(lhs = formula[[2L]], outcome_vars = outcome_vars, treatment = trt,
       predictors = preds, env = environment(formula) %||% globalenv())
}

## Detect the endpoint type and build the outcome part of the internal
## endpoint object.
build_outcome <- function(y, endpoint = NULL) {
  if (inherits(y, "Surv")) {
    if (!identical(attr(y, "type"), "right"))
      stop("Only right-censored survival outcomes are supported.", call. = FALSE)
    if (!is.null(endpoint) && endpoint != "tte")
      stop("A Surv() outcome requires endpoint = \"tte\".", call. = FALSE)
    return(list(type = "tte", time = as.numeric(y[, 1L]), status = as.numeric(y[, 2L]),
                event_level = NULL))
  }
  if (!is.null(endpoint) && endpoint == "tte")
    stop("endpoint = \"tte\" requires a Surv(time, status) outcome.", call. = FALSE)
  if (is.matrix(y) || is.data.frame(y))
    stop("The outcome must be a single variable or a Surv() object.", call. = FALSE)
  vals <- unique(y[!is.na(y)])
  is_two <- length(vals) == 2L
  if (is.logical(y) || is.factor(y) || is.character(y)) {
    if (!is.null(endpoint) && endpoint == "continuous")
      stop("A continuous endpoint requires a numeric outcome.", call. = FALSE)
    if (!is_two)
      stop("A binary outcome must have exactly two distinct non-missing values.", call. = FALSE)
    lev <- if (is.logical(y)) c(FALSE, TRUE) else if (is.factor(y)) levels(droplevels(y)) else
      sort(as.character(vals))
    event <- lev[2L]
    return(list(type = "binary", y = as.numeric(as.character(y) == as.character(event)),
                event_level = as.character(event)))
  }
  if (!is.numeric(y)) stop("Unsupported outcome type.", call. = FALSE)
  binary01 <- is_two && all(sort(vals) == c(0, 1))
  type <- endpoint %||% if (binary01) "binary" else "continuous"
  if (type == "binary") {
    if (!binary01) stop("A numeric binary outcome must be coded 0/1.", call. = FALSE)
    return(list(type = "binary", y = as.numeric(y), event_level = "1"))
  }
  list(type = "continuous", y = as.numeric(y), event_level = NULL)
}

## Treatment indicator z (1 = experimental, 0 = control).
build_treatment <- function(trt, control_arm = NULL, quiet = FALSE) {
  vals <- unique(trt[!is.na(trt)])
  if (length(vals) != 2L)
    stop(sprintf("The treatment variable must have exactly two distinct non-missing values (found %d).",
                 length(vals)), call. = FALSE)
  if (is.null(control_arm)) {
    if (is.logical(trt)) {
      control_arm <- FALSE
    } else if (is.numeric(trt) && all(sort(vals) == c(0, 1))) {
      control_arm <- 0
    } else {
      control_arm <- if (is.factor(trt)) levels(droplevels(trt))[1L] else sort(as.character(vals))[1L]
      if (!quiet)
        message(sprintf("Using '%s' as the control arm (set 'control_arm' to change).", control_arm))
    }
  }
  if (length(control_arm) != 1L || !as.character(control_arm) %in% as.character(vals))
    stop("'control_arm' must be one of the two values of the treatment variable: ",
         paste(sort(as.character(vals)), collapse = ", "), call. = FALSE)
  z <- as.numeric(as.character(trt) != as.character(control_arm))
  list(z = z, control_arm = as.character(control_arm),
       treatment_arm = setdiff(as.character(vals), as.character(control_arm)))
}

## Resolve the type of every predictor: "numeric" (threshold rules),
## "ordinal" (thresholds on the level order of an ordered factor), "binary"
## (one of two levels) or "categorical" (a subset of > 2 levels).
type_predictors <- function(X, predictor_types = NULL) {
  allowed <- c("numeric", "ordinal", "binary", "categorical")
  if (!is.null(predictor_types)) {
    if (is.null(names(predictor_types)) || !all(unlist(predictor_types) %in% allowed))
      stop("'predictor_types' must be a named vector with values in ",
           paste(allowed, collapse = ", "), ".", call. = FALSE)
    unknown <- setdiff(names(predictor_types), names(X))
    if (length(unknown))
      stop("'predictor_types' refers to unknown predictors: ", paste(unknown, collapse = ", "),
           call. = FALSE)
  }
  types <- setNames(character(ncol(X)), names(X))
  for (v in names(X)) {
    x <- X[[v]]
    type <- if (!is.null(predictor_types) && v %in% names(predictor_types)) {
      predictor_types[[v]]
    } else if (is.ordered(x)) {
      "ordinal"
    } else if (is.factor(x) || is.character(x) || is.logical(x)) {
      "factor"
    } else if (is.numeric(x)) {
      "numeric"
    } else {
      stop(sprintf("Predictor '%s' has unsupported class '%s'.", v, class(x)[1L]), call. = FALSE)
    }
    if (type == "numeric") {
      if (is.factor(x)) x <- if (is.ordered(x)) as.integer(x) else as.character(x)
      xn <- suppressWarnings(as.numeric(x))
      if (any(is.na(xn) & !is.na(x)))
        stop(sprintf("Predictor '%s' cannot be converted to numeric.", v), call. = FALSE)
      x <- xn
      usable <- length(unique(x[!is.na(x)])) >= 2L
    } else if (type == "ordinal") {
      x <- if (is.ordered(x)) droplevels(x) else
        factor(x, levels = sort(unique(x[!is.na(x)])), ordered = TRUE)
      usable <- nlevels(x) >= 2L
    } else {
      x <- if (is.logical(x)) droplevels(factor(x, levels = c(FALSE, TRUE))) else
        droplevels(as.factor(x))
      nl <- nlevels(x)
      if (type == "factor") type <- if (nl > 2L) "categorical" else "binary"
      if (type == "binary" && nl > 2L)
        stop(sprintf("Predictor '%s' has %d levels and cannot be treated as binary.", v, nl),
             call. = FALSE)
      if (type == "categorical" && nl == 2L) type <- "binary"
      usable <- nl >= 2L
    }
    X[[v]] <- x
    types[v] <- if (usable) type else "constant"
  }
  const <- names(types)[types == "constant"]
  if (length(const)) {
    warning("Dropping predictor(s) without variation: ", paste(const, collapse = ", "),
            call. = FALSE)
    X <- X[setdiff(names(X), const)]
    types <- types[setdiff(names(types), const)]
  }
  if (!length(types)) stop("No usable predictors.", call. = FALSE)
  list(X = X, types = types)
}

## Build the internal analysis object (endpoint + typed predictors).
prepare_siga_data <- function(formula, data, control_arm = NULL, endpoint = NULL,
                              direction = NULL, predictor_types = NULL, quiet = FALSE) {
  pe <- prepare_endpoint(formula, data, control_arm, endpoint, direction, quiet)
  tp <- type_predictors(as.data.frame(data)[pe$rows, pe$predictors, drop = FALSE], predictor_types)
  rownames(tp$X) <- NULL
  list(ep = pe$ep, X = tp$X, types = tp$types, info = pe$info, rows = pe$rows)
}

## Endpoint part only: outcome, treatment indicator and orientation.
prepare_endpoint <- function(formula, data, control_arm = NULL, endpoint = NULL,
                             direction = NULL, quiet = FALSE) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  data <- as.data.frame(data)
  if (!is.null(endpoint)) endpoint <- match.arg(endpoint, c("tte", "binary", "continuous"))
  pf <- parse_siga_formula(formula, data)
  env <- list2env(list(Surv = survival::Surv), parent = pf$env)
  y <- eval(pf$lhs, data, env)
  if (NROW(y) != nrow(data))
    stop("The outcome must have one value per row of 'data'.", call. = FALSE)
  out <- build_outcome(y, endpoint)
  trt_raw <- data[[pf$treatment]]

  ok <- !is.na(trt_raw) & if (out$type == "tte") !is.na(out$time) & !is.na(out$status) else !is.na(out$y)
  if (!all(ok) && !quiet)
    message(sprintf("Removing %d row(s) with a missing outcome or treatment.", sum(!ok)))
  if (sum(ok) < 10L) stop("Too few complete observations.", call. = FALSE)
  trt <- build_treatment(trt_raw[ok], control_arm, quiet = quiet)

  if (is.null(direction)) {
    direction <- if (out$type == "tte") "lower" else "higher"
  } else {
    direction <- match.arg(direction, c("lower", "higher"))
  }
  ep <- list(type = out$type, z = trt$z, n = sum(ok), direction = direction,
             sign = if (direction == "higher") 1 else -1)
  if (out$type == "tte") {
    ep$time <- out$time[ok]
    ep$status <- out$status[ok]
    if (any(ep$time < 0)) stop("Survival times must be non-negative.", call. = FALSE)
    if (!all(ep$status %in% c(0, 1))) stop("The event indicator must be coded 0/1.", call. = FALSE)
    ep$surv <- survival::Surv(ep$time, ep$status)
    ep$cox_ctl <- survival::coxph.control()
  } else {
    ep$y <- out$y[ok]
  }
  if (out$type == "binary") ep$family <- stats::binomial()
  if (out$type == "continuous") {
    y1 <- ep$y[ep$z == 1]; y0 <- ep$y[ep$z == 0]
    ep$sd_y <- sqrt(((length(y1) - 1) * stats::var(y1) + (length(y0) - 1) * stats::var(y0)) /
                      (length(y1) + length(y0) - 2))
    if (!is.finite(ep$sd_y) || ep$sd_y <= 0)
      stop("The continuous outcome has no within-arm variation.", call. = FALSE)
  }
  if (min(sum(ep$z == 1), sum(ep$z == 0)) < 2L)
    stop("Each treatment arm needs at least two observations.", call. = FALSE)
  info <- list(outcome = deparse1(pf$lhs), outcome_vars = pf$outcome_vars,
               treatment = pf$treatment, control_arm = trt$control_arm,
               treatment_arm = trt$treatment_arm, event_level = out$event_level,
               endpoint = out$type, direction = direction)
  list(ep = ep, info = info, rows = which(ok), predictors = pf$predictors)
}

## Subset the endpoint object to the observations `idx`.
ep_subset <- function(ep, idx) {
  ep$z <- ep$z[idx]
  ep$n <- length(ep$z)
  if (ep$type == "tte") {
    ep$time <- ep$time[idx]
    ep$status <- ep$status[idx]
    ep$surv <- ep$surv[idx]
  } else {
    ep$y <- ep$y[idx]
  }
  ep
}

## Response object passed to user-defined objective functions.
ep_response <- function(ep) if (ep$type == "tte") ep$surv else ep$y
