process_focal_and_estimand_w_targets <- function(focal, estimand, targets = NULL, treat) {

  if (is_null(estimand)) {
    return(list(focal = NULL,
                estimand = NULL,
                reported.estimand = "targets"))
  }

  if (is_not_null(targets)) {
    arg::wrn("{.arg targets} are not {.val NULL}; ignoring {.arg estimand}")

    return(list(focal = NULL,
                estimand = NULL,
                reported.estimand = "targets"))
  }

  arg::arg_string(estimand)
  estimand <- toupper(estimand)

  if (!has_treat_type(treat)) treat <- assign_treat_type(treat)
  treat.type <- get_treat_type(treat)

  treat.type[treat.type == "multinomial"] <- "multi-category"

  reported.estimand <- estimand

  AE <- list(binary =  c("ATT", "ATC", "ATE"),
             `multi-category` = c("ATT", "ATC", "ATE"),
             continuous = "ATE")

  if (estimand %nin% AE[[treat.type]]) {
    arg::err("{.val {estimand}} is not an allowable estimand with {treat.type} treatments. Only {.val {AE[[treat.type]]}} {?is/are} allowed")
  }

  if (treat.type == "continuous") {
    return(list(focal = NULL,
                estimand = "ATE",
                reported.estimand = "ATE"))
  }

  if (is_not_null(focal) && estimand %nin% c("ATT", "ATC")) {
    arg::wrn('{.code estimand = "{estimand}"} is not compatible with {.arg focal}. Setting {.arg estimand} to {.val {"ATT"}}')

    reported.estimand <- estimand <- "ATT"
  }

  if (treat.type == "multi-category") {
    unique.vals <- {
      if (is.factor(treat) || is.character(treat))
        levels(factor(treat))
      else
        sort(unique(treat))
    }

    #Check focal
    if (estimand %in% c("ATT", "ATC")) {
      if (is_null(focal)) {
        arg::err('when {.code estimand = "{estimand}"} for multi-category treatments, an argument must be supplied to {.arg focal}')
      }

      if (length(focal) > 1L || focal %nin% unique.vals) {
        arg::err("{.arg focal} must be the name of a level of treatment ({.or {.val unique.vals}})")
      }
    }
    else {
      focal <- NULL
    }

    return(list(focal = unname(focal),
                estimand = estimand,
                reported.estimand = reported.estimand))
  }

  ct <- .get_control_and_treated_levels(treat, estimand, focal)

  focal <- switch(estimand,
                  ATT = ct["treated"],
                  ATC = ct["control"],
                  NULL)

  treated <- ct["treated"]

  list(focal = unname(focal),
       estimand = estimand,
       reported.estimand = reported.estimand,
       treated = unname(treated))
}

.get_control_and_treated_levels <- function(treat, estimand, focal = NULL, treated = NULL) {

  if (is_not_null(.attr(treat, "control")) &&
      is_not_null(.attr(treat, "treated"))) {
    return(setNames(c(.attr(treat, "control"), .attr(treat, "treated")),
                    c("control", "treated")))
  }

  control <- NULL
  throw_message <- FALSE

  unique.vals <- {
    if (is.factor(treat) || is.character(treat))
      levels(factor(treat, nmax = 2L))
    else
      sort(unique(treat, nmax = 2L))
  }

  if (is_not_null(focal)) {
    if (length(focal) > 1L || focal %nin% unique.vals) {
      arg::err("{.arg focal} must be the name of a level of treatment ({.or {.val unique.vals}})")
    }

    if (estimand == "ATC") {
      control <- focal
      treated <- NULL
    }
    else {
      treated <- focal
    }
  }
  else if (is_not_null(.attr(treat, "treated", TRUE))) {
    treated <- .attr(treat, "treated", TRUE)
  }
  else if (is_not_null(.attr(treat, "control", TRUE))) {
    control <- .attr(treat, "control", TRUE)
  }
  else if (is_not_null(treated)) {
    if (length(treated) > 1L || treated %nin% unique.vals) {
      arg::err("{.arg treated} must be the name of a level of treatment ({.or {.val unique.vals}})")
    }
  }
  else if (is.logical(treat)) {
    treated <- TRUE
    control <- FALSE
  }
  else if (is.numeric(unique.vals)) {
    control <- unique.vals[unique.vals == 0]

    if (is_null(control)) {
      control <- unique.vals[1L]

      throw_message <- TRUE
    }
  }
  else if (can_str2num(unique.vals)) {
    unique.vals.numeric <- str2num(unique.vals)

    control <- unique.vals[unique.vals.numeric == 0]

    if (is_null(control)) {
      control <- unique.vals[which.min(unique.vals.numeric)]

      throw_message <- TRUE
    }
  }
  else {
    treated_options <- c("t", "tr", "treat", "treated", "exposed")
    control_options <- c("c", "co", "ctrl", "control", "unexposed")

    t_match <- which(unique.vals %in% treated_options)
    c_match <- which(unique.vals %in% control_options)

    if (length(t_match) == 1L) {
      treated <- unique.vals[t_match]
    }
    else if (length(c_match) == 1L) {
      control <- unique.vals[c_match]
    }
  }

  if (is_null(control) && is_null(treated)) {
    control <- unique.vals[1L]
    treated <- unique.vals[2L]

    throw_message <- TRUE
  }
  else if (is_null(control)) {
    control <- setdiff(unique.vals, treated)
  }
  else if (is_null(treated)) {
    treated <- setdiff(unique.vals, control)
  }

  if (throw_message) {
    if (estimand == "ATT") {
      tl <- add_quotes(treated, !is.numeric(unique.vals))
      arg::msg("assuming {.val {tl}} is the treated level. If not, supply an argument to {.arg focal}")
    }
    else if (estimand == "ATC") {
      cl <- add_quotes(control, !is.numeric(unique.vals))
      arg::msg("assuming {.val {cl}} is the control level. If not, supply an argument to {.arg focal}")
    }
    else {
      tl <- add_quotes(treated, !is.numeric(unique.vals))
      arg::msg("assuming {.val {tl}} is the treated level. If not, recode the treatment so that 1 is treated and 0 is control")
    }
  }

  setNames(c(control, treated),
           c("control", "treated"))
}

get_treated_level <- function(treat, estimand, focal = NULL) {
  ct <- .get_control_and_treated_levels(treat, estimand, focal)

  unname(ct["treated"])
}

process_s.weights <- function(s.weights, data = NULL) {
  #Process s.weights
  if (is_null(s.weights)) {
    return(NULL)
  }

  if (is.numeric(s.weights)) {
    return(as.numeric(s.weights))
  }

  if (!rlang::is_string(s.weights)) {
    arg::err("{.arg s.weights} must be a vector or data frame of sampling weights or a string corresponding to the name of the variable in {.arg data} that contains sampling weights")
  }

  if (is_null(data)) {
    arg::err("{.arg s.weights} was specified as a string but there was no argument to {.arg data}")
  }

  if (s.weights %nin% names(data)) {
    arg::err("the value supplied to {.arg s.weights} is not the name of a variable in {.arg data}")
  }

  as.numeric(data[[s.weights]])
}
process_b.weights <- function(b.weights, data = NULL) {
  #Process b.weights
  if (is_null(b.weights)) {
    return(NULL)
  }

  if (is.numeric(b.weights)) {
    return(as.numeric(b.weights))
  }

  if (!rlang::is_string(b.weights)) {
    arg::err("the argument to {.arg b.weights} must be a vector or data frame of base weights or a string corresponding to the name of the variable in {.arg data} that contains base weights")
  }

  if (is_null(data)) {
    arg::err("{.arg b.weights} was specified as a string but there was no argument to {.arg data}")
  }

  if (b.weights %nin% names(data)) {
    arg::err("the value supplied to {.arg b.weights} is not the name of a variable in {.arg data}")
  }

  as.numeric(data[[b.weights]])
}
process_norm <- function(norm, s.weights, b.weights) {
  norm <- arg::match_arg(norm, allowable_norms())

  if (norm == "linf" && !all_the_same(s.weights)) {
    arg::err("{.arg norm} cannot be {.val {norm}} when sampling weights are used")
  }

  if (norm %in% c("entropy", "log") && any(b.weights <= 0)) {
    arg::err('all base weights must be positive when {.code norm = "{norm}"}')
  }

  norm
}
process_min.w <- function(min.w, norm, b.weights) {
  arg::arg_number(min.w)
  arg::arg_lte(min.w, mean(b.weights))

  if (norm %in% c("entropy", "log")) {
    min.w <- max(min.w, .Machine$double.eps)
  }

  min.w
}

check_missing_covs <- function(covs) {
  k <- ncol(covs)
  for (i in seq_len(k)) {
    if (anyNA(covs[[i]]) || (is.numeric(covs[[i]]) && !all(is.finite(covs[[i]])))) {
      covariates.with.missingness <- names(covs)[i:k][vapply(i:k, function(j) anyNA(covs[[j]]) ||
                                                               (is.numeric(covs[[j]]) && !all(is.finite(covs[[j]]))),
                                                             logical(1L))]
      arg::err("missing and non-finite values are not allowed in the covariates. Covariates with missingness or non-finite values: {.var covariates.with.missingness}")
    }
  }
}

process_solver <- function(solver = NULL, norm, min.w) {
  allowable_solvers <- {
    if (norm %in% c("entropy", "log"))
      c("scs", "clarabel")
    else if (norm %in% c("l1", "linf"))
      c("highs", "osqp", "lpsolve"[min.w >= 0])
    else
      c("osqp", "highs")
  }

  solver_option <- sprintf("optweight_solver_%s", norm)

  if (is_null(solver)) {
    solver <- getOption(solver_option)
  }

  if (is_null(solver)) {
    solver_packages <- c(scs = "scs",
                         clarabel = "clarabel",
                         osqp = "osqp",
                         highs = "highs",
                         lpsolve = "lpSolve")

    for (i in allowable_solvers) {
      if (rlang::is_installed(solver_packages[i])) {
        return(i)
      }
    }

    return(allowable_solvers[1L])
  }

  arg::match_arg(solver, allowable_solvers)
}

process_duals <- function(d, tols) {
  for (time in seq_len(max(d$component))) {
    .i <- d$component == time & d$constraint %in% c("target", "balance")

    original.vars <- {
      if (is.list(tols)) .attr(tols[[time]], "original.vars")
      else .attr(tols, "original.vars")
    }

    na.cov <- is.na(d$cov[.i])

    d$cov[.i][!na.cov] <- original.vars[d$cov[.i][!na.cov]]

    if (is_not_null(d$treat) && !allNA(d$treat) && any(na.cov)) {
      d$cov[.i][na.cov] <- d$treat[.i][na.cov]
    }
  }

  d$treat <- NULL

  d <- collap(d, dual ~ component + constraint + cov,
              FUN = sum, sort = FALSE)

  rownames(d) <- NULL

  d
}

allowable_norms <- function() {
  c("l2", "l1", "linf", "entropy", "log")
}
