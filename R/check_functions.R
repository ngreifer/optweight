check.tols <- function(formula, data = NULL, tols, stop = FALSE) {

  if (!is.formula(formula)) stop("The argument to formula must a single formula with the covariates on the right side.", call. = FALSE)
  if (missing(tols)) tols <- NA_real_
  else if (!is.numeric(tols)) stop("tols must be a numeric vector.", call. = FALSE)

  #Process treat and covs from formula and data
  tt <- delete.response(terms(formula))
  t.c <- get.covs.and.treat.from.formula(tt, data, terms = TRUE, sep = "_")
  formula.covs <- t.c[["reported.covs"]]
  model.covs <- t.c[["model.covs"]]

  formula.vars <- attr(attr(formula.covs, "terms"), "term.labels")
  if (is_null(formula.vars)) {
    formula.vars <- "(Intercept)"
    attr(model.covs, "assign") <- 1
  }
  model.vars <- colnames(model.covs)

  if (length(tols) == 1) {
    onetol <- TRUE
    tols <- rep(tols, length(formula.vars))
  }
  else {
    onetol <- FALSE
  }

  if (is_null(model.covs)) stop("No covariates were specified.", call. = FALSE)

  if (any(bad.covs <- sapply(formula.covs, function(x) any(!is.finite(x)))) || nrow(formula.covs) != nrow(model.covs)) {
    stop(paste0("No missing or non-finite values are allowed in the covariates. Missing or non-finite values were found in the following covariates:\n", paste(names(formula.covs)[bad.covs], collapse = ", ")), call. = FALSE)
  }

  if (length(tols) != length(formula.vars)) {
    if (stop) {
      stop(paste0("tols must contain ", length(formula.vars), " numbers. ", length(tols),
                  " were included."), call. = FALSE)
    }
    message(paste0("tols must contain ", length(formula.vars), " numbers. ", length(tols),
                   " were included.\nAll values in the output have been assigned NA."))
    user.tols <- setNames(rep(NA_real_, length(formula.vars)), formula.vars)
    internal.tols <- setNames(rep(NA_real_, length(model.vars)), model.vars)
  }
  else {
    user.tols <- setNames(tols, formula.vars)
    internal.tols <- setNames(tols[attr(model.covs, "assign")], model.vars)
    if (!stop) {
      if (any(attr(terms(formula), "order") > 1)) {
        #message("tols look okay, but interactions were present in the formula, so make sure the order is correct.")
      }
      else {
        #message("tols look okay, but you should print them to make sure they're correct.")
      }
    }
  }

  out <- user.tols
  attr(out, "internal.tols") <- internal.tols
  attr(out, "original.vars") <- setNames(formula.vars[attr(model.covs, "assign")], model.vars)
  class(out) <- "optweight.tols"
  return(out)
}

print.optweight.tols <- function(x, internal = FALSE, digits = 5, ...) {
  tols <- x
  internal.tols <- attr(tols, "internal.tols")
  attributes(tols) <- NULL
  names(tols) <- names(x)
  if (all(is.na(tols))) {
    cat("- vars:\n\t")
    cat(paste(names(tols), collapse = "   "))
  }
  else {
    cat("- tols:\n")
    print(round(tols, digits))
    if (internal) {
      cat("\n- tols used internally by optweight:\n")
      print(round(internal.tols, digits))
    }
  }
  invisible(x)
}

check.targets <- function(formula, data = NULL, targets, stop = FALSE) {
  if (!is.formula(formula)) stop("The argument to formula must a single formula with the covariates on the right side.", call. = FALSE)

  #Process treat and covs from formula and data
  tt <- delete.response(terms(formula))
  t.c <- get.covs.and.treat.from.formula(tt, data, terms = TRUE, sep = "_")
  formula.covs <- t.c[["reported.covs"]]
  model.covs <- t.c[["model.covs"]]

  formula.vars <- attr(attr(formula.covs, "terms"), "term.labels")
  if (is_null(formula.vars)) {
    formula.vars <- "(Intercept)"
    attr(model.covs, "assign") <- 1
  }
  model.vars <- colnames(model.covs)

  if (is_null(model.covs)) stop("No covariates were specified.", call. = FALSE)

  if (any(bad.covs <- sapply(formula.covs, function(x) any(!is.finite(x)))) || nrow(formula.covs) != nrow(model.covs)) {
    stop(paste0("No missing or non-finite values are allowed in the covariates. Missing or non-finite values were found in the following covariates:\n", paste(names(formula.covs)[bad.covs], collapse = ", ")), call. = FALSE)
  }

  missing.targets <- missing(targets)
  if (missing.targets) targets <- colMeans(model.covs)
  else if (is_null(targets)) targets <- rep(NA_real_, length(model.vars))
  else if (!is.numeric(targets)) stop("targets must be a numeric vector.", call. = FALSE)

  if (length(targets) != length(model.vars)) {
    if (stop) {
      stop(paste0("targets must contain ", length(model.vars), " numbers. ", length(targets),
                  " were included."), call. = FALSE)
    }
    message(paste0("targets must contain ", length(model.vars), " numbers. ", length(targets),
                   " were included.\nAll values in the output have been assigned NA."))
    internal.targets <- setNames(rep(NA_real_, length(model.vars)), model.vars)
  }
  else {
    internal.targets <- setNames(targets, model.vars)
    original.variables <- setNames(formula.vars[attr(model.covs, "assign")], model.vars)
    for (v in formula.vars) {
      if (attr(terms(formula.covs), "order")[formula.vars == v] == 1 &&
          attr(terms(formula.covs), "dataClasses")[v] == "factor") {
        if (!any(is.na(internal.targets[original.variables == v])) &&
            !check_if_zero(sum(internal.targets[original.variables == v]) - 1)) {
          stop(paste("The target values for", v, "must add up to 1."), call. = FALSE)
        }
      }
    }
    if (!all(is.na(targets)) && any(sapply(formula.vars, function(v) {
      if (attr(terms(formula.covs), "order")[formula.vars == v] > 1) {
        vars.in.interaction <- rownames(attr(terms(formula.covs), "factors"))[attr(terms(formula.covs), "factors")[, v] == 1]
        if (sum(attr(terms(formula.covs), "dataClasses")[vars.in.interaction] == "factor") > 1) {
          return(TRUE)
        }
        else return(FALSE)
      }
      else return(FALSE)
    }))) warning("Interactions between factor variables were entered, but check.targets cannot verify whether the target values are suitable. See ?check.targets for details.", call. = FALSE)

    if (!stop) {
      if (any(attr(terms(formula), "order") > 1)) {
        #message("targets look okay, but interactions were present in the formula, so make sure the order is correct.")
      }
      else {
        #message("targets look okay, but you should print them to make sure they're correct.")
      }
    }
  }

  out <- internal.targets
  attr(out, "original.vars") <- setNames(formula.vars[attr(model.covs, "assign")], model.vars)
  attr(out, "ATE") <- (missing.targets || isTRUE(all(check_if_zero(colMeans(model.covs) - internal.targets))))
  class(out) <- "optweight.targets"
  return(out)
}

print.optweight.targets <- function(x, digits = 5, ...) {
  targets <- x
  attributes(targets) <- NULL
  names(targets) <- names(x)
  if (all(is.na(targets))) {
    cat("- vars:\n\t")
    cat(paste(names(targets), collapse = "   "))
  }
  else {
    cat("- targets:\n")
    print(round(targets, digits))
  }
  invisible(x)
}
