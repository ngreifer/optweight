check.tols <- function(formula, data = NULL, tols, stop = FALSE) {

  if (!is.formula(formula)) stop("The argument to formula must a single formula wit hthe covariates on the right side.", call. = FALSE)
  if (missing(tols)) tols <- NA_real_
  else if (!is.numeric(tols)) stop("tols must be a numeric vector.", call. = FALSE)

  #Process treat and covs from formula and data
  tt <- delete.response(terms(formula))
  t.c <- get.covs.and.treat.from.formula(tt, data)
  formula.covs <- t.c[["reported.covs"]]
  model.covs <- t.c[["model.covs"]]

  formula.vars <- attr(attr(formula.covs, "terms"), "term.labels")
  #formula.vars <- attr(tt, "term.labels")
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

  if (any(bad.covs <- !sapply(formula.covs, is.finite)) || nrow(formula.covs) != nrow(model.covs)) {
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
  attr(out, "original.vars") <- setNames(model.vars, formula.vars[attr(model.covs, "assign")])
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
