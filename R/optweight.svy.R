optweight.svy <- function(formula, data = NULL, tols = 0, targets = NULL, s.weights = NULL, verbose = FALSE, ...) {

  #Process targets
  tryCatch(targets <- check.targets(formula, data, targets, stop = TRUE),
           error = function(e) {
             e. <- conditionMessage(e)
             stop(e., call. = FALSE)})

  #Process treat and covs from formula and data
  t.c <- get.covs.and.treat.from.formula(formula, data)
  reported.covs <- t.c[["reported.covs"]]
  covs <- t.c[["model.covs"]]

  if (is_not_null(t.c[["treat"]])) {
    warning("the variable on the left side of the formula will be ignored.", call. = FALSE)
  }

  if (is_null(covs)) stop("No covariates were specified.", call. = FALSE)

  if (any(bad.covs <- !sapply(reported.covs, is.finite))) {
    stop(paste0("No missing or non-finite values are allowed in the covariates. Missing or non-finite values were found in the following covariates:\n", paste(names(reported.covs)[bad.covs], collapse = ", ")), call. = FALSE)
  }

  tryCatch(ct <- check.tols(formula, data, tols, stop = TRUE),
           error = function(e) {
             e. <- conditionMessage(e)
             stop(e., call. = FALSE)})
  tols <- attr(ct, "internal.tols")

  #Process s.weights
  sw <- process.s.weights(s.weights, data)

  ###Run optweight.fit

  fit_out <- optweight.svy.fit(covs = covs,
                               tols = tols,
                               targets = targets,
                               s.weights = sw,
                               verbose = verbose,
                               ...)

  #Check for convergence
  if (fit_out$info$status_val == -2) {
    warning(paste("The optimization failed to find a solution after", fit_out$info$iter, "iterations. The problem may be infeasible or more iterations may be required. Check the dual variables to see which constraints are likely causing this issue."), call. = FALSE)
  }
  else if (fit_out$info$status_val != 1) {
    warning("The optimization failed to find a stable solution.", call. = FALSE)
  }

  warn <- FALSE
  test.w <- if (is_null(sw)) fit_out$w else fit_out$w*sw
  if (any(is.na(test.w))) stop("Some weights are NA, which means something went wrong.", call. = FALSE)
  if (warn) warning("Some extreme weights were generated. Examine them with summary() and maybe relax the constraints.", call. = FALSE)
  call <- match.call()

  #Process duals
  original.vars <- attr(ct, "original.vars")
  d <- fit_out$duals
  d$cov <- vapply(d$cov, function(c) original.vars[names(original.vars) == c][1], character(1L))
  d$dual <- with(d, ave(dual, constraint, cov, FUN = sum))
  fit_out$duals <- unique(d)

  out <- list(weights = fit_out$w,
              covs = reported.covs,
              s.weights = sw,
              call = call,
              tols = tols,
              duals = fit_out$duals,
              info = fit_out$info)

  out[vapply(out, is_null, logical(1L))] <- NULL
  class(out) <- c("optweight.svy")

  return(out)
}

print.optweight.svy <- function(x, ...) {

  cat("An optweight.svy object\n")
  cat(paste0(" - number of obs.: ", length(x[["weights"]]), "\n"))
  cat(paste0(" - sampling weights: ", ifelse(all_the_same(x[["s.weights"]]),"none", "present"), "\n"))
  cat(paste0(" - covariates: ", ifelse(length(names(x[["covs"]])) > 60, "too many to name", paste(names(x[["covs"]]), collapse = ", ")), "\n"))
  invisible(x)
}
