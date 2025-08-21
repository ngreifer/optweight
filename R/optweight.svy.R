#' Estimate Targeting Weights Using Optimization
#'
#' Estimate targeting weights for covariates specified in `formula`. The
#' target means are specified with `targets` and the maximum distance
#' between each weighted covariate mean and the corresponding target mean is
#' specified by `tols`. See Zubizarreta (2015) for details of the
#' properties of the weights and the methods used to fit them.
#'
#' @inheritParams optweight.svy.fit
#' @inheritDotParams optweight.svy.fit norm min.w std.binary std.cont
#' @param formula A formula with nothing on the left hand side and the
#' covariates to be targeted on the right hand side. See [glm()] for
#' more details. Interactions and functions of covariates are allowed.
#' @param data An optional data set in the form of a data frame that contains
#' the variables in `formula`.
#' @param tols A vector of target balance tolerance values for each covariate.
#' The resulting weighted covariate means will be no further away from the
#' targets than the specified values. If only one value is supplied, it will be
#' applied to all covariates. Can also be the output of a call to [check_tols()].
#' @param s.weights A vector of sampling weights or the name of a variable in `data` that contains sampling weights. Optimization occurs on the product of the sampling weights and the estimated weights.
#' @param b.weights A vector of base weights or the name of a variable in `data` that contains base weights. If supplied, the desired norm of the distance between the estimated weights and the base weights is minimized.
#' @param verbose Whether information on the optimization problem solution
#' should be printed. This information contains how many iterations it took to
#' estimate the weights and whether the solution is optimal.
#'
#' @details
#' The optimization is performed by the lower-level function
#' [optweight.svy.fit()] using [osqp::solve_osqp()] in the
#' \pkg{osqp} package, which provides a straightforward interface to specifying
#' the constraints and objective function for quadratic optimization problems
#' and uses a fast and flexible solving algorithm.
#'
#' Weights are estimated so that the standardized differences between the
#' weighted covariate means and the corresponding targets are within the given
#' tolerance thresholds (unless `std.binary` or `std.cont` are
#' `FALSE`, in which case unstandardized mean differences are considered
#' for binary and continuous variables, respectively). For a covariate \eqn{x}
#' with specified tolerance \eqn{\delta}, the weighted mean will be within
#' \eqn{\delta} of the target. If standardized tolerance values are requested,
#' the standardization factor is the standard deviation of the covariate in the
#' whole sample. The standardization factor is always unweighted.
#'
#' See the [optweight()] help page for information on interpreting
#' dual variables and solving convergence failure.
#'
#' @returns
#' An `optweight.svy` object with the following elements:
#' \item{weights}{The estimated weights, one for each unit.}
#' \item{covs}{The
#' covariates used in the fitting. Only includes the raw covariates, which may
#' have been altered in the fitting process.}
#' \item{s.weights}{The provided
#' sampling weights.}
#' \item{call}{The function call.}
#' \item{tols}{The tolerance
#' values for each covariate.}
#' \item{duals}{A data.frame containing the dual
#' variables for each covariate. See Details for interpretation of these
#' values.}
#' \item{info}{The `info` component of the output of
#' [osqp::solve_osqp()], which contains information on the
#' performance of the optimization at termination.}
#'
#' @seealso
#' The OSQP [docs](https://osqp.org/docs/index.html) for more information on \pkg{osqp}, the underlying solver, and the options for [osqp::solve_osqp()]. [osqp::osqpSettings()] for details on options for `solve_osqp()`.
#'
#' [optweight.svy.fit()], the lower-level function that performs the fitting.
#'
#' [optweight()] for estimating weights that balance treatment groups.
#'
#' @references
#' Stellato, B., Banjac, G., Goulart, P., Boyd, S., & Bansal, V. (2024). *osqp: Quadratic Programming Solver using the 'OSQP' Library* R package version 0.6.3.3. \doi{10.32614/CRAN.package.osqp}
#'
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' cov.formula <- ~ age + educ + race + married +
#'                       nodegree
#'
#' targets <- check_targets(cov.formula, data = lalonde,
#'                         targets = c(23, 9, .3, .3, .4,
#'                                     .2, .5))
#'
#' tols <- check_tols(cov.formula, data = lalonde,
#'                    tols = 0)
#'
#' ows <- optweight.svy(cov.formula,
#'                      data = lalonde,
#'                      tols = tols,
#'                      targets = targets)
#' ows
#'
#' #Unweighted means
#' col_w_mean(ows$covs)
#'
#' #Weighted means; same as targets
#' col_w_mean(ows$covs, w = ows$weights)


#' @export
optweight.svy <- function(formula, data = NULL, tols = 0, targets = NULL, s.weights = NULL,
                          b.weights = NULL, verbose = FALSE, ...) {

  call <- match.call()

  #Process targets
  targets <- check_targets(formula, data, targets, stop = TRUE)

  #Process treat and covs from formula and data
  t.c <- get_covs_and_treat_from_formula(formula, data, sep = "_")
  reported.covs <- t.c[["reported.covs"]]
  covs <- t.c[["model.covs"]]

  if (is_not_null(t.c[["treat"]])) {
    .wrn("the variable on the left side of the formula will be ignored")
  }

  if (is_null(covs)) {
    .err("no covariates were specified")
  }

  check_missing_covs(reported.covs)
  ct <- check_tols(formula, data, tols, stop = TRUE)

  tols <- attr(ct, "internal.tols")

  #Process s.weights
  sw <- process.s.weights(s.weights, data)

  #Process b.weights
  bw <- process.b.weights(b.weights, data)

  ###Run optweight.fit
  fit_out <- optweight.svy.fit(covs = covs,
                               tols = tols,
                               targets = targets,
                               s.weights = sw,
                               b.weights = bw,
                               verbose = verbose,
                               ...)

  #Check for convergence
  status_val <- fit_out$info$status_val
  if (is_not_null(status_val) && chk::vld_number(status_val)) {
    if (status_val == -2) {
      .wrn(sprintf("the optimization failed to find a solution after %s iterations. The problem may be infeasible or more iterations may be required. Check the dual variables to see which constraints are likely causing this issue",
                   fit_out$info$iter))
    }
    else if (status_val != 1) {
      .wrn("the optimization failed to find a stable solution")
    }
  }

  test.w <- {
    if (is_null(sw)) fit_out$w
    else fit_out$w * sw
  }

  if (anyNA(test.w)) {
    .err("some weights are NA, which means something went wrong")
  }

  if (sd(test.w) / mean(test.w) > 4) {
    .wrn("some extreme weights were generated. Examine them with `summary()` and maybe relax the constraints")
  }

  #Process duals
  original.vars <- attr(ct, "original.vars")
  d <- fit_out$duals
  d$cov <- vapply(d$cov, function(c) original.vars[names(original.vars) == c][1], character(1L))
  d$dual <- with(d, ave(dual, constraint, cov, FUN = sum))
  fit_out$duals <- unique(d)

  original.vars <- attr(ct, "original.vars")
  d <- fit_out$duals
  d$cov <- vapply(d$cov, function(c) original.vars[names(original.vars) == c][1L], character(1L))
  d$dual <- with(d, ave(dual, constraint, cov, FUN = sum)) #Total effect of constraint on obj. fun. is sum of abs(duals)
  fit_out$duals <- unique(d)

  out <- list(weights = fit_out$w,
              covs = reported.covs,
              s.weights = sw,
              b.weights = bw,
              call = call,
              tols = tols,
              duals = fit_out$duals,
              info = fit_out$info)

  class(out) <- c("optweight.svy")

  out
}

#' @exportS3Method print optweight.svy
print.optweight.svy <- function(x, ...) {
  cat("An optweight.svy object\n")
  cat(sprintf(" - number of obs.: %s\n",
              length(x[["weights"]])))
  cat(sprintf(" - sampling weights: %s\n",
              if (all_the_same(x[["s.weights"]])) "none"
              else "present"))
  cat(sprintf(" - covariates: %s\n",
              if (length(names(x[["covs"]])) > 60L) "too many to name"
              else toString(names(x[["covs"]]))))

  invisible(x)
}
