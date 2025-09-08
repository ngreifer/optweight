#' Estimate Targeting Weights Using Optimization
#'
#' Estimate targeting weights for covariates specified in `formula`. The
#' target means are specified with `targets` and the maximum distance
#' between each weighted covariate mean and the corresponding target mean is
#' specified by `tols`. See Zubizarreta (2015) for details of the
#' properties of the weights and the methods used to fit them.
#'
#' @inheritParams optweight
#' @inheritParams optweight.svy.fit targets
#' @inheritDotParams optweight.svy.fit norm min.w std.binary std.cont
#' @param formula a formula with nothing on the left hand side and the
#' covariates to be targeted on the right hand side. See [glm()] for
#' more details. Interactions and functions of covariates are allowed.
#' @param tols a vector of target balance tolerance values for each covariate.
#' The resulting weighted covariate means will be no further away from the
#' targets than the specified values. If only one value is supplied, it will be
#' applied to all covariates. Can also be the output of a call to [process_tols()].
#'
#' @details
#' The optimization is performed by the lower-level function
#' [optweight.svy.fit()].
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
#' \item{covs}{The covariates used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
#' \item{s.weights}{The provided sampling weights.}
#' \item{call}{The function call.}
#' \item{tols}{The tolerance values for each covariate.}
#' \item{duals}{A data.frame containing the dual variables for each covariate. See [optweight()] for interpretation of these values.}
#' \item{info}{Information about the performance of the optimization at termination.}
#'
#' @seealso
#' [optweight.svy.fit()], the lower-level function that performs the fitting.
#'
#' [optweight.fit()] for more details about the optimization options.
#'
#' [optweight()] for estimating weights that balance treatment groups.
#'
#' @references
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' cov.formula <- ~ age + educ + race + married + nodegree
#'
#' targets <- process_targets(cov.formula, data = lalonde,
#'                            targets = c(23, 9, .3, .3, .4,
#'                                        .2, .5))
#'
#' ows <- optweight.svy(cov.formula,
#'                      data = lalonde,
#'                      tols = 0,
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
                          b.weights = NULL, norm = "l2", verbose = FALSE, ...) {

  mcall <- match.call()

  formula.present <- FALSE

  if (missing(formula) && is_not_null(data)) {
    chk::chk_data(data)
    formula <- reformulate(names(data))
  }
  else if (is.data.frame(formula) && is_null(data)) {
    data <- formula
    formula <- reformulate(names(data))
  }
  else if (rlang::is_formula(formula)) {
    formula.present <- TRUE
  }
  else {
    .err("the argument to `formula` must a single formula with the covariates on the right side")
  }

  #Process treat and covs from formula and data
  t.c <- get_covs_and_treat_from_formula2(formula, data, sep = "_")
  reported.covs <- t.c[["reported.covs"]]
  covs <- t.c[["model.covs"]]

  if (is_null(covs)) {
    .err("no covariates were specified")
  }

  if (is_not_null(t.c[["treat"]])) {
    .wrn("the variable on the left side of the formula will be ignored")
  }

  check_missing_covs(reported.covs)

  chk::chk_string(norm)
  norm <- tolower(norm)
  norm <- match_arg(norm, allowable_norms())

  #Process s.weights
  sw <- process_s.weights(s.weights, data)

  #Process b.weights
  bw <- process_b.weights(b.weights, data)

  #Process tols
  tols <- .process_tols_internal(covs, tols, reported.covs,
                                 if (formula.present) "formula" else "data")

  #Process targets
  targets <- .process_targets_internal(covs, targets, sw, reported.covs,
                                       if (formula.present) "formula" else "data")

  ###Run optweight.fit
  fit_out <- optweight.svy.fit(covs = covs,
                               targets = targets,
                               tols = tols,
                               s.weights = sw,
                               b.weights = bw,
                               norm = norm,
                               verbose = verbose,
                               ...)

  test.w <- {
    if (is_null(sw)) fit_out$w
    else fit_out$w * sw
  }

  if (anyNA(test.w)) {
    .err("some weights are NA, which means something went wrong")
  }

  #Process duals
  duals <- process_duals(fit_out$duals, tols)

  out <- list(weights = fit_out$w,
              covs = reported.covs,
              s.weights = sw,
              b.weights = bw,
              norm = norm,
              call = mcall,
              tols = tols,
              duals = duals,
              info = fit_out$info)

  class(out) <- "optweight.svy"

  out
}

#' @exportS3Method print optweight.svy
print.optweight.svy <- function(x, ...) {
  print.optweight(x, ...)
}
