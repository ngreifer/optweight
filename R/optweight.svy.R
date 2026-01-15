#' Stable Balancing Weights for Generalization
#'
#' Estimates stable balancing weights to generalize a sample characterized by supplied covariates to a given target population. The target means are specified with `targets` and the maximum distance between each weighted covariate mean. See Jackson et al. (2021) for details of the properties of the weights and the methods used to fit them.
#'
#' @inheritParams optweight
#' @param formula a formula with nothing on the left hand side and the covariates to be targeted on the right hand side. Interactions and functions of covariates are allowed. Can be omitted, in which case all variables in `data` are assumed targeted. If `data` is `NULL` and `formula` is a data.frame, `data` will be replaced with `formula`.
#' @param tols a vector of target balance tolerance values for each covariate. The resulting weighted covariate means will be no further away from the targets than the specified values. If only one value is supplied, it will be applied to all covariates. Can also be the output of a call to [process_tols()]. Default is 0 for all covariates.
#' @param targets a vector of target population mean values for each covariate. The resulting weights will yield sample means within `tols` units of the target values for each covariate. If any target values are `NA`, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance. To ensure the weighted mean for a covariate is equal to its unweighted mean (i.e., so that its original mean is its target mean), its original mean must be supplied as a target. For factor variables, a target value must be specified for each level of the factor, and these values must add up to 1. Can also be the output of a call to [process_targets()].
#' @param covs a numeric matrix of covariates to be targeted.
#'
#' @returns
#' For `optweight.svy()`, an `optweight.svy` object with the following elements:
#' \item{weights}{The estimated weights, one for each unit.}
#' \item{covs}{The covariates used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
#' \item{s.weights}{The provided sampling weights.}
#' \item{call}{The function call.}
#' \item{tols}{The tolerance values for each covariate.}
#' \item{duals}{A data.frame containing the dual variables for each covariate. See [optweight()] for interpretation of these values.}
#' \item{info}{A list containing information about the performance of the optimization at termination.}
#' \item{norm}{The `norm` used.}
#' \item{solver}{The `solver` used.}
#'
#' For `optweight.svy.fit()`, an `optweight.svy.fit` object with the following elements:
#' \item{w}{The estimated weights, one for each unit.}
#' \item{duals}{A data.frame containing the dual variables for each covariate.}
#' \item{info}{A list containing information about the performance of the optimization at termination.}
#' \item{norm}{The `norm` used.}
#' \item{solver}{The `solver` used.}
#'
#' @details
#' `optweight.svy()` is the primary user-facing function for estimating stable balancing weights for generalization to a target population. The optimization is performed by the lower-level function `optweight.svy.fit()`, which transforms the inputs into the required inputs for the optimization functions and then supplies the outputs (the weights, dual variables, and convergence information) back to `optweight.svy()`. Little processing of inputs is performed by `optweight.svy.fit()`, as this is normally handled by `optweight.svy()`.
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
#' Target constraints are applied to the product of the estimated weights and the sampling weights. In addition, sum of the product of the estimated weights and the sampling weights is constrained to be equal to the sum of the product of the base weights and sampling weights.
#'
#' See [optweight()] for information on `norm`, `solver`, and convergence failure.
#'
#' @seealso
#' [optweight()] for estimating weights that balance treatment groups.
#'
#' [process_targets()] for specifying the covariate target means supplied to `targets`.
#'
#' @references
#' Jackson, D., Rhodes, K., & Ouwens, M. (2021). Alternative weighting schemes when performing matching-adjusted indirect comparisons. *Research Synthesis Methods*, 12(3), 333–346. \doi{10.1002/jrsm.1466}
#'
#' @examplesIf rlang::is_installed("cobalt")
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' cov.names <- c("age", "educ", "race",
#'                "married", "nodegree")
#'
#' targets <- c(age = 23,
#'              educ = 9,
#'              race_black = .3,
#'              race_hispan = .3,
#'              race_white = .4,
#'              married = .2,
#'              nodegree = .5)
#'
#' ows <- optweight.svy(lalonde[cov.names],
#'                      targets = targets)
#' ows
#'
#' # Unweighted means
#' col_w_mean(lalonde[cov.names])
#'
#' # Weighted means; same as targets
#' col_w_mean(lalonde[cov.names],
#'            w = ows$weights)

#' @export
optweight.svy <- function(formula, data = NULL, tols = 0, targets = NULL, s.weights = NULL,
                          b.weights = NULL, norm = "l2", min.w = 1e-8, verbose = FALSE, ...) {

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
    .err("the argument to {.arg formula} must a single formula with the covariates on the right side")
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

  #Process s.weights
  sw <- process_s.weights(s.weights, data)

  #Process b.weights
  bw <- process_b.weights(b.weights, data)

  #Process tols
  tols <- .process_tols_internal(covs, tols, reported.covs,
                                 if (formula.present) "formula" else "data",
                                 tols_arg = "tols")

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
                               min.w = min.w,
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
              norm = fit_out$norm,
              call = mcall,
              tols = tols,
              duals = duals,
              info = fit_out$info,
              solver = fit_out$solver)

  class(out) <- "optweight.svy"

  out
}

#' @export
#' @rdname optweight.svy
optweight.svy.fit <- function(covs, targets, tols = 0, s.weights = NULL, b.weights = NULL,
                              norm = "l2", std.binary = FALSE, std.cont = TRUE,
                              min.w = 1e-8, verbose = FALSE, solver = NULL, ...) {

  chk::chk_not_missing(covs, "`covs`")
  chk::chk_not_missing(targets, "`targets`")

  if (!is.numeric(covs) && (!is.data.frame(covs) || !all(apply(covs, 2L, is.numeric)))) {
    .err("all covariates must be numeric")
  }

  covs <- as.matrix(covs)

  N <- nrow(covs)

  if (is_null(s.weights)) {
    sw <- alloc(1, N)
  }
  else {
    chk::chk_numeric(s.weights)
    chk::chk_length(s.weights, N)

    sw <- s.weights
  }

  if (is_null(b.weights)) {
    bw <- alloc(1, N)
  }
  else {
    chk::chk_numeric(b.weights)
    chk::chk_length(b.weights, N)

    bw <- b.weights
  }

  #Process targets
  if (!inherits(tols, "optweight.targets")) {
    targets <- .process_targets_internal(covs, targets = targets, sw = sw,
                                         targets_found_in = "covs")
  }

  #Process tols
  if (!inherits(tols, "optweight.tols") || is_null(.attr(tols, "internal.tols"))) {
    tols <- .process_tols_internal(covs, tols, tols_found_in = "covs",
                                   tols_arg = "tols")
  }

  tols <- tols |>
    .attr("internal.tols") |>
    abs()

  #Process norm
  norm <- process_norm(norm, sw, bw)

  #Process min.w
  min.w <- process_min.w(min.w, norm, bw)

  #Process solver
  solver <- process_solver(solver, norm, min.w)

  #Process args
  args <- make_process_opt_args(solver)(..., verbose = verbose)

  range_cons <- constraint_range_w(sw, min.w)

  constraint_df <- expand.grid(time = 0L,
                               type = "range_w",
                               constraint = list(range_cons),
                               stringsAsFactors = FALSE,
                               KEEP.OUT.ATTRS = FALSE) |>
    rbind(expand.grid(time = 1L,
                      type = c("mean_w", "target"),
                      constraint = list(NULL),
                      stringsAsFactors = FALSE,
                      KEEP.OUT.ATTRS = FALSE))

  bin.covs <- is_binary_col(covs)

  n <- sum(sw * bw)

  sds <- col.w.sd(covs, w = sw, bin.vars = bin.covs)

  targeted <- !is.na(targets)

  vars.to.standardize <- rep_with(FALSE, tols)
  if (std.binary) vars.to.standardize[bin.covs] <- TRUE
  if (std.cont) vars.to.standardize[!bin.covs] <- TRUE

  to_std <- which(vars.to.standardize & targeted & !is.na(sds) & !check_if_zero(sds))

  if (is_not_null(to_std)) {
    covs[, to_std] <- ss(covs, j = to_std) %r/% sds[to_std]
    targets[to_std] <- targets[to_std] / sds[to_std]
  }

  constraint_df[["constraint"]][constraint_df[["time"]] == 1L] <- list(
    mean_w = constraint_mean_w_svy(sw, n),
    target = constraint_target_svy(covs, sw, targets, targeted, tols, n)
  )

  constraint_df <- constraint_df |>
    prep_constraint_df(norm, bw, sw) |>
    prep_constraint_df_for_solver(solver)

  objective <- prep_objective(norm, bw, sw)

  opt_out <- opt_fit(constraint_df, objective, args, N,
                     solver = solver)

  w <- extract_weights(opt_out, N, min.w, range_cons)

  duals <- extract_duals(constraint_df, opt_out$dual_out)

  out <- list(w = w,
              duals = duals,
              info = opt_out$info_out,
              out = opt_out$out,
              norm = norm,
              solver = solver)

  class(out) <- "optweight.svy.fit"

  out
}

#' @exportS3Method print optweight.svy
print.optweight.svy <- print.optweight
