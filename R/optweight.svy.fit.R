#' Fitting Function for Optweight for Survey Weights
#'
#' `optweight.svy.fit()` performs the optimization for [optweight.svy()] and should, in most cases, not be
#' used directly. Little processing of inputs is performed, so they must be given
#' exactly as described below.
#'
#' @inheritParams optweight.fit
#' @param covs a numeric matrix of covariates to be targeted.
#' @param tols a vector of target balance tolerance values. Default is 0.
#' @param targets a vector of target population mean values for each covariate.
#' The resulting weights will yield sample means within `tols` units of
#' the target values for each covariate. If any target values are `NA`,
#' the corresponding variable will not be targeted and its weighted mean will
#' be wherever the weights yield the smallest variance. To ensure the weighted
#' mean for a covariate is equal to its unweighted mean (i.e., so that its
#' original mean is its target mean), its original mean must be supplied as a
#' target.
#' @param norm `character`; a string containing the name of the norm corresponding to the objective function to minimize. Allowable options include `"l1"` for the L1 norm, `"l2"` for the L2 norm (the default), `"linf"` for the L\eqn{\infty} norm, `"entropy"` for the negative entropy, and `"log"` for the sum of the negative logs. See Details at [optweight.fit()] for more information.
#' @param solver string; the name of the optimization solver to use. Allowable options depend on `norm`. Default is to use whichever eligible solver is installed, if any, or the default solver for the corresponding `norm`. See Details at [optweight.fit()] for information.
#'
#' @returns
#' An `optweight.svy.fit` object with the following elements:
#' \item{w}{The estimated weights, one for each unit.}
#' \item{duals}{A data.frame containing the dual variables for each covariate. See Zubizarreta
#' (2015) for interpretation of these values.}
#' \item{info}{A list containing information about the performance of the optimization at termination.}
#'
#' @details
#' `optweight.svy.fit()` transforms the inputs into the required inputs for the optimization functions, which are (sparse) matrices and vectors, and then supplies the outputs (the weights, dual variables, and convergence information) back to [optweight.svy()]. Little processing of inputs is performed, as this is normally handled by `optweight.svy()`.
#'
#' Target constraints are applied to the product of the estimated weights and the sampling weights. In addition, sum of the product of the estimated weights and the sampling weights is constrained to be equal to the sum of the product of the base weights and sampling weights.
#'
#' @seealso
#' [optweight.svy()] which you should use for estimating the
#' balancing weights, unless you know better.
#'
#' [optweight.fit()] for more details about the allowed norms and optimization.
#'
#' @references
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' covs <- splitfactor(lalonde[c("age", "educ", "race",
#'                               "married", "nodegree")],
#'                     drop.first = FALSE)
#'
#' targets <- c(23, 9, .3, .3, .4, .2, .5)
#'
#' ows.fit <- optweight.svy.fit(covs,
#'                              targets = targets,
#'                              norm = "l2")
#'
#' #Unweighted means
#' col_w_mean(covs)
#'
#' #Weighted means; same as targets
#' col_w_mean(covs, w = ows.fit$w)

#' @export
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
    sw <- rep.int(1, N)
  }
  else {
    chk::chk_numeric(s.weights)
    chk::chk_length(s.weights, N)

    sw <- s.weights
  }

  if (is_null(b.weights)) {
    bw <- rep.int(1, N)
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
  if (!inherits(tols, "optweight.tols") || is_null(attr(tols, "internal.tols"))) {
    tols <- .process_tols_internal(covs, tols, tols_found_in = "covs")
  }

  tols <- tols |>
    attr("internal.tols") |>
    abs()

  chk::chk_string(norm)
  norm <- tolower(norm)
  chk::chk_subset(norm, allowable_norms())

  if (norm == "linf" && !all_the_same(sw)) {
    .err('sampling weights cannot be used when `norm = "linf"`')
  }

  chk::chk_number(min.w)
  chk::chk_lt(min.w, mean(bw))

  if (norm %in% c("entropy", "log")) {
    if (any(bw <= 0)) {
      .err(sprintf("all base weights must be positive when `norm = %s`",
                   add_quotes(norm)))
    }

    min.w <- max(min.w, .Machine$double.eps)
  }

  solver <- process_solver(solver, norm, min.w)

  args <- make_process_opt_args(solver)(..., verbose = verbose)

  constraint_df <- expand.grid(time = 1L,
                               type = c("range_w", "mean_w", "target"),
                               constraint = list(NULL),
                               stringsAsFactors = FALSE,
                               KEEP.OUT.ATTRS = FALSE)

  range_cons <- constraint_range_w(sw, min.w)

  bin.covs <- is_binary_col(covs)

  n <- sum(sw * bw)

  sds <- sqrt(col.w.v(covs, w = sw, bin.vars = bin.covs))

  targeted <- !is.na(targets)

  vars.to.standardize <- rep_with(FALSE, tols)
  if (std.binary) vars.to.standardize[bin.covs] <- TRUE
  if (std.cont) vars.to.standardize[!bin.covs] <- TRUE

  to_std <- which(vars.to.standardize & targeted & !is.na(sds) & !check_if_zero(sds))

  if (is_not_null(to_std)) {
    covs[, to_std] <- mat_div(covs[, to_std, drop = FALSE], sds[to_std])
    targets[to_std] <- targets[to_std] / sds[to_std]
  }

  constraint_df[["constraint"]] <- list(
    range_w = range_cons,
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
              duals = duals[[1L]],
              info = opt_out$info_out,
              out = opt_out$out)

  class(out) <- "optweight.svy.fit"

  out
}
