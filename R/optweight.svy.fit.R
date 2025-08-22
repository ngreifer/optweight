#' Fitting Function for Optweight for Survey Weights
#'
#' `optweight.svy.fit()` performs the optimization (via \CRANpkg{osqp}) for [optweight.svy()] and should, in most cases, not be
#' used directly. Little processing of inputs is performed, so they must be given
#' exactly as described below.
#'
#' @inheritParams optweight.fit
#' @param covs A matrix of covariates to be targeted. Should must be numeric
#' but does not have to be full rank.
#' @param tols A vector of target balance tolerance values.
#' @param targets A vector of target population mean values for each covariate.
#' The resulting weights will yield sample means within `tols` units of
#' the target values for each covariate. If any target values are `NA`,
#' the corresponding variable will not be targeted and its weighted mean will
#' be wherever the weights yield the smallest variance. To ensure the weighted
#' mean for a covariate is equal to its unweighted mean (i.e., so that its
#' original mean is its target mean), its original mean must be supplied as a
#' target.
#'
#' @returns
#' An `optweight.svy.fit` object with the following elements:
#' \item{w}{The estimated weights, one for each unit.}
#' \item{duals}{A data.frame containing the dual variables for each covariate. See Zubizarreta
#' (2015) for interpretation of these values.}
#' \item{info}{The `info`
#' component of the output of [osqp::solve_osqp()], which contains
#' information on the performance of the optimization at termination.}
#'
#' @details
#' `optweight.svy.fit()` transforms the inputs into the required inputs for
#' [osqp::solve_osqp()], which are (sparse) matrices and vectors, and
#' then supplies the outputs (the weights, duals variables, and convergence
#' information) back to [optweight.svy()]. Little processing of inputs is
#' performed, as this is normally handled by `optweight.svy()`.
#'
#' The default values for some of the parameters sent to [osqp::solve_osqp()] are not the same as those in [osqp::osqpSettings()]. The following are the differences: `max_iter` is set to 20000, `eps_abs` and `eps_rel` are set to 1e-8 (i.e., \eqn{10^{-8}}), and `adaptive_rho_interval` is set to 10. All other values are the same.
#'
#' @seealso
#' [optweight.svy()] which you should use for estimating the
#' balancing weights, unless you know better.
#'
#' The OSQP [docs](https://osqp.org/docs/index.html) for more information on \pkg{osqp}, the underlying solver, and the options for [osqp::solve_osqp()]. [osqp::osqpSettings()] for details on options for `solve_osqp()`.
#'
#' @references
#' Wang, Y., & Zubizarreta, J. R. (2020). Minimal dispersion approximately balancing weights: Asymptotic properties and practical considerations. *Biometrika*, 107(1), 93–105. \doi{10.1093/biomet/asz050}
#'
#' Yiu, S., & Su, L. (2018). Covariate association eliminating weights: a unified weighting framework for causal effect estimation. *Biometrika*. \doi{10.1093/biomet/asy015}
#'
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' covs <- splitfactor(lalonde[c("age", "educ", "race",
#'                   "married", "nodegree")],
#'                   drop.first = FALSE)
#'
#' targets <- c(23, 9, .3, .3, .4, .2, .5)
#'
#' tols <- rep(0, 7)
#'
#' ows.fit <- optweight.svy.fit(covs,
#'                              tols = tols,
#'                              targets = targets,
#'                              norm = "l2")
#'
#' #Unweighted means
#' col_w_mean(covs)
#'
#' #Weighted means; same as targets
#' col_w_mean(covs, w = ows.fit$w)

#' @export
optweight.svy.fit <- function(covs, tols = 0, targets, s.weights = NULL, b.weights = NULL,
                              norm = "l2", std.binary = FALSE, std.cont = TRUE,
                              min.w = 1e-8, verbose = FALSE, ...) {

  chk::chk_not_missing(covs, "`covs`")
  chk::chk_not_missing(targets, "`targets`")

  if (!is.numeric(covs) && (!is.data.frame(covs) || !all(apply(covs, 2L, is.numeric)))) {
    .err("all covariates must be numeric")
  }

  covs <- as.matrix(covs)

  chk::chk_numeric(targets)
  chk::chk_numeric(tols)

  chk::chk_string(norm)
  norm <- tolower(norm)
  chk::chk_subset(norm, c("l2", "l1", "linf"))

  chk::chk_number(min.w)
  chk::chk_lt(min.w, 1)

  if (length(tols) == 1L) {
    tols <- rep.int(tols, ncol(covs))
  }

  if (length(targets) != ncol(covs)) {
    .err("`targets` must have the same number of values as there are baseline covariates")
  }

  N <- nrow(covs)

  sw <- {
    if (is_null(s.weights)) rep.int(1, N)
    else as.numeric(s.weights)
  }

  bw <- {
    if (is_null(b.weights)) rep.int(1, N)
    else as.numeric(b.weights)
  }

  if (norm == "linf" && !all_the_same(sw)) {
    .err("the L-inf norm cannot be used with sampling weights")
  }

  args <- ...mget(names(formals(osqp::osqpSettings)))

  for (e in c("eps_abs", "eps_rel")) {
    if (!utils::hasName(args, e)) {
      args[[e]] <- ...get("eps", 1e-5)
    }
  }

  if (!utils::hasName(args, "max_iter")) {
    args[["max_iter"]] <- 2e5
  }

  if (!utils::hasName(args, "adaptive_rho_interval")) {
    args[["adaptive_rho_interval"]] <- 10L
  }

  if (!utils::hasName(args, "polish")) {
    args[["polish"]] <- TRUE
  }

  args[["verbose"]] <- verbose

  constraint_df <- expand.grid(time = 0,
                               type = c("range_w", "mean_w", "target"),
                               constraint = list(NULL),
                               stringsAsFactors = FALSE,
                               KEEP.OUT.ATTRS = FALSE)

  bin.covs <- is_binary_col(covs)

  means <- col.w.m(covs, w = sw)
  sds <- sqrt(col.w.v(covs, w = sw, bin.vars = bin.covs))

  n <- sum(sw) #N

  targeted <- !is.na(targets)

  tols <- abs(tols)

  vars.to.standardize <- rep_with(FALSE, tols)
  if (std.binary) vars.to.standardize[bin.covs] <- TRUE
  if (std.cont) vars.to.standardize[!bin.covs] <- TRUE

  to_std <- which(vars.to.standardize & targeted & !is.na(sds) & !check_if_zero(sds))

  if (is_not_null(to_std)) {
    covs[, to_std] <- mat_div(covs[, to_std, drop = FALSE], sds[to_std])
    targets[to_std] <- targets[to_std] / sds[to_std]
  }

  constraint_df[["constraint"]] <- list(
    range_w = constraint_range_w(sw, min.w),
    mean_w = constraint_mean_w_svy(sw, n),
    target = constraint_target_svy(covs, sw, targets, targeted, tols)
  )

  if (norm == "l2") {
    objective <- objective_L2(bw, sw)
  }
  else if (norm == "l1") {
    objective <- objective_L1(bw, sw)

    constraint_df <- modify_constraints_L1(constraint_df, bw) |>
      rbind(expand.grid(time = 0, type = "conversion",
                        constraint = list(constraint_conversion_L1(bw, sw)),
                        stringsAsFactors = FALSE,
                        KEEP.OUT.ATTRS = FALSE))
  }
  else if (norm == "linf") {
    objective <- objective_Linf(bw, sw)

    constraint_df <- modify_constraints_Linf(constraint_df, bw) |>
      rbind(expand.grid(time = 0, type = "conversion",
                        constraint = list(constraint_conversion_Linf(bw, sw)),
                        stringsAsFactors = FALSE,
                        KEEP.OUT.ATTRS = FALSE))
  }

  constraint_df$nc <- lengths(grab(constraint_df[["constraint"]], "L"))
  constraint_df$nc_cum <- cumsum(constraint_df$nc)

  out <- osqp::solve_osqp(P = objective$P,
                          q = objective$q,
                          A = combine_constraints("A", constraint_df[["constraint"]]),
                          l = combine_constraints("L", constraint_df[["constraint"]]),
                          u = combine_constraints("U", constraint_df[["constraint"]]),
                          pars = do.call(osqp::osqpSettings, args))

  #Get dual vars for balance and target constraints
  target_indices <- unlist(lapply(which(constraint_df$type == "target"), function(i) {
    seq(constraint_df$nc_cum[i] - constraint_df$nc[i], constraint_df$nc_cum[i])[-1L]
  }))

  w <- out$x[seq_len(N)]

  if (abs(min.w) < .Machine$double.eps) {
    w[abs(w) < .Machine$double.eps] <- 0
  }

  w[w < min.w] <- min.w

  #Duals
  duals <- NULL

  ti <- which(constraint_df[["type"]] == "target")[1L]
  if (constraint_df[["nc"]][ti] > 0) {
    cons <- constraint_df[["constraint"]][[ti]][-(1:3)]
    target_ind <- seq(constraint_df$nc_cum[ti] - constraint_df$nc[ti],
                      constraint_df$nc_cum[ti])[-1L]

    duals <- data.frame(constraint = "target",
                        cov = if_null_then(cons$covs, NA_character_),
                        dual = abs(out$y[target_ind]))
  }

  opt_out <- list(w = w,
                  duals = duals,
                  info = out$info,
                  out = out)
  class(opt_out) <- "optweight.svy.fit"

  opt_out
}
