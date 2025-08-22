#' Fitting Function for Stable Balancing Weights
#'
#' `optweight.fit()` performs the optimization (via \CRANpkg{osqp}) for [optweight()] and should, in most cases, not be used directly. Little processing of inputs is performed, so they must be given exactly as described below.
#'
#' @param treat.list A list containing one vector of treatment statuses for each time point. Non-numeric (i.e., factor or character) vectors are allowed.
#' @param covs.list A list containing one matrix of covariates to be balanced for each time point. All matrices must be numeric but do not have to be full rank.
#' @param tols A list containing one vector of balance tolerance values for each time point.
#' @param estimand The desired estimand, which determines the target population. For binary treatments, can be "ATE", "ATT", "ATC", or `NULL`. For multi-category treatments, can be "ATE", "ATT", or `NULL`. For continuous treatments, can be "ATE" or `NULL`. The default for both is "ATE". For longitudinal treatments, only "ATE" is supported. `estimand` is ignored when `targets` is non-`NULL`. If both `estimand` and `targets` are `NULL`, no targeting will take place. See Details.
#' @param targets A vector of target population mean values for each baseline covariate. The resulting weights will yield sample means within `tols`/2 units of the target values for each covariate. If `NULL` or all `NA`, `estimand` will be used to determine targets. Otherwise, `estimand` is ignored. If any target values are `NA`, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance.
#' @param s.weights A vector of sampling weights. Optimization occurs on the product of the sampling weights and the estimated weights.
#' @param b.weights A vector of base weights. Default is a vector of 1s. The desired norm of the distance between the estimated weights and the base weights is minimized.
#' @param focal When multi-categorical treatments are used and the "ATT" is requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to resemble the focal group.
#' @param norm A string containing the name of the norm corresponding to the objective function to minimize. The options are `"l1"` for the L1 norm, `"l2"` for the L2 norm (the default), and `"linf"` for the L\eqn{\infty} norm. The L1 norm minimizes the average absolute distance between each weight and the base weights; the L2 norm minimizes the average squared distance between each weight and the base weights; the L\eqn{\infty} norm minimizes the largest absolute distance between each weight and the base weights. The L2 norm has a direct correspondence with the effective sample size, making it ideal if this is your criterion of interest.
#' @param std.binary,std.cont `logical`; whether the tolerances are in standardized mean units (`TRUE`) or raw units (`FALSE`) for binary variables and continuous variables, respectively. The default is `FALSE` for `std.binary` because raw proportion differences make more sense than standardized mean difference for binary variables. These arguments are analogous to the `binary` and `continuous` arguments in `bal.tab()` in \pkg{cobalt}.
#' @param min.w A single `numeric` value less than 1 for the smallest allowable weight. Some analyses require nonzero weights for all units, so a small, nonzero minimum may be desirable. Doing so will likely (slightly) increase the variance of the resulting weights depending on the magnitude of the minimum. The default is 1e-8, which does not materially change the properties of the weights from a minimum of 0 but prevents warnings in some packages that use weights to estimate treatment effects.
#' @param verbose Whether information on the optimization problem solution should be printed. This information contains how many iterations it took to estimate the weights and whether the solution is optimal.
#' @param force Stable balancing weights are currently not valid for use with longitudinal treatments, and will produce an error message if attempted. Set to `TRUE` to bypass this error message.
#' @param \dots Options that are passed to [osqp::osqpSettings()] for use in the `par` arguments of [osqp::solve_osqp()]. See Details for defaults.
#'
#' @returns
#' An `optweight.fit` object with the following elements:
#'   \item{w}{The estimated weights, one for each unit.}
#'   \item{duals}{A data.frame containing the dual variables for each covariate, or a list thereof. See Zubizarreta (2015) for interpretation of these values.}
#'   \item{info}{The `info` component of the output of [osqp::solve_osqp()], which contains information on the performance of the optimization at termination.}
#'
#' @details
#' `optweight.fit()` transforms the inputs into the required inputs for [osqp::solve_osqp()], which are (sparse) matrices and vectors, and then supplies the outputs (the weights, dual variables, and convergence information) back to [optweight()]. Little processing of inputs is performed, as this is normally handled by `optweight()`.
#'
#'   The default values for some of the parameters sent to [osqp::solve_osqp()] are not the same as those in [osqp::osqpSettings()]. The following are the differences: `max_iter` is set to 20000, `eps_abs` and `eps_rel` are set to 1e-8 (i.e., \eqn{10^{-8}}), and `adaptive_rho_interval` is set to 10. All other values are the same.
#'
#'   Note that stable balancing weights with longitudinal treatments are not valid and should not be used until further research is done.
#'
#' @references
#' Wang, Y., & Zubizarreta, J. R. (2020). Minimal dispersion approximately balancing weights: Asymptotic properties and practical considerations. *Biometrika*, 107(1), 93–105. \doi{10.1093/biomet/asz050}
#'
#' Yiu, S., & Su, L. (2018). Covariate association eliminating weights: a unified weighting framework for causal effect estimation. *Biometrika*. \doi{10.1093/biomet/asy015}
#'
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' @seealso
#' [optweight()] which you should use for estimating the balancing weights, unless you know better.
#'
#' The OSQP [docs](https://osqp.org/docs/index.html) for more information on \pkg{osqp}, the underlying solver, and the options for [osqp::solve_osqp()]. [osqp::osqpSettings()] for details on options for `solve_osqp()`.
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' treat.list <- list(lalonde$treat)
#' covs.list <- list(splitfactor(lalonde[2:8], drop.first = "if2"))
#' tols.list <- list(rep(.01, ncol(covs.list[[1]])))
#'
#' ow.fit <- optweight.fit(treat.list,
#'                         covs.list,
#'                         tols = tols.list,
#'                         estimand = "ATE",
#'                         norm = "l2")
#'
#' @export
optweight.fit <- function(treat.list, covs.list, tols, estimand = "ATE", targets = NULL,
                          s.weights = NULL, b.weights = NULL, focal = NULL, norm = "l2",
                          std.binary = FALSE, std.cont = TRUE, min.w = 1e-8, verbose = FALSE,
                          force = FALSE, ...) {

  chk::chk_not_missing(treat.list, "`treat.list`")
  chk::chk_not_missing(covs.list, "`covs.list`")
  chk::chk_not_missing(tols, "`tols`")

  chk::chk_list(treat.list)
  chk::chk_list(covs.list)
  chk::chk_list(tols)

  if (length(covs.list) > 1L) {
    chk::chk_flag(force)
    if (!force) {
      .err("optweights are currently not valid for longitudinal treatments. Set `force = TRUE` to bypass this message at your own risk")
    }
  }

  if (!all_apply(covs.list, function(c) all(apply(c, 2L, is.numeric)))) {
    .err("all covariates must be numeric")
  }

  chk::chk_string(norm)
  norm <- tolower(norm)
  chk::chk_subset(norm, c("l2", "l1", "linf"))

  times <- seq_along(covs.list)

  if (is_not_null(estimand)) {
    chk::chk_string(estimand)
    estimand <- toupper(estimand)
    chk::chk_subset(estimand, c("ATE", "ATT", "ATC"))

    if (length(times) > 1L && !identical(estimand, "ATE")) {
      .err("only the ATE or specified targets are compatible with longitudinal treatments")
    }
  }

  chk::chk_number(min.w)
  chk::chk_lt(min.w, 1)

  if (length(tols) == 1L) {
    tols <- rep_with(tols, times)
  }

  for (i in which(lengths(tols) == 1L)) {
    tols[[i]] <- rep.int(tols[[i]], ncol(covs.list[[i]]))
  }

  if (is_not_null(estimand) || is_null(targets) || all(is.na(targets))) {
    targets <- NULL
  }
  else if (!is.atomic(targets) || !is.numeric(targets)) {
    .err("`targets` must be a vector of target values for each baseline covariate")
  }
  else if (length(targets) != ncol(covs.list[[1L]])) {
    .err("`targets` must have the same number of values as there are baseline covariates")
  }

  N <- length(treat.list[[1L]])

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

  corr.type <- "pearson"

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

  treat.types <- character(length(times))

  unique.treats <- n <- bin.covs.list <- means <- sds <- targets.list <- targeted <-
    balanced <- treat.sds <- treat.means <- tols.list <- make_list(length(times))

  constraint_df <- expand.grid(time = times,
                               type = c("range_w", "mean_w", "balance", "target"),
                               constraint = list(NULL),
                               stringsAsFactors = FALSE,
                               KEEP.OUT.ATTRS = FALSE)

  for (i in times) {
    covs.list[[i]] <- as.matrix(covs.list[[i]])

    bin.covs.list[[i]] <- is_binary_col(covs.list[[i]])

    means[[i]] <- col.w.m(covs.list[[i]], w = sw)

    treat.types[i] <- {
      if (chk::vld_character_or_factor(treat.list[[i]]) || is_binary(treat.list[[i]])) "cat"
      else "cont"
    }

    if (treat.types[i] == "cat") {
      treat.list[[i]] <- as.character(treat.list[[i]])

      unique.treats[[i]] <- sort(unique(treat.list[[i]]))

      n[[i]] <- vapply(unique.treats[[i]],
                       function(t) sum(sw[treat.list[[i]] == t]),
                       numeric(1L))

      in_focal <- {
        if (is_not_null(estimand) && estimand %in% c("ATT", "ATC"))
          in_focal <- which(treat.list[[i]] == focal)
        else
          NULL
      }

      sds[[i]] <- {
        if (is_not_null(estimand) && estimand %in% c("ATT", "ATC"))
          sqrt(col.w.v(covs.list[[i]][in_focal, , drop = FALSE],
                       w = sw[in_focal],
                       bin.vars = bin.covs.list[[i]]))
        else
          sqrt(colMeans(do.call("rbind", lapply(unique.treats[[i]], function(t) {
            in_treat <- which(treat.list[[i]] == t)

            col.w.v(covs.list[[i]][in_treat, , drop = FALSE],
                    w = sw[in_treat], bin.vars = bin.covs.list[[i]])
          }))))
      }
    }
    else {
      treat.list[[i]] <- as.numeric(treat.list[[i]])

      n[[i]] <- sum(sw) #N

      sds[[i]] <- sqrt(col.w.v(covs.list[[i]], w = sw,
                               bin.vars = bin.covs.list[[i]]))
    }

    targets.list[[i]] <- {
      if (i > 1L || (is_null(targets) && is_null(estimand)))
        rep.int(NA_real_, ncol(covs.list[[i]]))
      else if (is_not_null(targets))
        targets
      else if (estimand %in% c("ATT", "ATC"))
        col.w.m(covs.list[[i]][in_focal, , drop = FALSE],
                w = sw[in_focal])
      else
        means[[i]]
    }

    targeted[[i]] <- !is.na(targets.list[[i]])

    names(tols[[i]]) <- colnames(covs.list[[i]])

    tols.list[[i]] <- abs(tols[[i]])

    if (treat.types[i] == "cat") {
      balanced[[i]] <- !targeted[[i]]

      treat.sds[[i]] <- NA_real_
      treat.means[[i]] <- NA_real_

      vars.to.standardize <- rep_with(FALSE, tols.list[[i]])
      if (std.binary) vars.to.standardize[bin.covs.list[[i]]] <- TRUE
      if (std.cont) vars.to.standardize[!bin.covs.list[[i]]] <- TRUE

      to_std <- which(vars.to.standardize & !check_if_zero(sds[[i]]))
    }
    else {
      balanced[[i]] <- rep_with(TRUE, targeted[[i]])

      targets_i <- ifelse(targeted[[i]], targets.list[[i]], means[[i]])

      covs.list[[i]] <- center(covs.list[[i]], at = targets_i) #center covs at targets (which will be eventual means)

      treat.sds[[i]] <- sqrt(col.w.v(treat.list[[i]], w = sw))
      treat.means[[i]] <- w.m(treat.list[[i]], w = sw)

      treat.list[[i]] <- (treat.list[[i]] - treat.means[[i]]) / treat.sds[[i]]

      to_std <- !check_if_zero(sds[[i]])
    }

    if (is_not_null(to_std)) {
      covs.list[[i]][, to_std] <- mat_div(covs.list[[i]][, to_std, drop = FALSE], sds[[i]][to_std])
      targets.list[[i]][to_std] <- targets.list[[i]][to_std] / sds[[i]][to_std]
    }

    constraint_df[["constraint"]][constraint_df[["time"]] == i] <- list(
      range_w = if (i == 1L) constraint_range_w(sw, min.w, focal, treat.list[[i]]),
      mean_w = switch(treat.types[i],
                      cat = constraint_mean_w_cat(treat.list[[i]], unique.treats[[i]], sw, n[[i]]),
                      cont = constraint_mean_w_cont(sw, n[[i]])),
      balance = switch(treat.types[i],
                       cat = constraint_balance_cat(covs.list[[i]], treat.list[[i]], sw, tols.list[[i]],
                                                    balanced[[i]], unique.treats[[i]], n[[i]]),
                       cont = constraint_balance_cont(covs.list[[i]], treat.list[[i]], sw, tols.list[[i]],
                                                      balanced[[i]], corr.type)),
      target = switch(treat.types[i],
                      cat = constraint_target_cat(covs.list[[i]], treat.list[[i]], sw, targets.list[[i]],
                                                  tols[[i]], targeted[[i]], unique.treats[[i]], n[[i]], focal),
                      cont = constraint_target_cont(covs.list[[i]], treat.list[[i]], sw, targeted[[i]]))
    )
  }

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
  balance_indices <- unlist(lapply(which(constraint_df$type == "balance"), function(i) {
    seq(constraint_df$nc_cum[i] - constraint_df$nc[i], constraint_df$nc_cum[i])[-1L]
  }))

  target_indices <- unlist(lapply(which(constraint_df$type == "target"), function(i) {
    seq(constraint_df$nc_cum[i] - constraint_df$nc[i], constraint_df$nc_cum[i])[-1L]
  }))

  w <- out$x[seq_len(N)]

  if (abs(min.w) < .Machine$double.eps) {
    w[abs(w) < .Machine$double.eps] <- 0
  }

  w[w < min.w] <- min.w

  #Duals
  duals <- make_list(length(times))

  for (i in times) {
    td <- bd <- NULL

    ti <- which(constraint_df[["time"]] == i & constraint_df[["type"]] == "target")[1L]
    if (constraint_df[["nc"]][ti] > 0) {
      cons <- constraint_df[["constraint"]][[ti]][-(1:3)]
      target_ind <- seq(constraint_df$nc_cum[ti] - constraint_df$nc[ti],
                        constraint_df$nc_cum[ti])[-1L]

      td <- data.frame(constraint = "target",
                       cov = if_null_then(cons$covs, NA_character_),
                       treat = if_null_then(cons$treat, NA_character_),
                       dual = abs(out$y[target_ind]))
    }

    bi <- which(constraint_df[["time"]] == i & constraint_df[["type"]] == "balance")[1L]
    if (constraint_df[["nc"]][bi] > 0) {
      cons <- constraint_df[["constraint"]][[bi]][-(1:3)]
      balance_ind <- seq(constraint_df$nc_cum[bi] - constraint_df$nc[bi],
                         constraint_df$nc_cum[bi])[-1L]

      bd <- data.frame(constraint = "balance",
                       cov = if_null_then(cons$covs, NA_character_),
                       treat = if_null_then(cons$treat.comb, NA_character_),
                       dual = abs(out$y[balance_ind]))
    }

    if (is_not_null(td) || is_not_null(bd)) {
      duals[[i]] <- rbind(td, bd)
    }
  }

  opt_out <- list(w = w,
                  duals = duals,
                  info = out$info,
                  out = out,
                  constraint_df = constraint_df)
  class(opt_out) <- "optweight.fit"

  opt_out
}
