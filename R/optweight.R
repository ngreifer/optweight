#' Estimate Stable Balancing Weights
#'
#' Estimate stable balancing weights for treatments and covariates specified in `formula`. The degree of balance for each covariate is specified by `tols` and the target population can be specified with `targets` or `estimand`. See Zubizarreta (2015), Wang & Zubizarreta (2019), and Yiu & Su (2018) for details of the properties of the weights and the methods used to fit them.
#'
#' @inheritParams optweight.fit
#' @inheritDotParams optweight.fit norm min.w force std.binary std.cont
#' @param formula A formula with a treatment variable on the left hand side and the covariates to be balanced on the right hand side, or a list thereof. See [glm()] for more details. Interactions and functions of covariates are allowed.
#' @param data An optional data set in the form of a data frame that contains the variables in `formula`.
#' @param tols A vector of balance tolerance values for each covariate, or a list thereof. The resulting weighted balance statistics will be at least as small as these values. If only one value is supplied, it will be applied to all covariates. Can also be the output of a call to [check_tols()] for point treatments. See Details.
#' @param estimand The desired estimand, which determines the target population. For binary treatments, can be "ATE", "ATT", "ATC", or `NULL`. For multi-category treatments, can be "ATE", "ATT", or `NULL`. For continuous treatments, can be "ATE" or `NULL`. The default for both is "ATE". For longitudinal treatments, only "ATE" is supported. `estimand` is ignored when `targets` is non-`NULL`. If both `estimand` and `targets` are `NULL`, no targeting will take place. See Details.
#' @param targets A vector of target population mean values for each baseline covariate. The resulting weights will yield sample means within `tols`/2 units of the target values for each covariate. If `NULL` or all `NA`, `estimand` will be used to determine targets. Otherwise, `estimand` is ignored. If any target values are `NA`, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance. Can also be the output of a call to [check_targets()]. See Details.
#' @param s.weights A vector of sampling weights or the name of a variable in `data` that contains sampling weights. Optimization occurs on the product of the sampling weights and the estimated weights.
#' @param b.weights A vector of base weights or the name of a variable in `data` that contains base weights. If supplied, the desired norm of the distance between the estimated weights and the base weights is minimized.
#' @param focal When multi-category treatments are used and the "ATT" is requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group. If specified, `estimand` will automatically be set to `"ATT"`.
#'
#' @returns
#' If only one time point is specified, an `optweight` object with the following elements:
#' \item{weights}{The estimated weights, one for each unit.}
#' \item{treat}{The values of the treatment variable.}
#' \item{covs}{The covariates used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
#' \item{s.weights}{The provided sampling weights.}
#' \item{b.weights}{The provided base weights.}
#' \item{estimand}{The estimand requested.}
#' \item{focal}{The focal variable if the ATT was requested with a multi-category treatment.}
#' \item{call}{The function call.}
#' \item{tols}{The tolerance values for each covariate.}
#' \item{duals}{A data.frame containing the dual variables for each covariate. See Details for interpretation of these values.}
#' \item{info}{The `info` component of the output of [osqp::solve_osqp()], which contains information on the performance of the optimization at termination.}
#'
#' Otherwise, if multiple time points are specified, an `optweightMSM` object with the following elements:
#' \item{weights}{The estimated weights, one for each unit.}
#' \item{treat.list}{A list of the values of the treatment variables at each time point.}
#' \item{covs.list}{A list of the covariates at each time point used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
#' \item{s.weights}{The provided sampling weights.}
#' \item{b.weights}{The provided base weights.}
#' \item{call}{The function call.}
#' \item{tols}{A list of tolerance values for each covariate at each time point.}
#' \item{duals}{A list of data.frames containing the dual variables for each covariate at each time point. See Details for interpretation of these values.}
#' \item{info}{The `info` component of the output of [osqp::solve_osqp()], which contains information on the performance of the optimization at termination.}
#'
#' @details
#' The optimization is performed by the lower-level function [optweight.fit()] using [osqp::solve_osqp()] in the \pkg{osqp} package, which provides a straightforward interface to specifying the constraints and objective function for quadratic optimization problems and uses a fast and flexible solving algorithm.
#'
#' For binary and multi-category treatments, weights are estimated so that the weighted mean differences of the covariates are within the given tolerance thresholds (unless `std.binary` or `std.cont` are `TRUE`, in which case standardized mean differences are considered for binary and continuous variables, respectively). For a covariate \eqn{x} with specified tolerance \eqn{\delta}, the weighted means of each each group will be within \eqn{\delta} of each other. Additionally, when the ATE is specified as the estimand or a target population is specified, the weighted means of each group will each be within \eqn{\delta/2} of the target means; this ensures generalizability to the same population from which the original sample was drawn.
#'
#' If standardized tolerance values are requested, the standardization factor corresponds to the estimand requested: when the ATE is requested or a target population specified, the standardization factor is the square root of the average variance for that covariate across treatment groups, and when the ATT or ATC are requested, the standardization factor is the standard deviation of the covariate in the focal group. The standardization factor is always unweighted.
#'
#' For continuous treatments, weights are estimated so that the weighted correlation between the treatment and each covariate is within the specified tolerance threshold. If the ATE is requested or a target population is specified, the means of the weighted covariates and treatment are restricted to be equal to those of the target population to ensure generalizability to the desired target population. The weighted correlation is computed as the weighted covariance divided by the product of the *unweighted* standard deviations. The means used to center the variables in computing the covariance are those specified in the target population.
#'
#' For longitudinal treatments, only "wide" data sets, where each row corresponds to a unit's entire variable history, are supported. You can use [reshape()] or other functions to transform your data into this format; see example in the documentation for `weightitMSM()` in the \pkg{WeightIt} package. Currently, longitudinal treatments are not recommended as the use as stable balancing weights with them has not been validated.
#'
#' ## Dual Variables
#'
#' Two types of constraints may be associated with each covariate: target constraints and balance constraints. Target constraints require the mean of the covariate to be at (or near) a specific target value in each treatment group (or for the whole group when treatment is continuous). Balance constraints require the means of the covariate in pairs of treatments to be near each other. For binary and multi-category treatments, balance constraints are redundant if target constraints are provided for a variable. For continuous variables, balance constraints refer to the correlation between treatment and the covariate and are not redundant with target constraints. In the `duals` component of the output, each covariate has a dual variable for each nonredundant constraint placed on it.
#'
#'   The dual variable for each constraint is the instantaneous rate of change of the objective function at the optimum due to a change in the constraint. Because this relationship is not linear, large changes in the constraint will not exactly map onto corresponding changes in the objective function at the optimum, but will be close for small changes in the constraint. For example, for a covariate with a balance constraint of .01 and a corresponding dual variable of .4, increasing (i.e., relaxing) the constraint to .025 will decrease the value of the objective function at the optimum by approximately \eqn{(.025 - .01) * .4 = .006}. When the L2 norm is used, this change corresponds to a change in the variance of the weights, which directly affects the effective sample size (though the magnitude of this effect depends on the original value of the effective sample size).
#'
#'   For factor variables, `optweight()` takes the sum of the absolute dual variables for the constraints for all levels and reports it as the the single dual variable for the variable itself. This summed dual variable works the same way as dual variables for continuous variables do.
#'
#' ## Solving Convergence Failure
#'
#'   Sometimes the optimization will fail to converge at a solution. There are a variety of reasons why this might happen, which include that the constraints are nearly impossible to satisfy or that the optimization surface is relatively flat. It can be hard to know the exact cause or how to solve it, but this section offers some solutions one might try.
#'
#'   Rarely is the problem too few iterations, though this is possible. Most problems can be solved in the default 200,000 iterations, but sometimes it can help to increase this number with the `max_iter` argument. Usually, though, this just ends up taking more time without a solution found.
#'
#'   If the problem is that the constraints are too tight, it can be helpful to loosen the constraints. Sometimes examining the dual variables of a solution that has failed to converge can reveal which constraints are causing the problem.
#'
#'   Sometimes a suboptimal solution is possible; such a solution does not satisfy the constraints exactly but will come pretty close. To allow these solutions, the arguments `eps_abs` and `eps_rel` can be increased from `1e-8` to larger values. These should be adjusted together since they both must be satisfied for convergence to occur; this can be done easily using the shortcut argument `eps`, which changes both `eps_abs` and `eps_rel` to the set value.
#'
#'   With continuous treatments, solutions that failed to converge may still be useable. Make sure to assess balance and examine the weights even after a optimal solution is not found, because the solution that is found may be good enough.
#'
#' @references
#' Wang, Y., & Zubizarreta, J. R. (2020). Minimal dispersion approximately balancing weights: Asymptotic properties and practical considerations. *Biometrika*, 107(1), 93–105. \doi{10.1093/biomet/asz050}
#'
#' Yiu, S., & Su, L. (2018). Covariate association eliminating weights: a unified weighting framework for causal effect estimation. *Biometrika*. \doi{10.1093/biomet/asy015}
#'
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' @seealso
#' The OSQP [docs](https://osqp.org/docs/index.html) for more information on \pkg{osqp}, the underlying solver, and the options for [osqp::solve_osqp()]. [osqp::osqpSettings()] for details on options for `solve_osqp()`.
#'
#' [optweight.fit()], the lower-level function that performs the fitting.
#'
#' \CRANpkg{sbw}, which was the inspiration for this package and provides some additional functionality for binary treatments.
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (ow1 <- optweight(treat ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   tols = c(.01, .02, .03, .04, .05),
#'                   estimand = "ATE"))
#' bal.tab(ow1)
#'
#' #Exactly alancing covariates with respect to race (multi-category)
#' (ow2 <- optweight(race ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   tols = 0, estimand = "ATT", focal = "black"))
#' bal.tab(ow2)
#'
#' # #Balancing covariates with longitudinal treatments
#' # #NOT VALID; DO NOT DO THIS.
#' # library("twang")
#' # data("iptwExWide")
#' #
#' # ##Weighting more recent covariates more strictly
#' # (ow3 <- optweight(list(tx1 ~ use0 + gender + age,
#' #                        tx2 ~ tx1 + use1 + use0  + gender +
#' #                          age,
#' #                        tx3 ~ tx2 + use2 + tx1 + use1 +
#' #                          use0 + gender + age),
#' #                   data = iptwExWide,
#' #                   tols = list(c(.001, .001, .001),
#' #                               c(.001, .001, .01, .01, .01),
#' #                               c(.001, .001, .01, .01,
#' #                                 .1, .1, .1))))
#' # bal.tab(ow3)
#'
#' #Balancing covariates between treatment groups (binary)
#' #and requesting a specified target population
#' (ow4a <- optweight(treat ~ age + educ + married +
#'                      nodegree + re74, data = lalonde,
#'                    tols = 0,
#'                    targets = c(26, 12, .4, .5, 1000),
#'                    estimand = NULL))
#' bal.tab(ow4a, disp.means = TRUE)
#'
#' #Balancing covariates between treatment groups (binary)
#' #and not requesting a target population
#' (ow4b <- optweight(treat ~ age + educ + married +
#'                      nodegree + re74, data = lalonde,
#'                    tols = 0,
#'                    targets = NULL,
#'                    estimand = NULL))
#' bal.tab(ow4b, disp.means = TRUE)
#'

#' @export
optweight <- function(formula, data = NULL, tols = 0, estimand = "ATE",
                      targets = NULL, s.weights = NULL, b.weights = NULL, focal = NULL,
                      verbose = FALSE, ...) {

  call <- match.call()

  formula.list <- {
    if (is.list(formula)) formula
    else list(formula)
  }

  times <- seq_along(formula.list)
  onetime <- length(times) == 1L

  if (!onetime) {
    force <- ...get("force", FALSE)
    chk::chk_flag(force)
    if (!force) {
      .err("optweights are currently not valid for longitudinal treatments. Set `force = TRUE` to bypass this message at your own risk")
    }
  }

  tols.list <- {
    if (is.list(tols)) tols
    else list(tols)
  }

  if (length(tols.list) == 1L) {
    tols.list <- replicate(length(times), tols.list[[1L]], simplify = FALSE)
  }

  #Process targets
  targets <- check_targets(formula.list[[1L]], data, targets, stop = TRUE)

  if (isTRUE(attr(targets, "ATE"))) {
    estimand <- "ATE"
    targets <- rep_with(NA_real_, targets)
  }

  reported.covs.list <- simple.covs.list <- covs.list <- treat.list <- ct <- make_list(length(formula.list))
  n <- rep_with(NA_integer_, times)

  for (i in times) {
    #Process treat and covs from formula and data
    t.c <- get_covs_and_treat_from_formula2(formula.list[[i]], data)
    simple.covs.list[[i]] <- t.c[["simple.covs"]]
    reported.covs.list[[i]] <- t.c[["reported.covs"]]

    covs.list[[i]] <- t.c[["model.covs"]]
    treat.list[[i]] <- t.c[["treat"]]
    #treat.name <- t.c[["treat.name"]]

    #Get treat type
    treat.list[[i]] <- assign_treat_type(treat.list[[i]])
    treat.type <- attr(treat.list[[i]], "treat.type")

    if (onetime) {
      if (is_null(covs.list[[i]])) {
        .err("no covariates were specified")
      }

      if (is_null(treat.list[[i]])) {
        .err("no treatment variable was specified")
      }

      if (anyNA(treat.list[[i]]) || !all(is.finite(treat.list[[i]]))) {
        .err("no missing or non-finite values are allowed in the treatment variable")
      }

      #Process estimand and focal
      f.e.r <- process.focal.and.estimand(focal, estimand, targets, treat.list[[i]], treat.type)
      focal <- f.e.r[["focal"]]
      estimand <- f.e.r[["estimand"]]
      reported.estimand <- f.e.r[["reported.estimand"]]
    }
    else {
      if (is_null(covs.list[[i]])) {
        .err(sprintf("no covariates were specified in the %s formula", ordinal(i)))
      }

      if (is_null(treat.list[[i]])) {
        .err(sprintf("no treatment variable was specified in the %s formula", ordinal(i)))
      }

      if (anyNA(treat.list[[i]]) || !all(is.finite(treat.list[[i]]))) {
        .err(sprintf("no missing or non-finite values are allowed in the treatment variable. Missing or non-finite values were found in treatment %s", i))
      }

      if (is_not_null(estimand) && toupper(estimand) != "ATE") {
        .err("the only estimand allowed with multiple treatments is the ATE")
      }

      focal <- NULL
      estimand <- toupper(estimand)
      reported.estimand <- estimand
    }

    check_missing_covs(reported.covs.list[[i]])

    n[i] <- length(treat.list[[i]])

    tryCatch({
      ct[[i]] <- check_tols(formula.list[[i]], data, tols.list[[i]], stop = TRUE)
    },
    error = function(e) {

      if (onetime) {
        .err(conditionMessage(e))
      }

      .err(sprintf("for treatment %s, %s",
                   i, conditionMessage(e)))
    })

    tols.list[[i]] <- attr(ct[[i]], "internal.tols")
  }

  if (!all_the_same(n)) {
    .err("the same number of units must be present in each time point")
  }

  #Process s.weights
  sw <- process.s.weights(s.weights, data)

  #Process b.weights
  bw <- process.b.weights(b.weights, data)

  ###Run optweight.fit
  fit_out <- optweight.fit(treat.list = treat.list,
                           covs.list = covs.list,
                           tols = tols.list,
                           estimand = estimand,
                           focal = focal,
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

  warn <- FALSE
  test.w <- {
    if (is_null(sw)) fit_out$w
    else fit_out$w * sw
  }

  if (anyNA(test.w)) {
    .err("some weights are NA, which means something went wrong")
  }

  for (ti in treat.list) {
    warn <- switch(get_treat_type(ti),
                   continuous = rms_dev(test.w) > 4,
                   any_apply(unique(ti), function(x) rms_dev(test.w[ti == x]) > 4)
    )

    if (warn) {
      .wrn("some extreme weights were generated. Examine them with `summary()` and maybe relax the constraints")
      break
    }
  }

  #Process duals
  for (i in times) {
    original.vars <- attr(ct[[i]], "original.vars")
    d <- fit_out$duals[[i]]
    d$cov <- vapply(d$cov, function(c) original.vars[names(original.vars) == c][1L], character(1L))
    d$dual <- ave(d$dual, d$constraint, d$cov, FUN = sum) #Total effect of constraint on obj. fun. is sum of abs(duals)
    fit_out$duals[[i]] <- unique(d[names(d) != "treat"])
    rownames(fit_out$duals[[i]]) <- NULL
  }

  if (onetime) {
    out <- list(weights = fit_out$w,
                treat = treat.list[[1L]],
                covs = simple.covs.list[[1L]],
                s.weights = sw,
                b.weights = bw,
                estimand = switch(treat.type, continuous = NULL, reported.estimand),
                focal = focal,
                call = call,
                tols = tols.list[[1L]],
                duals = fit_out$duals[[1L]],
                info = fit_out$info)

    class(out) <- "optweight"
  }
  else {
    out <- list(weights = fit_out$w,
                treat.list = treat.list,
                covs.list = simple.covs.list,
                s.weights = sw,
                b.weights = bw,
                #estimand = reported.estimand,
                call = call,
                tols = tols.list,
                duals = fit_out$duals,
                info = fit_out$info)

    class(out) <- c("optweightMSM", "optweight")
  }

  out
}

#' @exportS3Method print optweight
print.optweight <- function(x, ...) {
  treat.type <- attr(x[["treat"]], "treat.type")

  cat("An optweight object\n")
  cat(sprintf(" - number of obs.: %s\n",
              length(x[["weights"]])))
  cat(sprintf(" - sampling weights: %s\n",
              if (all_the_same(x[["s.weights"]])) "none" else "present"))
  cat(sprintf(" - treatment: %s\n",
              switch(treat.type,
                     continuous = "continuous",
                     `multi-category` = sprintf("%s-category (%s)",
                                           nunique(x[["treat"]]),
                                           toString(levels(x[["treat"]]))),
                     "2-category")))

  if (is_not_null(x[["estimand"]])) {
    cat(sprintf(" - estimand: %s%s\n",
                x[["estimand"]],
                if (is_not_null(x[["focal"]])) sprintf(" (focal: %s)", x[["focal"]]) else ""))
  }

  cat(sprintf(" - covariates: %s\n",
              if (length(names(x[["covs"]])) > 60L) "too many to name"
              else toString(names(x[["covs"]]))))

  invisible(x)
}

#' @exportS3Method print optweightMSM
print.optweightMSM <- function(x, ...) {
  treat.types <- vapply(x[["treat.list"]], attr, character(1L), "treat.type")

  cat("An optweightMSM object\n")
  cat(sprintf(" - number of obs.: %s\n",
              length(x[["weights"]])))
  cat(sprintf(" - sampling weights: %s\n",
              if (all_the_same(x[["s.weights"]])) "none" else "present"))
  cat(sprintf(" - number of time points: %s\n",
              length(x[["treat.list"]])))

  cat(sprintf(" - treatment: \n%s",
              do.call("paste0", lapply(seq_along(x$covs.list), function(i) {
                sprintf("    + time %s: %s\n",
                        i,
                        switch(treat.types[i],
                               continuous = "continuous",
                               `multi-category` = sprintf("%s-category (%s)",
                                                     nunique(x[["treat.list"]][[i]]),
                                                     toString(levels(x[["treat.list"]][[i]])))),
                        "2-category")
              }))))

  cat(sprintf(" - covariates: \n%s",
              do.call("paste0", lapply(seq_along(x$covs.list), function(i) {
                if (i == 1L) {
                  sprintf("    + baseline: %s\n", toString(names(x$covs.list[[i]])))
                }
                else {
                  sprintf("    + after time %s: %s\n", i - 1, toString(names(x$covs.list[[i]])))
                }
              }))))

  invisible(x)
}
