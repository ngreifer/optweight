#' Estimate Stable Balancing Weights
#'
#' Estimate stable balancing weights for treatments and covariates specified in `formula`. The degree of balance for each covariate is specified by `tols` and the target population can be specified with `targets` or `estimand`. See Zubizarreta (2015) and Wang & Zubizarreta (2019) for details of the properties of the weights and the methods used to fit them.
#'
#' @inheritParams optweight.fit
#' @inheritDotParams optweight.fit min.w std.binary std.cont
#' @inheritDotParams optweightMV.fit
#' @param formula A formula with a treatment variable on the left hand side and the covariates to be balanced on the right hand side, or a list thereof. See [glm()] for more details. Interactions and functions of covariates are allowed.
#' @param formula.list A list of formulas, each with a treatment variable on the left hand side and the covariates to be balanced on the right hand side.
#' @param data An optional data set in the form of a data frame that contains the variables in `formula`.
#' @param tols A vector of balance tolerance values for each covariate, or a list thereof. The resulting weighted balance statistics will be at least as small as these values. If only one value is supplied, it will be applied to all covariates. Can also be the output of a call to [process_tols()]. See Details.
#' @param tols.list A list of vectors of balance tolerance values for each covariate for each treatment. The resulting weighted balance statistics will be at least as small as these values. If only one value is supplied, it will be applied to all covariates. See Details.
#' @param estimand The desired estimand, which determines the target population. For binary treatments, can be "ATE", "ATT", "ATC", or `NULL`. For multi-category treatments, can be "ATE", "ATT", or `NULL`. For continuous treatments, can be "ATE" or `NULL`. The default for both is "ATE". For `optweightMV()`, only "ATE" or `NULL` are supported. `estimand` is ignored when `targets` is non-`NULL`. If both `estimand` and `targets` are `NULL`, no targeting will take place. See Details.
#' @param targets A vector of target population mean values for each baseline covariate. The resulting weights will yield sample means within `tols`/2 units of the target values for each covariate. If `NULL` or all `NA`, `estimand` will be used to determine targets. Otherwise, `estimand` is ignored. If any target values are `NA`, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance. Can also be the output of a call to [process_targets()]. See Details.
#' @param s.weights A vector of sampling weights or the name of a variable in `data` that contains sampling weights.
#' @param b.weights A vector of base weights or the name of a variable in `data` that contains base weights. If supplied, the desired norm of the distance between the estimated weights and the base weights is minimized.
#' @param norm `character`; a string containing the name of the norm corresponding to the objective function to minimize. Allowable options include `"l1"` for the L1 norm, `"l2"` for the L2 norm (the default), `"linf"` for the L\eqn{\infty} norm, `"entropy"` for the negative entropy, and `"log"` for the sum of the logs. See [optweight.fit()] for details.
#' @param focal When multi-category treatments are used and `estimand = "ATT"`, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group. If specified, `estimand` will automatically be set to `"ATT"`.
#'
#' @returns
#' For `optweight()`, an `optweight` object with the following elements:
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
#' \item{info}{Information about the performance of the optimization at termination.}
#'
#' For `optweightMV()`, an `optweightMV` object with the following elements:
#' \item{weights}{The estimated weights, one for each unit.}
#' \item{treat.list}{A list of the values of the treatment variables.}
#' \item{covs.list}{A list of the covariates for each treatment used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
#' \item{s.weights}{The provided sampling weights.}
#' \item{b.weights}{The provided base weights.}
#' \item{call}{The function call.}
#' \item{tols}{A list of tolerance values for each covariate for each treatment.}
#' \item{duals}{A list of data.frames containing the dual variables for each covariate for each treatment. See Details for interpretation of these values.}
#' \item{info}{Information about the performance of the optimization at termination.}
#'
#' @details
#' The optimization is performed by the lower-level function [optweight.fit()] (for `optweight()`) or [optweightMV.fit()] (for `optweightMV()`).
#'
#' For binary and multi-category treatments, weights are estimated so that the weighted mean differences of the covariates are within the given tolerance thresholds (unless `std.binary` or `std.cont` are `TRUE`, in which case standardized mean differences are considered for binary and continuous variables, respectively). For a covariate \eqn{x} with specified tolerance \eqn{\delta}, the weighted means of each each group will be within \eqn{\delta} of each other. Additionally, when the ATE is specified as the estimand or a target population is specified, the weighted means of each group will each be within \eqn{\delta/2} of the target means; this ensures generalizability to the same population from which the original sample was drawn.
#'
#' If standardized tolerance values are requested, the standardization factor corresponds to the estimand requested: when the ATE is requested or a target population specified, the standardization factor is the square root of the average variance for that covariate across treatment groups, and when the ATT or ATC are requested, the standardization factor is the standard deviation of the covariate in the focal group. The standardization factor is computed accounting for `s.weights`.
#'
#' For continuous treatments, weights are estimated so that the weighted correlation between the treatment and each covariate is within the specified tolerance threshold. If the ATE is requested or a target population is specified, the means of the weighted covariates and treatment are restricted to be equal to those of the target population to ensure generalizability to the desired target population. The weighted correlation is computed as the weighted covariance divided by the product of the *unweighted* standard deviations. The means used to center the variables in computing the covariance are those specified in the target population.
#'
#' ## Dual Variables
#'
#' Two types of constraints may be associated with each covariate: target constraints and balance constraints. Target constraints require the mean of the covariate to be at (or near) a specific target value in each treatment group (or for the whole group when treatment is continuous). Balance constraints require the means of the covariate in pairs of treatments to be near each other. For binary and multi-category treatments, balance constraints are redundant if target constraints are provided for a variable. For continuous variables, balance constraints refer to the correlation between treatment and the covariate and are not redundant with target constraints. In the `duals` component of the output, each covariate has a dual variable for each nonredundant constraint placed on it.
#'
#'   The dual variable for each constraint is the instantaneous rate of change of the objective function at the optimum corresponding to a change in the constraint. Because this relationship is not linear, large changes in the constraint will not exactly map onto corresponding changes in the objective function at the optimum, but will be close for small changes in the constraint. For example, for a covariate with a balance constraint of .01 and a corresponding dual variable of 40, increasing (i.e., relaxing) the constraint to .025 will decrease the value of the objective function at the optimum by approximately \eqn{(.025 - .01) * 40 = .6}.
#'
#'   For factor variables, `optweight()` takes the sum of the absolute dual variables for the constraints for all levels and reports it as the the single dual variable for the variable itself. This summed dual variable works the same way as dual variables for continuous variables do.
#'
#'
#' @references
#' Chattopadhyay, A., Cohn, E. R., & Zubizarreta, J. R. (2024). One-Step Weighting to Generalize and Transport Treatment Effect Estimates to a Target Population. *The American Statistician*, 78(3), 280–289. \doi{10.1080/00031305.2023.2267598}
#'
#' Källberg, D., & Waernbaum, I. (2023). Large Sample Properties of Entropy Balancing Estimators of Average Causal Effects. *Econometrics and Statistics*. \doi{10.1016/j.ecosta.2023.11.004}
#'
#' Wang, Y., & Zubizarreta, J. R. (2020). Minimal dispersion approximately balancing weights: Asymptotic properties and practical considerations. *Biometrika*, 107(1), 93–105. \doi{10.1093/biomet/asz050}
#'
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' @seealso
#' [optweight.fit()], the lower-level function that performs the fitting. Links on that page can help with diagnosing and fixing more subtle issues with the optimization.
#'
#' \CRANpkg{sbw}, which was the inspiration for this package and provides some additional functionality for binary treatments.
#'
#' \CRANpkg{WeightIt}, which provides a simplified interface to `optweight()` and a more efficient implementation of entropy balancing.
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' # Balancing covariates between treatment groups (binary)
#' (ow1 <- optweight(treat ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   tols = c(.01, .02, .03, .04, .05),
#'                   estimand = "ATE"))
#' bal.tab(ow1)
#'
#' # Exactly alancing covariates with respect to race (multi-category)
#' (ow2 <- optweight(race ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   tols = 0, estimand = "ATT",
#'                   focal = "black"))
#' bal.tab(ow2)
#'
#' # Balancing covariates between treatment groups (binary)
#' # and requesting a specified target population
#' targets <- process_targets(~ age + educ + married +
#'                              nodegree + re74,
#'                            data = lalonde,
#'                            targets = c(26, 12, .4, .5,
#'                                        1000))
#'
#' (ow3a <- optweight(treat ~ age + educ + married +
#'                      nodegree + re74, data = lalonde,
#'                    targets = targets,
#'                    estimand = NULL))
#'
#' bal.tab(ow3a, disp.means = TRUE)
#'
#' # Balancing covariates between treatment groups (binary)
#' # and not requesting a target population
#' (ow3b <- optweight(treat ~ age + educ + married +
#'                      nodegree + re74, data = lalonde,
#'                    targets = NULL,
#'                    estimand = NULL))
#'
#' bal.tab(ow3b, disp.means = TRUE)
#'
#' # Balancing two treatments
#' (ow4 <- optweightMV(list(treat ~ age + educ + race + re74,
#'                          re75 ~ age + educ + race + re74),
#'                     data = lalonde))
#'
#' summary(ow4)
#'
#' bal.tab(ow4)
#'
#' # Using a different norm
#' (ow1b <- optweight(treat ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   tols = c(.01, .02, .03, .04, .05),
#'                   estimand = "ATE",
#'                   norm = "l1"))
#'
#' summary(ow1b, weight.range = FALSE)
#' summary(ow1, weight.range = FALSE)
#'
#' # Allowing for negative weights
#' ow5 <- optweight(treat ~ age + educ + married + race +
#'                    nodegree + re74 + re75,
#'                  data = lalonde,
#'                  estimand = "ATE",
#'                  min.w = -Inf)
#'
#' summary(ow5)

#' @export
optweight <- function(formula, data = NULL, tols = 0, estimand = "ATE",
                      targets = NULL, s.weights = NULL, b.weights = NULL, focal = NULL,
                      norm = "l2", verbose = FALSE, ...) {

  mcall <- match.call()

  chk::chk_string(norm)
  norm <- tolower(norm)
  norm <- match_arg(norm, allowable_norms())

  #Process treat and covs from formula and data
  t.c <- get_covs_and_treat_from_formula2(formula, data, sep = "_")
  simple.covs <- t.c[["simple.covs"]]
  reported.covs <- t.c[["reported.covs"]]

  covs <- t.c[["model.covs"]]
  treat <- t.c[["treat"]]

  #Get treat type
  treat <- assign_treat_type(treat)
  treat.type <- attr(treat, "treat.type")

  if (is_null(covs)) {
    .err("no covariates were specified")
  }

  if (is_null(treat)) {
    .err("no treatment variable was specified")
  }

  if (anyNA(treat) || !all(is.finite(treat))) {
    .err("no missing or non-finite values are allowed in the treatment variable")
  }

  #Process estimand and focal
  f.e.r <- process_focal_and_estimand_w_targets(focal, estimand, targets, treat)
  focal <- f.e.r[["focal"]]
  estimand <- f.e.r[["estimand"]]
  reported.estimand <- f.e.r[["reported.estimand"]]

  check_missing_covs(reported.covs)

  #Process s.weights
  sw <- process_s.weights(s.weights, data)

  #Process b.weights
  bw <- process_b.weights(b.weights, data)

  #Process tols
  tols <- .process_tols_internal(covs, tols, reported.covs)

  #Process targets
  if (is_null(estimand) || is_not_null(targets)) {
    if (is_null(estimand) && is_null(targets)) {
      targets <- NA_real_
    }
    else if (is_not_null(estimand) && is_not_null(targets)) {
      .wrn("`targets` are not `NULL`; ignoring `estimand`")
      estimand <- NULL
    }

    targets <- .process_targets_internal(covs, targets, sw, reported.covs)
  }

  ###Run optweight.fit
  fit_out <- optweight.fit(covs = covs,
                           treat = treat,
                           tols = tols,
                           estimand = estimand,
                           focal = focal,
                           targets = targets,
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
              treat = treat,
              covs = simple.covs,
              s.weights = sw,
              b.weights = bw,
              estimand = switch(treat.type, continuous = NULL, reported.estimand),
              focal = focal,
              norm = norm,
              call = mcall,
              tols = tols,
              duals = duals,
              info = fit_out$info)

  class(out) <- "optweight"

  out
}

#' @export
#' @rdname optweight
optweightMV <- function(formula.list, data = NULL, tols.list = list(0), estimand = "ATE",
                        targets = NULL, s.weights = NULL, b.weights = NULL, focal = NULL,
                        norm = "l2", verbose = FALSE, ...) {

  mcall <- match.call()

  formula.list <- {
    if (is.list(formula.list)) formula.list
    else list(formula.list)
  }

  focal <- NULL

  if (is_not_null(estimand)) {
    estimand <- toupper(estimand)

    if (estimand != "ATE") {
      .err("the only estimand allowed with multivariate treatments is the ATE")
    }
  }

  times <- seq_along(formula.list)

  chk::chk_string(norm)
  norm <- tolower(norm)
  norm <- match_arg(norm, allowable_norms())

  reported.covs.list <- simple.covs.list <- covs.list <- treat.list <- make_list(length(formula.list))
  treat.names <- rep_with(NA_character_, times)
  n <- rep_with(NA_integer_, times)

  for (i in times) {
    #Process treat and covs from formula and data
    t.c <- get_covs_and_treat_from_formula2(formula.list[[i]], data, sep = "_")
    simple.covs.list[[i]] <- t.c[["simple.covs"]]
    reported.covs.list[[i]] <- t.c[["reported.covs"]]

    covs.list[[i]] <- t.c[["model.covs"]]
    treat.list[[i]] <- t.c[["treat"]]

    #Get treat type
    treat.list[[i]] <- assign_treat_type(treat.list[[i]])

    if (is_null(covs.list[[i]])) {
      .err(sprintf("no covariates were specified in the %s formula", ordinal(i)))
    }

    if (is_null(treat.list[[i]])) {
      .err(sprintf("no treatment variable was specified in the %s formula", ordinal(i)))
    }

    treat.names[i] <- if_null_then(attr(treat.list[[i]], "treat.name"),
                                   sprintf("treatment %s", i))

    if (anyNA(treat.list[[i]]) || !all(is.finite(treat.list[[i]]))) {
      .err(sprintf("no missing or non-finite values are allowed in the treatment variable. Missing or non-finite values were found in %s",
                   treat.names[i]))
    }

    check_missing_covs(reported.covs.list[[i]])

    n[i] <- length(treat.list[[i]])
  }

  if (!all_the_same(n)) {
    .err("the same number of units must be present for each treatment")
  }

  #Process s.weights
  sw <- process_s.weights(s.weights, data)

  #Process b.weights
  bw <- process_b.weights(b.weights, data)

  #Process tols
  tols.list <- {
    if (is.list(tols.list)) tols.list
    else list(tols.list)
  }

  if (length(tols.list) == 1L) {
    tols.list <- tols.list[rep_with(1L, times)]
  }

  for (i in times) {
    tryCatch({
      tols.list[[i]] <- .process_tols_internal(covs.list[[i]], tols.list[[i]],
                                               reported.covs.list[[i]],
                                               tols_found_in = "formula.list")
    },
    error = function(e) {
      .err(sprintf("For %s, %s",
                   treat.names[i], conditionMessage(e)),
           tidy = FALSE)
    })
  }

  #Process targets
  if (is_null(estimand) || is_not_null(targets)) {
    if (is_null(estimand) && is_null(targets)) {
      targets <- NA_real_
    }
    else if (is_not_null(estimand) && is_not_null(targets)) {
      .wrn("`targets` are not `NULL`; ignoring `estimand`")
      estimand <- NULL
    }
    else {
      chk::chk_string(estimand)
      estimand <- toupper(estimand)

      if (estimand != "ATE") {
        .err(sprintf("`estimand` cannot be %s with multivariate treatments",
                     add_quotes(estimand)))
      }
    }

    targets <- .process_targets_internal(cbind_distinct(covs.list), targets = targets, sw = sw,
                                         cbind_distinct(reported.covs.list),
                                         targets_found_in = "formula.list")
  }

  ###Run optweight.fit
  fit_out <- optweightMV.fit(treat.list = treat.list,
                             covs.list = covs.list,
                             tols.list = tols.list,
                             estimand = estimand,
                             focal = focal,
                             targets = targets,
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
  duals <- lapply(times, function(i) process_duals(fit_out$duals[[i]], tols.list[[i]]))

  out <- list(weights = fit_out$w,
              treat.list = treat.list,
              covs.list = simple.covs.list,
              s.weights = sw,
              b.weights = bw,
              norm = norm,
              call = mcall,
              tols = tols.list,
              duals = duals,
              info = fit_out$info)

  class(out) <- c("optweightMV", "optweight")

  out
}

#' @exportS3Method print optweight
print.optweight <- function(x, ...) {
  treat.type <- attr(x[["treat"]], "treat.type")

  if (is_not_null(treat.type)) {
    treat.type[treat.type == "multinomial"] <- "multi-category"
  }

  cat(sprintf("An %s object\n", class(x)[1L]))

  cat(sprintf(" - number of obs.: %s\n",
              length(x[["weights"]])))

  cat(sprintf(" - norm minimized: %s\n",
              add_quotes(x[["norm"]])))

  cat(sprintf(" - sampling weights: %s\n",
              if (is_not_null(x[["s.weights"]]) && all_the_same(x[["s.weights"]])) "none" else "present"))

  cat(sprintf(" - base weights: %s\n",
              if (is_not_null(x[["b.weights"]]) && all_the_same(x[["b.weights"]])) "none" else "present"))

  if (is_not_null(x[["treat"]])) {
    cat(sprintf(" - treatment: %s\n",
                switch(treat.type,
                       continuous = "continuous",
                       binary = "2-category",
                       sprintf("%s-category (%s)",
                               nunique(x[["treat"]]),
                               toString(levels(x[["treat"]]))))))
  }

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

#' @exportS3Method print optweightMV
print.optweightMV <- function(x, ...) {
  treat.types <- vapply(x[["treat.list"]], attr, character(1L), "treat.type")
  treat.types[treat.types == "multinomial"] <- "multi-category"

  treat.names <- vapply(x[["treat.list"]], attr, character(1L), "treat.name")

  cat(sprintf("An %s object\n", class(x)[1L]))

  cat(sprintf(" - number of obs.: %s\n",
              length(x[["weights"]])))

  cat(sprintf(" - norm minimized: %s\n",
              add_quotes(x[["norm"]])))

  cat(sprintf(" - sampling weights: %s\n",
              if (is_not_null(x[["s.weights"]]) && all_the_same(x[["s.weights"]])) "none" else "present"))

  cat(sprintf(" - base weights: %s\n",
              if (is_not_null(x[["b.weights"]]) && all_the_same(x[["b.weights"]])) "none" else "present"))

  cat(sprintf(" - number of treatments: %s\n%s",
              length(x[["treat.list"]]),
              do.call("paste0", lapply(seq_along(x$covs.list), function(i) {
                sprintf("    %s: %s\n",
                        treat.names[i],
                        switch(treat.types[i],
                               continuous = "continuous",
                               binary = "2-category",
                               sprintf("%s-category (%s)",
                                       nunique(x[["treat.list"]][[i]]),
                                       toString(levels(x[["treat.list"]][[i]])))))
              }))))

  cat(sprintf(" - covariates: \n%s",
              do.call("paste0", lapply(seq_along(x$covs.list), function(i) {
                sprintf("    + for %s: %s\n",
                        treat.names[i],
                        if (length(names(x[["covs.list"]][[i]])) > 60L) "too many to name"
                        else toString(names(x[["covs.list"]][[i]])))
              }))))

  invisible(x)
}
