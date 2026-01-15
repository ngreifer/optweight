#' Stable Balancing Weights
#'
#' Estimates stable balancing weights for the supplied treatments and covariates. The degree of balance for each covariate is specified by `tols` and the target population can be specified with `targets` or `estimand`. See Zubizarreta (2015) and Wang & Zubizarreta (2020) for details of the properties of the weights and the methods used to fit them.
#'
#' @param formula a formula with a treatment variable on the left hand side and the covariates to be balanced on the right hand side, or a list thereof. Interactions and functions of covariates are allowed.
#' @param data an optional data set in the form of a data frame that contains the variables in `formula`.
#' @param tols a vector of balance tolerance values for each covariate. The resulting weighted balance statistics will be at least as small as these values. If only one value is supplied, it will be applied to all covariates. Can also be the output of a call to [process_tols()]. See Details. Default is 0 for all covariates.
#' @param estimand a string containing the desired estimand, which determines the target population. For binary treatments, can be "ATE", "ATT", "ATC", or `NULL`. For multi-category treatments, can be "ATE", "ATT", or `NULL`. For continuous treatments, can be "ATE" or `NULL`. The default for both is "ATE". `estimand` is ignored when `targets` is non-`NULL`. If both `estimand` and `targets` are `NULL`, no targeting will take place. See Details.
#' @param targets an optional vector of target population mean values for each covariate. The resulting weights ensure the midpoint between group means are within `target.tols` units of the target values for each covariate. If `NULL` or all `NA`, `estimand` will be used to determine targets. Otherwise, `estimand` is ignored. If any target values are `NA`, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest value of the objective function; this is only allowed for binary and multi-category treatments. Can also be the output of a call to [process_targets()]. See Details.
#' @param target.tols a vector of target balance tolerance values for each covariate. For binary and multi-category treatments, the average of each pair of means will be at most as far from the target means as these values. Can also be the output of a call to [process_tols()]. See Details. Default is 0 for all covariates. Ignored with continuous treatments and when `estimand` is `"ATT"` or `"ATC"`.
#' @param s.weights a vector of sampling weights. For `optweight()`, can also be the name of a variable in `data` that contains sampling weights.
#' @param b.weights a vector of base weights. If supplied, the desired norm of the distance between the estimated weights and the base weights is minimized. For `optweight()`, can also the name of a variable in `data` that contains base weights.
#' @param norm `character`; a string containing the name of the norm corresponding to the objective function to minimize. Allowable options include `"l1"` for the \eqn{L_1} norm, `"l2"` for the \eqn{L_2} norm (the default), `"linf"` for the \eqn{L_\infty} norm, `"entropy"` for the relative entropy, and `"log"` for the sum of the negative logs. See Details.
#' @param focal when multi-category treatments are used and `estimand = "ATT"`, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group. If specified, `estimand` will automatically be set to `"ATT"`.
#' @param covs a numeric matrix of covariates to be balanced.
#' @param treat a vector of treatment statuses. Non-numeric (i.e., factor or character) vectors are allowed.
#' @param std.binary,std.cont `logical`; whether the tolerances are in standardized mean units (`TRUE`) or raw units (`FALSE`) for binary variables and continuous variables, respectively. The default is `FALSE` for `std.binary` because raw proportion differences make more sense than standardized mean difference for binary variables. These arguments are analogous to the `binary` and `continuous` arguments in \pkgfun{cobalt}{bal.tab}.
#' @param min.w `numeric`; a single value less than 1 for the smallest allowable weight. Some analyses require nonzero weights for all units, so a small, nonzero minimum may be desirable. The default is `1e-8` (\eqn{10^{-8}}), which does not materially change the properties of the weights from a minimum of 0 but prevents warnings in some packages that use weights in model fitting. When `norm` is `"entropy"` or `"log"` and `min.w <= 0`, `min.w` will be set to the smallest nonzero value.
#' @param verbose `logical`; whether information on the optimization problem solution should be printed. Default is `FALSE`.
#' @param solver string; the name of the optimization solver to use. Allowable options depend on `norm`. Default is to use whichever eligible solver is installed, if any, or the default solver for the corresponding `norm`. See Details for information.
#' @param \dots for `optweight()`, additional arguments passed to `optweight.fit()`, including options that are passed to the settings function corresponding to `solver`.
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
#' \item{tols}{The balance tolerance values for each covariate.}
#' \item{target.tols}{The target balance tolerance values for each covariate.}
#' \item{duals}{A data.frame containing the dual variables for each covariate. See Details for interpretation of these values.}
#' \item{info}{A list containing information about the performance of the optimization at termination.}
#' \item{norm}{The `norm` used.}
#' \item{solver}{The `solver` used.}
#'
#' For `optweight.fit()`, an `optweight.fit` object with the following elements:
#' \item{w}{The estimated weights, one for each unit.}
#' \item{duals}{A data.frame containing the dual variables for each covariate.}
#' \item{info}{A list containing information about the performance of the optimization at termination.}
#' \item{norm}{The `norm` used.}
#' \item{solver}{The `solver` used.}
#'
#' @details
#' `optweight()` is the primary user-facing function for estimating stable balancing weights. The optimization is performed by the lower-level function `optweight.fit()`, which transforms the inputs into the required inputs for the optimization functions and then supplies the outputs (the weights, dual variables, and convergence information) back to `optweight()`. Little processing of inputs is performed by `optweight.fit()`, as this is normally handled by `optweight()`.
#'
#' For binary and multi-category treatments, weights are estimated so that the weighted mean differences of the covariates are within the given tolerance thresholds controlled by `tols` and `target.tols` (unless `std.binary` or `std.cont` are `TRUE`, in which case standardized mean differences are considered for binary and continuous variables, respectively). For a covariate \eqn{x} with specified balance tolerance \eqn{\delta} and target tolerance \eqn{\varepsilon}, the weighted means of each each group will be within \eqn{\delta} of each other, and the midpoint between the weighted group means will be with \eqn{\varepsilon} of the target means. More specifically, the constraints are specified as follows:
#' \deqn{
#' \left| \bar{x}^w_1 - \bar{x}^w_0 \right| \le \delta \\
#' \left| \frac{\bar{x}^w_1 + \bar{x}^w_0}{2} - \bar{x}^* \right| \le \varepsilon
#' }
#' where \eqn{\bar{x}^w_1} and \eqn{\bar{x}^w_0} are the weighted means of covariate \eqn{x} for treatment groups 1 and 0, respectively, and \eqn{\bar{x}^*} is the target mean for that covariate. \eqn{\delta} corresponds to `tols`, and \eqn{\varepsilon} corresponds to `target.tols`. Setting a covariate's value of `target.tols` to `Inf` or its `target` to `NA` both serve to remove the second constraint, as is done in Barnard et al. (2025).
#'
#' If standardized tolerance values are requested, the standardization factor corresponds to the estimand requested: when the ATE is requested or a target population specified, the standardization factor is the square root of the average variance for that covariate across treatment groups, and when the ATT or ATC are requested, the standardization factor is the standard deviation of the covariate in the focal group. The standardization factor is computed accounting for `s.weights`.
#'
#' Target and balance constraints are applied to the product of the estimated weights and the sampling weights. In addition, the sum of the product of the estimated weights and the sampling weights is constrained to be equal to the sum of the product of the base weights and sampling weights. For binary and multi-category treatments, these constraints apply within each treatment group.
#'
#' ## Continuous treatments
#'
#' For continuous treatments, weights are estimated so that the weighted correlation between the treatment and each covariate is within the specified tolerance threshold. The means of the weighted covariates and treatment are restricted to be exactly equal to those of the target population to ensure generalizability to the desired target population, regardless of `tols` or `target.tols`. The weighted correlation is computed as the weighted covariance divided by the product of the *unweighted* standard deviations. The means used to center the variables in computing the covariance are those specified in the target population.
#'
#' ## `norm`
#'
#' The objective function for the optimization problem is \eqn{f\left(\mathbf{w}, \mathbf{b},\mathbf{s}\right)}, where \eqn{\mathbf{w}=\{w_1, \dots, w_n\}} are the estimated weights, \eqn{\mathbf{s}=\{s_1, \dots, s_n\}} are sampling weights (supplied by `s.weights`), and \eqn{\mathbf{b}=\{b_1, \dots, b_n\}} are base weights (supplied by `b.weights`). The `norm` argument determines \eqn{f(.,.,.)}, as detailed below:
#'
#' * when `norm = "l2"`, \eqn{f\left(\mathbf{w}, \mathbf{b},\mathbf{s}\right) = \frac{1}{n} \sum_i {s_i(w_i - b_i)^2}}
#' * when `norm = "l1"`, \eqn{f\left(\mathbf{w}, \mathbf{b},\mathbf{s}\right) = \frac{1}{n} \sum_i {s_i \vert w_i - b_i \vert}}
#' * when `norm = "linf"`, \eqn{f\left(\mathbf{w}, \mathbf{b},\mathbf{s}\right) = \max_i {\vert w_i - b_i \vert}}
#' * when `norm = "entropy"`, \eqn{f\left(\mathbf{w}, \mathbf{b},\mathbf{s}\right) = \frac{1}{n} \sum_i {s_i w_i \log \frac{w_i}{b_i}}}
#' * when `norm = "log"`, \eqn{f\left(\mathbf{w}, \mathbf{b},\mathbf{s}\right) = \frac{1}{n} \sum_i {-s_i \log \frac{w_i}{b_i}}}
#'
#' By default, `s.weights` and `b.weights` are set to 1 for all units unless supplied. `b.weights` must be positive when `norm` is `"entropy"` or `"log"`, and `norm = "linf"` cannot be used when `s.weights` are supplied.
#'
#' When `norm = "l2"` and both `s.weights` and `b.weights` are `NULL`, weights are estimated to maximize the effective sample size. When `norm = "entropy"`, the estimated weights are equivalent to entropy balancing weights (Källberg & Waernbaum, 2023). When `norm = "log"`, `b.weights` are ignored in the optimization, as they do not affect the estimated weights.
#'
#' ## Dual Variables
#'
#' Two types of constraints may be associated with each covariate: target constraints and balance constraints, controlled by `target.tols` and `tols`, respectively. In the `duals` component of the output, each covariate has a dual variable for each constraint placed on it. The dual variable for each constraint is the instantaneous rate of change of the objective function at the optimum corresponding to a change in the constraint. Because this relationship is not linear, large changes in the constraint will not exactly map onto corresponding changes in the objective function at the optimum, but will be close for small changes in the constraint. For example, for a covariate with a balance constraint of .01 and a corresponding dual variable of 40, increasing (i.e., relaxing) the constraint to .025 will decrease the value of the objective function at the optimum by approximately \eqn{(.025 - .01) * 40 = .6}.
#'
#' For factor variables, `optweight()` takes the sum of the absolute dual variables for the constraints for all levels and reports it as the the single dual variable for the variable itself. This summed dual variable works the same way as dual variables for continuous variables do.
#'
#' An addition dual variable is computed for the constraint on the range of the weights, controlled by `min.w`. A high dual variable for this constraint implies that decreasing `min.w` will decrease the value of the objective function at the optimum.
#'
#' ## `solver`
#'
#' The `solver` argument controls which optimization solver is used. Different solvers are compatible with each `norm`. See the table below for allowable options, which package they require, which function does the solving, and which function controls the settings.
#'
#' | `solver`     | `norm`                   | Package            | Solver function | Settings function |
#' |--------------|--------------------------|--------------------|-----------------|-------------------|
#' | `"osqp"`     | `"l2"`, `"l1"`, `"linf"` | \CRANpkg{osqp}     | [osqp::solve_osqp()] | [osqp::osqpSettings()] |
#' | `"highs"`    | `"l2"`, `"l1"`, `"linf"` | \CRANpkg{highs}    | \pkgfun{highs}{highs_solve} | \pkgfun{highs}{highs_control} / \pkgfun{highs}{highs_available_solver_options} |
#' | `"lpsolve"`  | `"l1"`, `"linf"`         | \CRANpkg{lpSolve}  | \pkgfun{lpSolve}{lp} | . |
#' | `"scs"`      | `"entropy"`, `"log"`     | \CRANpkg{scs}      | \pkgfun{scs}{scs} | \pkgfun{scs}{scs_control} |
#' | `"clarabel"` | `"entropy"`, `"log"`     | \CRANpkg{clarabel} | \pkgfun{clarabel}{clarabel} | \pkgfun{clarabel}{clarabel_control} |
#'
#' Note that `"lpsolve"` can only be used when `min.w` is nonnegative.
#'
#' The default `solver` for each `norm` is as follows:
#'
#' | `norm`      | Default `solver` |
#' |-------------|------------------|
#' | `"l2"`      | `"osqp"`         |
#' | `"l1"`      | `"highs"`        |
#' | `"linf"`    | `"highs"`        |
#' | `"entropy"` | `"scs"`          |
#' | `"log"`     | `"scs"`          |
#'
#' If the package corresponding to a default `solver` is not installed but the package for a different eligible solver is, that will be used. Otherwise, you will be asked to install the required package. \pkg{osqp} is required for \pkg{optweight}, and so will be the default for the `"l1"` and `"linf"` norms if \pkg{highs} is not installed. The default package is the one has shown good performance for the given norm in informal testing; generally, all eligible solvers perform about equally well in terms of accuracy but differ in time taken.
#'
#' ## Solving Convergence Failure
#'
#' Sometimes the optimization will fail to converge at a solution. There are a variety of reasons why this might happen, which include that the constraints are nearly impossible to satisfy or that the optimization surface is relatively flat. It can be hard to know the exact cause or how to solve it, but this section offers some solutions one might try. Typically, solutions can be found most easily when using the `"l2"` norm; other norms, especially `"linf"` and `"l1"`, are more likely to see problems.
#'
#' Rarely is the problem too few iterations, though this is possible. Most problems can be solved in the default 200,000 iterations, but sometimes it can help to increase this number with the `max_iter` argument. Usually, though, this just ends up taking more time without a solution found.
#'
#' If the problem is that the constraints are too tight, it can be helpful to loosen the constraints. Sometimes examining the dual variables of a solution that has failed to converge can reveal which constraints are causing the problem. An extreme value of a dual variable typically suggests that its corresponding constraint is one cause of the failure to converge.
#'
#' Sometimes a suboptimal solution is possible; such a solution does not satisfy the constraints exactly but will come pretty close. To allow these solutions, the argument `eps` can be increased to larger values. This is more likely to occur when `s.weights` are supplied.
#'
#' Sometimes using a different solver can improve performance. Using the default `solver` for each `norm`, as described above, can reduce the probability of convergence failures.
#'
#' @references
#' Barnard, M., Huling, J. D., & Wolfson, J. (2025). Partially Retargeted Balancing Weights for Causal Effect Estimation Under Positivity Violations (No. arXiv:2510.22072). arXiv. \doi{10.48550/arXiv.2510.22072}
#'
#' Chattopadhyay, A., Cohn, E. R., & Zubizarreta, J. R. (2024). One-Step Weighting to Generalize and Transport Treatment Effect Estimates to a Target Population. *The American Statistician*, 78(3), 280–289. \doi{10.1080/00031305.2023.2267598}
#'
#' de los Angeles Resa, M., & Zubizarreta, J. R. (2020). Direct and Stable Weight Adjustment in Non-Experimental Studies With Multivalued Treatments: Analysis of the Effect of an Earthquake on Post-Traumatic Stress. *Journal of the Royal Statistical Society Series A: Statistics in Society*, 183(4), 1387–1410. \doi{10.1111/rssa.12561}
#'
#' Källberg, D., & Waernbaum, I. (2023). Large Sample Properties of Entropy Balancing Estimators of Average Causal Effects. *Econometrics and Statistics*. \doi{10.1016/j.ecosta.2023.11.004}
#'
#' Wang, Y., & Zubizarreta, J. R. (2020). Minimal dispersion approximately balancing weights: Asymptotic properties and practical considerations. *Biometrika*, 107(1), 93–105. \doi{10.1093/biomet/asz050}
#'
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' @seealso
#' [optweightMV()] for estimating stable balancing weights for multivariate (i.e., multiple) treatments simultaneously.
#'
#' \CRANpkg{sbw}, which was the inspiration for this package and provides some additional functionality for binary treatments.
#'
#' \CRANpkg{WeightIt}, which provides a simplified interface to `optweight()` and a more efficient implementation of entropy balancing.
#'
#' @examplesIf rlang::is_installed("cobalt")
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
#' ow4 <- optweight(treat ~ age + educ + married + race +
#'                    nodegree + re74 + re75,
#'                  data = lalonde,
#'                  estimand = "ATE",
#'                  min.w = -Inf)
#'
#' summary(ow4)
#'
#' # Using `optweight.fit()`
#' treat <- lalonde$treat
#' covs <- splitfactor(lalonde[2:8], drop.first = "if2")
#'
#' ow.fit <- optweight.fit(covs,
#'                         treat = treat,
#'                         tols = .02,
#'                         estimand = "ATE",
#'                         norm = "l2")

#' @export
optweight <- function(formula, data = NULL, tols = 0,
                      estimand = "ATE", targets = NULL, target.tols = 0,
                      s.weights = NULL, b.weights = NULL, focal = NULL,
                      norm = "l2", min.w = 1e-8, verbose = FALSE, ...) {

  mcall <- match.call()

  #Process treat and covs from formula and data
  t.c <- get_covs_and_treat_from_formula2(formula, data, sep = "_")
  simple.covs <- t.c[["simple.covs"]]
  reported.covs <- t.c[["reported.covs"]]

  covs <- t.c[["model.covs"]]
  treat <- t.c[["treat"]]

  #Get treat type
  treat <- assign_treat_type(treat)
  treat.type <- .attr(treat, "treat.type")

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
  tols <- .process_tols_internal(covs, tols, reported.covs,
                                 tols_arg = "tols")

  #Process tols
  target.tols <- .process_tols_internal(covs, target.tols, reported.covs,
                                        tols_arg = "target.tols")

  #Process targets
  if (is_null(estimand) || is_not_null(targets)) {
    if (is_null(estimand) && is_null(targets)) {
      targets <- NA_real_
    }
    else if (is_not_null(estimand) && is_not_null(targets)) {
      .wrn("{.arg targets} are not {.val NULL}; ignoring {.arg estimand}")
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
                           target.tols = target.tols,
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
    .err("some weights are {.val NA}, which means something went wrong")
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
              norm = fit_out$norm,
              call = mcall,
              tols = tols,
              target.tols = target.tols,
              duals = duals,
              info = fit_out$info,
              solver = fit_out$solver)

  class(out) <- "optweight"

  out
}

#' @export
#' @rdname optweight
optweight.fit <- function(covs, treat, tols = 0,
                          estimand = "ATE", targets = NULL, target.tols = 0,
                          s.weights = NULL, b.weights = NULL, focal = NULL,
                          norm = "l2", std.binary = FALSE, std.cont = TRUE,
                          min.w = 1e-8, verbose = FALSE, solver = NULL, ...) {

  if ((!missing(covs) && is.list(covs) && !is.data.frame(covs)) || is_not_null(...get("covs.list")) ||
      (!missing(treat) && is.list(treat)) || is_not_null(...get("treat.list"))) {
    .err("{.fn optweight.fit} was called with list arguments; perhaps you meant to call {.fn optweightMV.fit}")
  }

  chk::chk_not_missing(covs, "`covs`")
  chk::chk_not_missing(treat, "`treat`")

  if (!is.numeric(covs) && (!is.data.frame(covs) || !all(apply(covs, 2L, is.numeric)))) {
    .err("all covariates must be numeric")
  }

  covs <- as.matrix(covs)

  treat.name <- .attr(treat, "treat.name") %or% "treat"

  treat.type <- {
    if (chk::vld_character_or_factor(treat) || is_binary(treat)) "cat"
    else "cont"
  }

  N <- length(treat)

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

  #Process tols and target.tols
  if (!inherits(tols, "optweight.tols") || is_null(.attr(tols, "internal.tols"))) {
    tols <- .process_tols_internal(covs, tols = tols, tols_found_in = "covs",
                                   tols_arg = "tols")
  }

  tols <- tols |>
    .attr("internal.tols") |>
    abs()

  if (!inherits(target.tols, "optweight.tols") || is_null(.attr(target.tols, "internal.tols"))) {
    target.tols <- .process_tols_internal(covs, tols = target.tols, tols_found_in = "covs",
                                          tols_arg = "target.tols")
  }

  target.tols <- target.tols |>
    .attr("internal.tols") |>
    abs()

  #Process estimand and targets
  if (treat.type == "cat") {
    treat <- qF(treat)

    if (is_null(estimand)) {
      if (is_null(targets)) {
        targets <- NA_real_
      }
    }
    else if (is_not_null(targets)) {
      .wrn("{.arg targets} are not {.val NULL}; ignoring {.arg estimand}")
      estimand <- focal <- NULL
    }
    else {
      chk::chk_string(estimand)
      estimand <- toupper(estimand)

      chk::chk_subset(estimand, c("ATE", "ATT", "ATC"))

      if (estimand %in% c("ATT", "ATC")) {
        if (is_null(focal)) {
          .err('{.arg focal} must be supplied when {.code estimand = "{estimand}"}')
        }

        focal <- as.character(focal)

        in_focal <- whichv(treat, focal)

        if (is_null(in_focal)) {
          .err("{.arg focal} must be the name of a level of treatment")
        }

        targets <- fmean(ss(covs, in_focal), w = sw[in_focal])

        if (!all(target.tols == 0) && !all(target.tols == Inf)) {
          .wrn('{.arg target.tols} is ignored when {.code estimand = "{estimand}"}')
        }

        target.tols[] <- tols
      }
      else {
        targets <- NULL # calculated automatically for ATE
        focal <- in_focal <- NULL
      }
    }

    if (!inherits(targets, "optweight.targets")) {
      targets <- .process_targets_internal(covs, targets = targets, sw = sw,
                                           targets_found_in = "covs")
    }
  }
  else {
    treat <- as.numeric(treat)

    if (is_null(estimand)) {
      if (is_null(targets)) {
        targets <- NA_real_
      }
    }
    else if (is_not_null(targets)) {
      .wrn("{.arg targets} are not {.val NULL}; ignoring {.arg estimand}")
      estimand <- NULL
    }
    else {
      chk::chk_string(estimand)
      estimand <- toupper(estimand)

      if (estimand != "ATE") {
        .err('{.arg estimand} cannot be "{estimand}" with continuous treatments')
      }

      targets <- NULL # calculated automatically for ATE
    }

    if (!allv(target.tols, 0)) {
      .wrn("{.arg target.tols} is ignored for continuous treatments. Setting all {.arg target.tols} to 0")
    }

    target.tols[] <- 0

    if (!inherits(targets, "optweight.targets")) {
      targets <- .process_targets_internal(covs, targets = targets, sw = sw,
                                           targets_found_in = "covs")
    }

    if (anyNA(targets)) {
      .err("all covariates must have a target when continuous treatments are used")
    }

    focal <- NULL
  }

  #Process norm
  norm <- process_norm(norm, sw, bw)

  #Process min.w
  min.w <- process_min.w(min.w, norm, bw)

  #Process solver
  solver <- process_solver(solver, norm, min.w)

  #Process args
  args <- make_process_opt_args(solver)(..., verbose = verbose)

  range_cons <- constraint_range_w(sw, min.w, focal, treat)

  constraint_df <- expand.grid(time = 0L,
                               type = "range_w",
                               constraint = list(range_cons),
                               stringsAsFactors = FALSE,
                               KEEP.OUT.ATTRS = FALSE) |>
    rbind(expand.grid(time = 1L,
                      type = c("mean_w", "balance", "target"),
                      constraint = list(NULL),
                      stringsAsFactors = FALSE,
                      KEEP.OUT.ATTRS = FALSE))

  bin.covs <- is_binary_col(covs)

  targeted <- !is.na(targets)

  if (treat.type == "cat") {
    n <- fsum(sw * bw, g = treat)

    sds <- {
      if (is_not_null(focal))
        col.w.sd(ss(covs, in_focal), w = sw[in_focal], bin.vars = bin.covs)
      else
        sqrt(colMeans(col.w.v(covs, g = treat, w = sw, bin.vars = bin.covs)))
    }

    vars.to.standardize <- rep_with(FALSE, tols)
    if (std.binary) vars.to.standardize[bin.covs] <- TRUE
    if (std.cont) vars.to.standardize[!bin.covs] <- TRUE

    to_std <- which(vars.to.standardize & !check_if_zero(sds))

    if (is_not_null(to_std)) {
      covs[, to_std] <- ss(covs, j = to_std) %r/% sds[to_std]
      targets[to_std] <- targets[to_std] / sds[to_std]
    }

    constraint_df[["constraint"]][constraint_df[["time"]] == 1L] <- list(
      mean_w = constraint_mean_w_cat(treat, sw, n),
      balance = constraint_balance_cat(covs, treat, sw, tols, n),
      target = constraint_target_cat(covs, treat, sw, targets,
                                     target.tols, targeted, n, focal)
    )
  }
  else {
    n <- sum(sw * bw)

    sds <- col.w.sd(covs, w = sw, bin.vars = bin.covs)

    covs <- covs %r-% targets #center covs at targets (which will be eventual means)

    treat <- fscale(treat, w = sw)

    to_std <- !check_if_zero(sds)

    if (is_not_null(to_std)) {
      covs[, to_std] <- ss(covs, j = to_std) %r/% sds[to_std]
      targets[to_std] <- targets[to_std] / sds[to_std]
    }

    constraint_df[["constraint"]][constraint_df[["time"]] == 1L] <- list(
      mean_w = constraint_mean_w_cont(sw, n),
      balance = constraint_balance_cont(covs, treat, sw, tols, n),
      target = constraint_target_cont(covs, treat, sw, n, treat.name)
    )
  }

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

  class(out) <- "optweight.fit"

  out
}

#' @exportS3Method print optweight
print.optweight <- function(x, ...) {
  treat.type <- .attr(x[["treat"]], "treat.type")

  if (is_not_null(treat.type)) {
    treat.type[treat.type == "multinomial"] <- "multi-category"
  }

  cat(sprintf("An %s object\n", .it(class(x)[1L])))

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
                               fnunique(x[["treat"]]),
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
