#' Fitting Function for Stable Balancing Weights
#'
#' `optweight.fit()` and `optweightMV.fit()` perform the optimization for [optweight()] and [optweightMV()] and should, in most cases, not be used directly. Little processing of inputs is performed, so they must be given exactly as described below.
#'
#' @param covs a numeric matrix of covariates to be balanced.
#' @param treat a vector of treatment statuses. Non-numeric (i.e., factor or character) vectors are allowed.
#' @param covs.list a list containing one numeric matrix of covariates to be balanced for each treatment.
#' @param treat.list a list containing one vector of treatment statuses for each treatment.
#' @param tols a vector of balance tolerance values for each covariate. Default is 0.
#' @param tols.list a list of balance tolerance vectors, one for each treatment, each with a value for each covariate.
#' @param estimand the desired estimand, which determines the target population. For binary treatments, can be "ATE", "ATT", "ATC", or `NULL`. For multi-category treatments, can be "ATE", "ATT", or `NULL`. For continuous treatments, can be "ATE" or `NULL`. The default for both is "ATE". For `optweightMV.fit()`, only "ATE" or `NULL` are supported. `estimand` is ignored when `targets` is non-`NULL`. If both `estimand` and `targets` are `NULL`, no targeting will take place.
#' @param targets an optional vector of target population mean values for each baseline covariate. The resulting weights will yield sample means within `tols`/2 units of the target values for each covariate. If `NULL` or all `NA`, `estimand` will be used to determine targets. Otherwise, `estimand` is ignored. If any target values are `NA`, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance.
#' @param s.weights an optional vector of sampling weights. Default is a vector of 1s.
#' @param b.weights an optional vector of base weights. Default is a vector of 1s.
#' @param focal when multi-categorical treatments are used and the `estimand = "ATT"`, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to resemble the focal group.
#' @param norm `character`; a string containing the name of the norm corresponding to the objective function to minimize. Allowable options include `"l1"` for the L1 norm, `"l2"` for the L2 norm (the default), `"linf"` for the L\eqn{\infty} norm, `"entropy"` for the negative entropy, and `"log"` for the sum of the negative logs. See Details.
#' @param std.binary,std.cont `logical`; whether the tolerances are in standardized mean units (`TRUE`) or raw units (`FALSE`) for binary variables and continuous variables, respectively. The default is `FALSE` for `std.binary` because raw proportion differences make more sense than standardized mean difference for binary variables. These arguments are analogous to the `binary` and `continuous` arguments in `bal.tab()` in \pkg{cobalt}.
#' @param min.w `numeric`; a single value less than 1 for the smallest allowable weight. Some analyses require nonzero weights for all units, so a small, nonzero minimum may be desirable. The default is `1e-8` (\eqn{10^{-8}}), which does not materially change the properties of the weights from a minimum of 0 but prevents warnings in some packages that use weights to estimate treatment effects. When `norm` is `"entropy"` or `"log"` and `min.w <= 0`, `min.w` will be set to the smallest nonzero value.
#' @param verbose `logical`; whether information on the optimization problem solution should be printed. Default is `FALSE`.
#' @param solver string; the name of the optimization solver to use. Allowable options depend on `norm`. Default is to use whichever eligible solver is installed, if any, or the default solver for the corresponding `norm`. See Details for information.
#' @param \dots Options that are passed to the settings function corresponding to `solver`.
#'
#' @returns
#' An `optweight.fit` or `optweightMV.fit` object with the following elements:
#'   \item{w}{The estimated weights, one for each unit.}
#'   \item{duals}{A data.frame containing the dual variables for each covariate (for `optweight.fit()`), or a list thereof (for `optweightMV.fit()`). See `vignette("optweight")` for interpretation of these values.}
#'   \item{info}{A list containing information about the performance of the optimization at termination.}
#'
#' @details
#' `optweight.fit()` and `optweightMV.fit()` transform the inputs into the required inputs for the optimization functions, which are (sparse) matrices and vectors, and then supplies the outputs (the weights, dual variables, and convergence information) back to [optweight()] or [optweightMV()]. Little processing of inputs is performed, as this is normally handled by `optweight()` or `optweightMV()`.
#'
#' Target and balance constraints are applied to the product of the estimated weights and the sampling weights. In addition,the  sum of the product of the estimated weights and the sampling weights is constrained to be equal to the sum of the product of the base weights and sampling weights. For binary and multi-category treatments, these constraints apply within each treatment group.
#'
#' ## `norm`
#'
#' The objective function for the optimization problem is \eqn{f\left(w_i, b_i, s_i\right)}, where \eqn{w_i} is the estimated weight for unit \eqn{i}, \eqn{s_i} is the sampling weight for unit \eqn{i} (supplied by `s.weights`) and \eqn{b_i} is the base weight for unit \eqn{i} (supplied by `b.weights`). The `norm` argument determines \eqn{f(.,.,.)}, as detailed below:
#'
#' * when `norm = "l2"`, \eqn{f\left(w_i, b_i, s_i\right) = \frac{1}{n} \sum_i {s_i(w_i - b_i)^2}}
#' * when `norm = "l1"`, \eqn{f\left(w_i, b_i, s_i\right) = \frac{1}{n} \sum_i {s_i \vert w_i - b_i \vert}}
#' * when `norm = "linf"`, \eqn{f\left(w_i, b_i, s_i\right) = \max_i {\vert w_i - b_i \vert}}
#' * when `norm = "entropy"`, \eqn{f\left(w_i, b_i, s_i\right) = \frac{1}{n} \sum_i {s_i w_i \log \frac{w_i}{b_i}}}
#' * when `norm = "log"`, \eqn{f\left(w_i, b_i, s_i\right) = \frac{1}{n} \sum_i {-s_i \log \frac{w_i}{b_i}}}
#'
#' By default, `s.weights` and `b.weights` are set to 1 for all units unless supplied. `b.weights` must be positive when `norm` is `"entropy"` or `"log"`, and `norm = "linf"` cannot be used when `s.weights` are supplied.
#'
#' When `norm = "l2"` and both `s.weights` and `b.weights` are `NULL`, weights are estimated to maximize the effective sample size. When `norm = "entropy"`, the estimated weights are equivalent to entropy balancing weights (Källberg & Waernbaum, 2023). When `norm = "log"`, `b.weights` are ignored in the optimization, as they do not affect the estimated weights.
#'
#' ## `solver`
#'
#' The `solver` argument controls which optimization solver is used. Different solvers are compatible with each `norm`. See the table below for allowable options, which package they require, which function does the solving, and which function controls the settings.
#'
#' | `solver` | `norm` | Package | Solver function | Settings function |
#' |----------|--------|---------|-----------------|-------------------|
#' | `"osqp"` | `"l2"`, `"l1"`, `"linf"` | \CRANpkg{osqp} | [osqp::solve_osqp()] | [osqp::osqpSettings()] |
#' | `"highs"` | `"l2"`, `"l1"`, `"linf"` | \CRANpkg{highs} | \pkgfun{highs}{highs_solve} | \pkgfun{highs}{highs_control} / \pkgfun{highs}{highs_available_solver_options} |
#' | `"lpsolve"` | `"l1"`, `"linf"` | \CRANpkg{lpSolve} | \pkgfun{lpSolve}{lp} | . |
#' | `"scs"` | `"entropy"`, `"log"` | \CRANpkg{scs} | \pkgfun{scs}{scs} | \pkgfun{scs}{scs_control} |
#' | `"clarabel"` | `"entropy"`, `"log"` | \CRANpkg{clarabel} | \pkgfun{clarabel}{clarabel} | \pkgfun{clarabel}{clarabel_control} |
#'
#' Note that `"lpsolve"` can only be used when `min.w` is nonnegative.
#'
#' The default `solver` for each `norm` is as follows:
#'
#' | `norm` | Default `solver` |
#' |--------|------------------|
#' | `"l2"` | `"osqp"` |
#' | `"l1"` | `"highs"` |
#' | `"linf"` | `"highs"` |
#' | `"entropy"` | `"scs"` |
#' | `"log"` | `"scs"` |
#'
#' If the package corresponding to a default `solver` is not installed but the package for a different eligible solver is, that will be used. Otherwise, you will be asked to install the required package. \pkg{osqp} is required for \pkg{optweight}, and so will be the default for the `"l1"` and `"linf"` norms if \pkg{highs} is not installed. The default package is the one has shown good performance for the given norm; generally, all eligible solvers perform about equally well in terms of accuracy but differ in time taken.
#'
#' ## Solving Convergence Failure
#'
#' Sometimes the optimization will fail to converge at a solution. There are a variety of reasons why this might happen, which include that the constraints are nearly impossible to satisfy or that the optimization surface is relatively flat. It can be hard to know the exact cause or how to solve it, but this section offers some solutions one might try. Typically, solutions can be found most easily when using the `"l2"` norm; other norms, especially `"linf"` and `"l1"`, are more likely to see problems.
#'
#' Rarely is the problem too few iterations, though this is possible. Most problems can be solved in the default 200,000 iterations, but sometimes it can help to increase this number with the `max_iter` argument. Usually, though, this just ends up taking more time without a solution found.
#'
#' If the problem is that the constraints are too tight, it can be helpful to loosen the constraints. Sometimes examining the dual variables of a solution that has failed to converge can reveal which constraints are causing the problem.
#'
#' Sometimes a suboptimal solution is possible; such a solution does not satisfy the constraints exactly but will come pretty close. To allow these solutions, the argument `eps` can be increased to larger values.
#'
#' Sometimes using a different solver can improve performance. Using the default `solver` for each `norm`, as described above, can reduce the probability of convergence failures.
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
#' [optweight()] and [optweightMV()] which you should use for estimating the balancing weights, unless you know better.
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' treat <- lalonde$treat
#' covs <- splitfactor(lalonde[2:8], drop.first = "if2")
#'
#' ow.fit <- optweight.fit(covs,
#'                         treat,
#'                         tols = .02,
#'                         estimand = "ATE",
#'                         norm = "l2")

#' @export
optweight.fit <- function(covs, treat, tols = 0, estimand = "ATE", targets = NULL,
                          s.weights = NULL, b.weights = NULL, focal = NULL, norm = "l2",
                          std.binary = FALSE, std.cont = TRUE, min.w = 1e-8, verbose = FALSE,
                          solver = NULL, ...) {

  if ((!missing(covs) && is.list(covs) && !is.data.frame(covs)) || is_not_null(...get("covs.list")) ||
      (!missing(treat) && is.list(treat)) || is_not_null(...get("treat.list"))) {
    .err("`optweight.fit()` was called with list arguments; perhaps you meant to call `optweightMV.fit()")
  }

  chk::chk_not_missing(covs, "`covs`")
  chk::chk_not_missing(treat, "`treat`")

  if (!is.numeric(covs) && (!is.data.frame(covs) || !all(apply(covs, 2L, is.numeric)))) {
    .err("all covariates must be numeric")
  }

  covs <- as.matrix(covs)

  treat.name <- if_null_then(attr(treat, "treat.name"), "treat")

  treat.type <- {
    if (chk::vld_character_or_factor(treat) || is_binary(treat)) "cat"
    else "cont"
  }

  N <- length(treat)

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

  #Process tols
  if (!inherits(tols, "optweight.tols") || is_null(attr(tols, "internal.tols"))) {
    tols <- .process_tols_internal(covs, tols = tols, tols_found_in = "covs")
  }

  tols <- tols |>
    attr("internal.tols") |>
    abs()

  #Process estimand and targets
  if (treat.type == "cat") {
    treat <- as.character(treat)

    if (is_null(estimand)) {
      if (is_null(targets)) {
        targets <- NA_real_
      }
    }
    else if (is_not_null(targets)) {
      .wrn("`targets` are not `NULL`; ignoring `estimand`")
      estimand <- focal <- NULL
    }
    else {
      chk::chk_string(estimand)
      estimand <- toupper(estimand)

      chk::chk_subset(estimand, c("ATE", "ATT", "ATC"))

      if (estimand %in% c("ATT", "ATC")) {
        if (is_null(focal)) {
          .err(sprintf("`focal` must be supplied when `estimand = %s`",
                       add_quotes(estimand)))
        }

        focal <- as.character(focal)

        in_focal <- which(treat == focal)

        if (is_null(in_focal)) {
          .err("`focal` must be the name of a level of treatment")
        }

        targets <- col.w.m(covs[in_focal, , drop = FALSE],
                           w = sw[in_focal])
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
      .wrn("`targets` are not `NULL`; ignoring `estimand`")
      estimand <- NULL
    }
    else {
      chk::chk_string(estimand)
      estimand <- toupper(estimand)

      if (estimand != "ATE") {
        .err(sprintf("`estimand` cannot be %s with continuous treatments",
                     add_quotes(estimand)))
      }

      targets <- NULL # calculated automatically for ATE
    }

    if (!inherits(targets, "optweight.targets")) {
      targets <- .process_targets_internal(covs, targets = targets, sw = sw,
                                           targets_found_in = "covs")
    }

    if (anyNA(targets)) {
      .err("all covariates must have a target when continuous treatments are used")
    }

    focal <- NULL
  }

  chk::chk_string(norm)
  norm <- tolower(norm)
  chk::chk_subset(norm, allowable_norms())

  if (norm == "linf" && !all_the_same(sw)) {
    .err("the L-inf norm cannot be used with sampling weights")
  }

  chk::chk_number(min.w)
  chk::chk_lte(min.w, mean(bw))

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
                               type = c("range_w", "mean_w", "balance", "target"),
                               constraint = list(NULL),
                               stringsAsFactors = FALSE,
                               KEEP.OUT.ATTRS = FALSE)

  bin.covs <- is_binary_col(covs)

  range_cons <- constraint_range_w(sw, min.w, focal, treat)

  if (treat.type == "cat") {
    unique.treats <- sort(unique(treat))

    n <- vapply(unique.treats,
                function(t) sum(sw[treat == t] * bw[treat == t]),
                numeric(1L))

    sds <- {
      if (is_not_null(focal))
        sqrt(col.w.v(covs[in_focal, , drop = FALSE],
                     w = sw[in_focal],
                     bin.vars = bin.covs))
      else
        sqrt(colMeans(do.call("rbind", lapply(unique.treats, function(t) {
          in_treat <- which(treat == t)

          col.w.v(covs[in_treat, , drop = FALSE],
                  w = sw[in_treat], bin.vars = bin.covs)
        }))))
    }

    targeted <- !is.na(targets)

    balanced <- !targeted

    treat.sd <- NA_real_
    treat.mean <- NA_real_

    vars.to.standardize <- rep_with(FALSE, tols)
    if (std.binary) vars.to.standardize[bin.covs] <- TRUE
    if (std.cont) vars.to.standardize[!bin.covs] <- TRUE

    to_std <- which(vars.to.standardize & !check_if_zero(sds))

    if (is_not_null(to_std)) {
      covs[, to_std] <- mat_div(covs[, to_std, drop = FALSE], sds[to_std])
      targets[to_std] <- targets[to_std] / sds[to_std]
    }

    constraint_df[["constraint"]] <- list(
      range_w = range_cons,
      mean_w = constraint_mean_w_cat(treat, unique.treats, sw, n),
      balance = constraint_balance_cat(covs, treat, sw, tols,
                                       balanced, unique.treats, n),
      target = constraint_target_cat(covs, treat, sw, targets,
                                     tols, targeted, unique.treats, n, focal)
    )
  }
  else {
    n <- sum(sw * bw)

    sds <- sqrt(col.w.v(covs, w = sw, bin.vars = bin.covs))

    targeted <- !is.na(targets)

    balanced <- rep_with(TRUE, targeted)

    covs <- center(covs, at = targets) #center covs at targets (which will be eventual means)

    treat.sd <- sqrt(col.w.v(treat, w = sw))
    treat.mean <- w.m(treat, w = sw)

    treat <- (treat - treat.mean) / treat.sd

    to_std <- !check_if_zero(sds)

    if (is_not_null(to_std)) {
      covs[, to_std] <- mat_div(covs[, to_std, drop = FALSE], sds[to_std])
      targets[to_std] <- targets[to_std] / sds[to_std]
    }

    constraint_df[["constraint"]] <- list(
      range_w = range_cons,
      mean_w = constraint_mean_w_cont(sw, n),
      balance = constraint_balance_cont(covs, treat, sw, tols, balanced, n),
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
              duals = duals[[1L]],
              info = opt_out$info_out,
              out = opt_out$out)

  class(out) <- "optweight.fit"

  out
}

#' @export
#' @rdname optweight.fit
optweightMV.fit <- function(covs.list, treat.list, tols.list = list(0), estimand = "ATE", targets = NULL,
                            s.weights = NULL, b.weights = NULL, norm = "l2",
                            std.binary = FALSE, std.cont = TRUE, min.w = 1e-8, verbose = FALSE,
                            solver = NULL, ...) {

  chk::chk_not_missing(treat.list, "`treat.list`")
  chk::chk_not_missing(covs.list, "`covs.list`")

  chk::chk_list(treat.list)
  chk::chk_list(covs.list)

  times <- seq_along(covs.list)

  treat.types <- treat.names <- character(length(times))

  for (i in times) {
    if (!is.numeric(covs.list[[i]]) && (!is.data.frame(covs.list[[i]]) || !all(apply(covs.list[[i]], 2L, is.numeric)))) {
      .err("all covariates must be numeric")
    }

    covs.list[[i]] <- as.matrix(covs.list[[i]])

    treat.types[i] <- {
      if (chk::vld_character_or_factor(treat.list[[i]]) || is_binary(treat.list[[i]])) "cat"
      else "cont"
    }

    treat.names[i] <- if_null_then(attr(treat.list[[i]], "treat.name"),
                                   names(treat.list)[i],
                                   sprintf("treatment %s", i))
  }

  N <- length(treat.list[[1L]])

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

  #Process tols
  chk::chk_not_missing(tols.list, "`tols.list`")
  chk::chk_list(tols.list)

  if (length(tols.list) == 1L) {
    tols.list <- tols.list[rep_with(1, times)]
  }

  for (i in times) {
    if (!inherits(tols.list[[i]], "optweight.tols") || is_null(attr(tols.list[[i]], "internal.tols"))) {
      tols.list[[i]] <- .process_tols_internal(covs.list[[i]], tols.list[[i]], tols_found_in = "covs.list")
    }

    tols.list[[i]] <- tols.list[[i]] |>
      attr("internal.tols") |>
      abs()
  }

  #Process estimand and targets
  if (is_null(estimand)) {
    if (is_null(targets)) {
      targets <- NA_real_
    }
  }
  else if (is_not_null(targets)) {
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

    targets <- NULL # calculated automatically for ATE
  }

  if (!inherits(targets, "optweight.targets")) {
    targets <- .process_targets_internal(cbind_distinct(covs.list), targets = targets, sw = sw,
                                         targets_found_in = "covs.list")
  }

  if (any_apply(times, function(i) {
    treat.types[i] == "cont" && anyNA(targets[colnames(covs.list[[i]])])
  })) {
    .err("all covariates associated with a continuous treatment must have a target")
  }

  chk::chk_string(norm)
  norm <- tolower(norm)
  chk::chk_subset(norm, allowable_norms())

  if (norm == "linf" && !all_the_same(sw)) {
    .err("the L-inf norm cannot be used with sampling weights")
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

  constraint_df <- expand.grid(time = times,
                               type = c("range_w", "mean_w", "balance", "target"),
                               constraint = list(NULL),
                               stringsAsFactors = FALSE,
                               KEEP.OUT.ATTRS = FALSE)

  range_cons <- constraint_range_w(sw, min.w)

  unique.treats <- n <- bin.covs.list <- sds <- targets.list <- targeted <-
    balanced <- treat.sds <- treat.means <- make_list(length(times))

  for (i in times) {
    bin.covs.list[[i]] <- is_binary_col(covs.list[[i]])

    targets.list[[i]] <- targets[colnames(covs.list[[i]])]

    targeted[[i]] <- !is.na(targets.list[[i]])

    if (treat.types[i] == "cat") {
      treat.list[[i]] <- as.character(treat.list[[i]])

      unique.treats[[i]] <- sort(unique(treat.list[[i]]))

      n[[i]] <- vapply(unique.treats[[i]],
                       function(t) sum(sw[treat.list[[i]] == t] * bw[treat.list[[i]] == t]),
                       numeric(1L))

      sds[[i]] <- sqrt(colMeans(do.call("rbind", lapply(unique.treats[[i]], function(t) {
        in_treat <- which(treat.list[[i]] == t)

        col.w.v(covs.list[[i]][in_treat, , drop = FALSE],
                w = sw[in_treat], bin.vars = bin.covs.list[[i]])
      }))))

      balanced[[i]] <- !targeted[[i]]

      treat.sds[[i]] <- NA_real_
      treat.means[[i]] <- NA_real_

      vars.to.standardize <- rep_with(FALSE, tols.list[[i]])
      if (std.binary) vars.to.standardize[bin.covs.list[[i]]] <- TRUE
      if (std.cont) vars.to.standardize[!bin.covs.list[[i]]] <- TRUE

      to_std <- which(vars.to.standardize & !check_if_zero(sds[[i]]))

      if (is_not_null(to_std)) {
        covs.list[[i]][, to_std] <- mat_div(covs.list[[i]][, to_std, drop = FALSE], sds[[i]][to_std])
        targets.list[[i]][to_std] <- targets.list[[i]][to_std] / sds[[i]][to_std]
      }

      constraint_df[["constraint"]][constraint_df[["time"]] == i] <- list(
        range_w = if (i == 1L) range_cons,
        mean_w = constraint_mean_w_cat(treat.list[[i]], unique.treats[[i]], sw, n[[i]]),
        balance = constraint_balance_cat(covs.list[[i]], treat.list[[i]], sw, tols.list[[i]],
                                         balanced[[i]], unique.treats[[i]], n[[i]]),
        target = constraint_target_cat(covs.list[[i]], treat.list[[i]], sw, targets.list[[i]],
                                       tols.list[[i]], targeted[[i]], unique.treats[[i]], n[[i]])
      )
    }
    else {
      treat.list[[i]] <- as.numeric(treat.list[[i]])

      n[[i]] <- sum(sw * bw)

      sds[[i]] <- sqrt(col.w.v(covs.list[[i]], w = sw,
                               bin.vars = bin.covs.list[[i]]))

      balanced[[i]] <- rep_with(TRUE, targeted[[i]])

      covs.list[[i]][] <- center(covs.list[[i]], at = targets.list[[i]]) #center covs at targets (which will be eventual means)

      treat.sds[[i]] <- sqrt(col.w.v(treat.list[[i]], w = sw))
      treat.means[[i]] <- w.m(treat.list[[i]], w = sw)

      treat.list[[i]][] <- (treat.list[[i]] - treat.means[[i]]) / treat.sds[[i]]

      to_std <- which(!check_if_zero(sds[[i]]))

      if (is_not_null(to_std)) {
        covs.list[[i]][, to_std] <- mat_div(covs.list[[i]][, to_std, drop = FALSE], sds[[i]][to_std])
        targets.list[[i]][to_std] <- targets.list[[i]][to_std] / sds[[i]][to_std]
      }

      constraint_df[["constraint"]][constraint_df[["time"]] == i] <- list(
        range_w = if (i == 1L) range_cons,
        mean_w = constraint_mean_w_cont(sw, n[[i]]),
        balance = constraint_balance_cont(covs.list[[i]], treat.list[[i]], sw, tols.list[[i]],
                                          balanced[[i]], n[[i]]),
        target = constraint_target_cont(covs.list[[i]], treat.list[[i]], sw, n[[i]],
                                        treat.names[i])
      )
    }
  }

  constraint_df <- constraint_df |>
    prep_constraint_df(norm, bw, sw) |>
    prep_constraint_df_for_solver(solver)

  objective <- prep_objective(norm, bw, sw)

  opt_out <- opt_fit(constraint_df, objective, args, N,
                     solver = solver)

  w <- extract_weights(opt_out, N, min.w, range_cons)

  duals <- extract_duals(constraint_df, opt_out$dual_out, times)

  out <- list(w = w,
              duals = duals,
              info = opt_out$info_out,
              out = opt_out$out)

  class(out) <- "optweightMV.fit"

  out
}
