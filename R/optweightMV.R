#' Stable Balancing Weights for Multivariate Treatments
#'
#' Estimates stable balancing weights for the supplied multivariate (i.e., multiple) treatments and covariates. The degree of balance for each covariate is specified by `tols.list`. See Zubizarreta (2015) and Wang & Zubizarreta (2020) for details of the properties of the weights and the methods used to fit them.
#'
#' @inheritParams optweight
#' @param formula.list a list of formulas, each with a treatment variable on the left hand side and the covariates to be balanced on the right hand side.
#' @param data an optional data set in the form of a data frame that contains the variables in `formula.list`.
#' @param tols.list a list of vectors of balance tolerance values for each covariate for each treatment. The resulting weighted balance statistics will be at least as small as these values. If only one value is supplied, it will be applied to all covariates. See Details. Default is 0 for all covariates.
#' @param estimand the desired estimand, which determines the target population. Only "ATE" or `NULL` are supported. `estimand` is ignored when `targets` is non-`NULL`. If both `estimand` and `targets` are `NULL`, no targeting will take place.
#' @param targets an optional vector of target population mean values for each covariate. The resulting weights will yield sample means within `tols`/2 units of the target values for each covariate. If `NULL` or all `NA`, `estimand` will be used to determine targets. Otherwise, `estimand` is ignored. If any target values are `NA`, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance; this is only allowed if all treatments are binary or multi-category. Can also be the output of a call to [process_targets()]. See Details.
#' @param s.weights a vector of sampling weights. For `optweightMV()`, can also be the name of a variable in `data` that contains sampling weights.
#' @param b.weights a vector of base weights. If supplied, the desired norm of the distance between the estimated weights and the base weights is minimized. For `optweightMV()`, can also the name of a variable in `data` that contains base weights.
#' @param covs.list a list containing one numeric matrix of covariates to be balanced for each treatment.
#' @param treat.list a list containing one vector of treatment statuses for each treatment.
#' @param \dots for `optweightMV()`, additional arguments passed to `optweightMV.fit()`, including options that are passed to the settings function corresponding to `solver`.
#'
#' @returns
#' For `optweightMV()`, an `optweightMV` object with the following elements:
#' \item{weights}{The estimated weights, one for each unit.}
#' \item{treat.list}{A list of the values of the treatment variables.}
#' \item{covs.list}{A list of the covariates for each treatment used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
#' \item{s.weights}{The provided sampling weights.}
#' \item{b.weights}{The provided base weights.}
#' \item{call}{The function call.}
#' \item{tols}{A list of tolerance values for each covariate for each treatment.}
#' \item{duals}{A list of data.frames containing the dual variables for each covariate for each treatment. See [optweight()] for interpretation of these values.}
#' \item{info}{A list containing information about the performance of the optimization at termination.}
#' \item{norm}{The `norm` used.}
#' \item{solver}{The `solver` used.}
#'
#' For `optweightMV.fit()`, an `optweightMV.fit` object with the following elements:
#' \item{w}{The estimated weights, one for each unit.}
#' \item{duals}{A data.frame containing the dual variables for each covariate.}
#' \item{info}{A list containing information about the performance of the optimization at termination.}
#' \item{norm}{The `norm` used.}
#' \item{solver}{The `solver` used.}
#'
#' @details
#' `optweightMV()` is the primary user-facing function for estimating stable balancing weights for multivariate treatments. The optimization is performed by the lower-level function `optweightMV.fit()`, which transforms the inputs into the required inputs for the optimization functions and then supplies the outputs (the weights, dual variables, and convergence information) back to `optweightMV()`. Little processing of inputs is performed by `optweightMV.fit()`, as this is normally handled by `optweightMV()`.
#'
#' See [optweight()] for more information about balance tolerances (i.e., those specified in `tols.list`), `targets`, `norm`, `solver`, and convergence failure.
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
#' [optweight()] for more information on the optimization, specifications, and options.
#'
#' @examplesIf rlang::is_installed("cobalt")
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' # Balancing two treatments
#' (ow1 <- optweightMV(list(treat ~ age + educ + race + re74,
#'                          re75 ~ age + educ + race + re74),
#'                     data = lalonde))
#'
#' summary(ow1)
#'
#' bal.tab(ow1)

#' @export
optweightMV <- function(formula.list, data = NULL, tols.list = list(0), estimand = "ATE",
                        targets = NULL, s.weights = NULL, b.weights = NULL,
                        norm = "l2", min.w = 1e-8, verbose = FALSE, ...) {

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
  duals <- lapply(times, function(i) process_duals(fit_out$duals[[i]], tols.list[[i]]))

  out <- list(weights = fit_out$w,
              treat.list = treat.list,
              covs.list = simple.covs.list,
              s.weights = sw,
              b.weights = bw,
              norm = fit_out$norm,
              call = mcall,
              tols = tols.list,
              duals = duals,
              info = fit_out$info,
              solver = fit_out$solver)

  class(out) <- c("optweightMV", "optweight")

  out
}

#' @export
#' @rdname optweightMV
optweightMV.fit <- function(covs.list, treat.list, tols.list = list(0), estimand = "ATE",
                            targets = NULL, s.weights = NULL, b.weights = NULL,
                            norm = "l2", std.binary = FALSE, std.cont = TRUE, min.w = 1e-8, verbose = FALSE,
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

  #Process norm
  norm <- process_norm(norm, sw, bw)

  #Process min.w
  min.w <- process_min.w(min.w, norm, bw)

  #Process solver
  solver <- process_solver(solver, norm, min.w)

  #Process args
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
              out = opt_out$out,
              norm = norm,
              solver = solver)

  class(out) <- "optweightMV.fit"

  out
}

#' @exportS3Method print optweightMV
print.optweightMV <- function(x, ...) {
  treat.types <- vapply(x[["treat.list"]], attr, character(1L), "treat.type")
  treat.types[treat.types == "multinomial"] <- "multi-category"

  treat.names <- vapply(x[["treat.list"]], attr, character(1L), "treat.name")

  cat(sprintf("An %s object\n", .it(class(x)[1L])))

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
