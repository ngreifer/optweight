#' Construct and Check Targets Input
#'
#' Checks whether proposed target population means values for `targets` are suitable in number and order for submission to [optweight()], [optweightMV()], and [optweight.svy()], and returns an object that can supplied to the `targets` argument of these functions.
#'
#' @inheritParams optweight.svy
#' @inheritParams process_tols
#' @param targets a vector of target population mean values for each covariate. These should be in the order corresponding to the order of the corresponding variable in `formula`, except for interactions, which will appear after all lower-order terms. For factor variables, a target value must be specified for each level of the factor, and these values must add up to 1. If `NULL`, the current sample means will be produced (weighted by `s.weights`). If `NA`, an `NA` vector named with the covariate names will be produced.
#' @param x an `optweight.targets` object; the output of a call to `process_targets()`.
#'
#' @returns
#' An `optweight.targets` object, which is a named vector of target population mean values, one for each (expanded) covariate specified in `formula`. This should be used as an input to the `targets` argument of [optweight()], [optweightMV()], and [optweight.svy()].
#'
#' @details
#' The purpose of `process_targets()` is to allow users to ensure that their proposed input to `targets` in [optweight()], [optweightMV()], and [optweight.svy()] is correct both in the number of entries and their order. This is especially important when factor variables and interactions are included in the formula because factor variables are split into several dummies and interactions are moved to the end of the variable list, both of which can cause some confusion and potential error when entering `targets` values.
#'
#' Factor variables are internally split into a dummy variable for each level, so the user must specify a target population mean value for each level of the factor. These must add up to 1, and an error will be displayed if they do not. These values represent the proportion of units in the target population with each factor level.
#'
#' Interactions (e.g., `a:b` or `a*b` in the `formula` input) are always sent to the end of the variable list even if they are specified elsewhere in the `formula`. It is important to run `process_targets()` to ensure the order of the proposed `targets` corresponds to the represented order of covariates used in the formula. You can run `process_targets(., targets = NA)` to see the order of covariates that is required without specifying any targets.
#'
#' @seealso [process_tols()]
#'
#' @examplesIf rlang::is_installed("cobalt")
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' # Generating targets; means by default
#' targets <- process_targets(~ age + race + married +
#'                              nodegree + re74,
#'                            data = lalonde)
#'
#' # Notice race is split into three values
#' targets
#'
#' # Generating targets; NA by default
#' targets <- process_targets(~ age + race + married +
#'                              nodegree + re74,
#'                            data = lalonde,
#'                            targets = NA)
#' targets
#'
#' # Can also supply just a dataset
#' covs <- lalonde |>
#'   subset(select = c(age, race, married,
#'                     nodegree, re74))
#'
#' targets <- process_targets(covs)
#'
#' targets

#' @export
process_targets <- function(formula, data = NULL, targets = NULL, s.weights = NULL) {

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
  t.c <- formula |>
    terms(data = data) |>
    delete.response() |>
    get_covs_and_treat_from_formula2(data, terms = TRUE, sep = "_")

  reported.covs <- t.c[["reported.covs"]]
  covs <- t.c[["model.covs"]]

  if (is_null(covs)) {
    .err("no covariates were specified")
  }

  check_missing_covs(reported.covs)

  #Process s.weights
  sw <- process_s.weights(s.weights, data)

  .process_targets_internal(covs, targets, sw, reported.covs,
                            if (formula.present) "formula" else "data")
}

.process_targets_internal <- function(model.covs, targets = NULL, sw = NULL, formula.covs = NULL,
                                      targets_found_in = "formula") {

  model.vars <- colnames(model.covs)

  if (is_null(targets)) {
    internal.targets <- col.w.m(model.covs, sw)
    ATE <- TRUE
  }
  else if (all(is.na(targets)) && is_null(names(targets))) {
    internal.targets <- setNames(rep_with(NA_real_, model.vars),
                                 model.vars)
    ATE <- FALSE
  }
  else {
    chk::chk_numeric(targets)

    if (is_null(names(targets)) || !all(nzchar(names(targets)))) {
      if (length(targets) != length(model.vars)) {
        .err(sprintf("`targets` must contain %s values, but %s were included",
                     length(model.vars), length(targets)))
      }

      names(targets) <- model.vars
    }

    if (!all(names(targets) %in% model.vars)) {
      .err(sprintf("all variables named in `targets` must be present in `%s`",
                   targets_found_in))
    }

    model.covs.means <- col.w.m(model.covs, sw)

    overlap <- intersect(names(targets), model.vars)

    internal.targets <- model.covs.means
    internal.targets[overlap] <- targets[overlap]

    ATE <- !anyNA(internal.targets) && all(check_if_zero(internal.targets - model.covs.means))
  }

  if (is_not_null(formula.covs)) {
    formula.vars <- attr(attr(formula.covs, "terms"), "term.labels")

    if (is_null(formula.vars)) {
      formula.vars <- "(Intercept)"
      attr(model.covs, "assign") <- 1
    }

    original.variables <- formula.vars[attr(model.covs, "assign")] |> setNames(model.vars)

    for (v in formula.vars) {
      if (attr(terms(formula.covs), "order")[formula.vars == v] == 1 &&
          attr(terms(formula.covs), "dataClasses")[formula.vars == v] == "factor" &&
          !anyNA(internal.targets[original.variables == v]) &&
          !check_if_zero(sum(internal.targets[original.variables == v]) - 1)) {
        .err(sprintf("the target values for %s must add up to 1",
                     add_quotes(v)))
      }
    }

    if (!all(is.na(internal.targets)) && any_apply(formula.vars, function(v) {
      if (attr(terms(formula.covs), "order")[formula.vars == v] <= 1) {
        return(FALSE)
      }

      vars.in.interaction <- rownames(attr(terms(formula.covs), "factors"))[attr(terms(formula.covs), "factors")[, v] == 1]

      sum(attr(terms(formula.covs), "dataClasses")[vars.in.interaction] == "factor") > 1
    })) {
      .wrn("interactions between factor variables were entered, but `process_targets()` cannot verify whether the target values are suitable. See `?check_targets` for details")
    }

    attr(internal.targets, "original.vars") <- formula.vars[attr(model.covs, "assign")] |> setNames(model.vars)
  }

  attr(internal.targets, "ATE") <- ATE

  class(internal.targets) <- "optweight.targets"

  internal.targets
}

#' @export
#' @rdname process_targets
check.targets <- function(...) {
  process_targets(...)
}

#' @exportS3Method print optweight.targets
#' @rdname process_targets
print.optweight.targets <- function(x, digits = 5, ...) {
  targets <- x
  attributes(targets) <- NULL
  names(targets) <- names(x)

  if (all(is.na(targets))) {
    cat0(" - ", .it("variables"), ":\n\t",
         paste(names(targets), collapse = space(3L)))
  }
  else {
    cat0(" - ", .it("targets"), ":\n")
    print(round(targets, digits))
  }

  invisible(x)
}
