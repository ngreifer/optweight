#' Construct and Check Targets Input
#'
#' Checks whether proposed target population means values for `targets` are suitable in number and order for submission to [optweight()] and [optweight.svy()]. Users should include one value per variable in `formula`. For factor variables, one value per level of the variable is required. The output of `check_targets()` can also be used as an input to `targets` in [optweight()] and [optweight.svy()].
#'
#' @param formula A formula with the covariates to be balanced with [optweight()] on the right-hand side. See [glm()] for more details. Interactions and functions of covariates are allowed.
#' @param data An optional data set in the form of a data frame that contains the variables in `formula`.
#' @param targets A vector of target population mean values for each covariate. These should be in the order corresponding to the order of the corresponding variable in `formula`, except for interactions, which will appear after all lower-order terms. For factor variables, a target value must be specified for each level of the factor, and these values must add up to 1. If empty, the current sample means will be produced. If `NULL`, an `NA` vector named with the covariate names will be produced.
#' @param stop Logical; if `TRUE`, an error will be thrown if the number of values in `targets` is not equal to the correct number of (expanded) covariates in `formula`, and no messages will be displayed if the `targets` input is satisfactory. If `FALSE`, a message will be displayed if the number of values in `targets` is not equal to the correct number of covariates in `formula`, and other messages will be displayed.
#' @param x An `optweight.targets` object; the output of a call to `check_targets()`.
#' @param digits How many digits to print.
#' @param ... Ignored.
#'
#' @returns
#' An `optweight.targets` object, which is a named vector of target population mean values, one for each (expanded) covariate specified in `formula`. This should be used as user inputs to [optweight()] and [optweight.svy()].
#'
#' @details
#' The purpose of `check_targets()` is to allow users to ensure that their proposed input to `targets` in [optweight()] and [optweight.svy()] is correct both in the number of entries and their order. This is especially important when factor variables and interactions are included in the formula because factor variables are split into several dummies and interactions are moved to the end of the variable list, both of which can cause some confusion and potential error when entering `targets` values.
#'
#' Factor variables are internally split into a dummy variable for each level, so the user must specify a target population mean value for each level of the factor. These must add up to 1, and an error will be displayed if they do not. These values represent the proportion of units in the target population with each factor level.
#'
#' Interactions (e.g., `a:b` or `a*b` in the `formula` input) are always sent to the end of the variable list even if they are specified elsewhere in the `formula`. It is important to run `check_targets()` to ensure the order of the proposed `targets` corresponds to the represented order of covariates used in the formula. You can run `check_targets(targets = NULL)` to see the order of covariates that is required without specifying any targets.
#'
#' @seealso [check_tols()]
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' # Checking if the correct number of entries are included:
#' check_targets(treat ~ age + race + married +
#'                 nodegree + re74,
#'               data = lalonde,
#'               targets = c(25, .4, .1, .5, .3,
#'                           .5, 4000))
#' # Notice race is split into three values (.4, .1, and .5)
#'

#' @export
check_targets <- function(formula, data = NULL, targets, stop = FALSE) {
  if (!rlang::is_formula(formula)) {
    .err("the argument to formula must a single formula with the covariates on the right side")
  }

  #Process treat and covs from formula and data
  tt <- delete.response(terms(formula))
  t.c <- get_covs_and_treat_from_formula2(tt, data, terms = TRUE, sep = "_")
  formula.covs <- t.c[["reported.covs"]]
  model.covs <- t.c[["model.covs"]]

  formula.vars <- attr(attr(formula.covs, "terms"), "term.labels")
  if (is_null(formula.vars)) {
    formula.vars <- "(Intercept)"
    attr(model.covs, "assign") <- 1
  }
  model.vars <- colnames(model.covs)

  if (is_null(model.covs)) {
    .err("no covariates were specified")
  }

  check_missing_covs(formula.covs)

  missing.targets <- missing(targets)

  if (missing.targets) {
    targets <- colMeans(model.covs)
  }
  else if (is_null(targets)) {
    targets <- rep_with(NA_real_, model.vars)
  }
  else if (!is.numeric(targets)) {
    .err("`targets` must be a numeric vector")
  }

  if (length(targets) != length(model.vars)) {
    if (stop) {
      .err(sprintf("`targets` must contain %s numbers, but %s were included",
                   length(model.vars), length(targets)))
    }

    .msg(sprintf("`targets` must contain %s numbers, but %s were included.\nAll values in the output have been assigned NA",
                 length(model.vars), length(targets)))

    internal.targets <- rep_with(NA_real_, model.vars)
  }

  internal.targets <- targets |> setNames(model.vars)
  original.variables <- formula.vars[attr(model.covs, "assign")] |> setNames(model.vars)

  for (v in formula.vars) {
    if (attr(terms(formula.covs), "order")[formula.vars == v] == 1 &&
        attr(terms(formula.covs), "dataClasses")[formula.vars == v] == "factor" &&
        !any(is.na(internal.targets[original.variables == v])) &&
        !check_if_zero(sum(internal.targets[original.variables == v]) - 1)) {
      .err(sprintf("the target values for %s must add up to 1",
                   add_quotes(v)))
    }
  }

  if (!all(is.na(targets)) && any_apply(formula.vars, function(v) {
    if (attr(terms(formula.covs), "order")[formula.vars == v] <= 1) {
      return(FALSE)
    }

    vars.in.interaction <- rownames(attr(terms(formula.covs), "factors"))[attr(terms(formula.covs), "factors")[, v] == 1]

    sum(attr(terms(formula.covs), "dataClasses")[vars.in.interaction] == "factor") > 1
  })) {
    .wrn("interactions between factor variables were entered, but `check_targets()` cannot verify whether the target values are suitable. See `?check_targets` for details")
  }

  if (!stop) {
    if (any(attr(terms(formula), "order") > 1)) {
      #message("targets look okay, but interactions were present in the formula, so make sure the order is correct.")
    }
    else {
      #message("targets look okay, but you should print them to make sure they're correct.")
    }
  }

  out <- internal.targets

  attr(out, "original.vars") <- formula.vars[attr(model.covs, "assign")] |> setNames(model.vars)

  attr(out, "ATE") <- (missing.targets || isTRUE(all(check_if_zero(colMeans(model.covs) - internal.targets))))

  class(out) <- "optweight.targets"

  out
}

#' @export
#' @rdname check_targets
check.targets <- function(...) {
  check_targets(...)
}

#' @rdname check_targets
#' @exportS3Method print optweight.targets
print.optweight.targets <- function(x, digits = 5, ...) {
  targets <- x
  attributes(targets) <- NULL
  names(targets) <- names(x)

  if (all(is.na(targets))) {
    cat(sprintf("- vars:\n\t%s",
                paste(names(targets), collapse = space(3L))))
  }
  else {
    cat("- targets:\n")
    print(round(targets, digits))
  }

  invisible(x)
}
