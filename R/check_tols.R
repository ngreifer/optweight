#' Construct and Check Tolerance Input
#'
#' Checks whether proposed tolerance values for `tols` are suitable in number and order for submission to [optweight()]. Users should include one value per item in `formula`. The output can also be used as an input to `tols` in [optweight()].
#'
#' @param formula A formula with the covariates to be balanced with [optweight()] on the right-hand side. See [glm()] for more details. Interactions and functions of covariates are allowed. Lists of formulas are not allowed; multiple formulas must be checked one at a time.
#' @param data An optional data set in the form of a data frame that contains the variables in `formula`.
#' @param tols A vector of balance tolerance values in standardized mean difference units for each covariate. These should be in the order corresponding to the order of the corresponding variable in `formula`, except for interactions, which will appear after all lower-order terms. If only one value is supplied, it will be applied to all covariates.
#' @param stop Logical; if `TRUE`, an error will be thrown if the number of values in `tols` is not equal to the correct number of covariates in `formula`, and no messages will be displayed if the `tols` input is satisfactory. If `FALSE`, a message will be displayed if the number of values in `tols` is not equal to the correct number of covariates in `formula`, and other messages will be displayed.
#' @param x An `optweight.tols` object; the output of a call to `check_tols()`.
#' @param internal Logical; whether to print the tolerance values that are to be used internally by [optweight()]. See Value section.
#' @param digits How many digits to print.
#' @param ... Ignored.
#'
#' @returns
#' An `optweight.tols` object, which is a named vector of tolerance values, one for each variable specified in `formula`. This should be used as user inputs to [optweight()]. The `"internal.tols"` attribute contains the tolerance values to be used internally by [optweight()]. These will differ from the vector values when there are factor variables that are split up; the user only needs to submit one tolerance per factor variable, but separate tolerance values are produced for each new dummy created.
#'
#' @details
#' The purpose of `check_tols()` is to allow users to ensure that their proposed input to `tols` in [optweight()] is correct both in the number of entries and their order. This is especially important when factor variables and interactions are included in the formula because factor variables are split into several dummies and interactions are moved to the end of the variable list, both of which can cause some confusion and potential error when entering `tols` values.
#'
#' Factor variables are internally split into a dummy variable for each level, but the user only needs to specify one tolerance value per original variable; `check_tols()` automatically expands the `tols` input to match the newly created variables.
#'
#' Interactions (e.g., `a:b` or `a*b` in the `formula` input) are always sent to the end of the variable list even if they are specified elsewhere in the `formula`. It is important to run `check_tols()` to ensure the order of the proposed `tols` corresponds to the represented order of covariates used in [optweight()]. You can run `check_tols()` with no `tols` input to see the order of covariates that is required.
#'
#' `check_tols()` was designed to be used primarily for its message printing and `print` method, but you can also assign its output to an object for use as an input to `tols` in [optweight()].
#'
#' Note that only one formula and vector of tolerance values can be assessed at a time; for multiple treatment periods, each formula and tolerance vector must be entered separately.
#'
#' @seealso [check_targets()]
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' # Checking if the correct number of entries are included:
#' check_tols(treat ~ age + educ + married +
#'                 nodegree + re74, data = lalonde,
#'                 tols = c(.01, .02, .03, .04))
#'
#' # Checking the order of interactions; notice they go
#' # at the end even if specified at the beginning. The
#' # .09 values are where the interactions might be expected
#' # to be, but they are in fact not.
#' c <- check_tols(treat ~ age:educ + married*race +
#'                 nodegree + re74, data = lalonde,
#'                 tols = c(.09, .01, .01, .09, .01, .01))
#'
#' print(c, internal = TRUE)

#' @export
check_tols <- function(formula, data = NULL, tols, stop = FALSE) {

  if (!rlang::is_formula(formula)) {
    .err("the argument to `formula` must a single formula with the covariates on the right side")
  }

  if (missing(tols)) {
    tols <- NA_real_
  }
  else if (!is.numeric(tols)) {
    .err("`tols` must be a numeric vector")
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

  if (length(tols) == 1) {
    onetol <- TRUE
    tols <- rep_with(tols, formula.vars)
  }
  else {
    onetol <- FALSE
  }

  if (is_null(model.covs)) {
    .err("no covariates were specified")
  }

  check_missing_covs(formula.covs)

  if (length(tols) != length(formula.vars)) {
    if (stop) {
      .err(sprintf("`tols` must contain %s numbers; %s were included",
                   length(formula.vars), length(tols)))
    }

    .msg(sprintf("`tols` must contain %s numbers; %s were included.\nAll values in the output have been assigned NA",
                 length(formula.vars), length(tols)))

    user.tols <- rep_with(NA_real_, formula.vars)
    internal.tols <- rep_with(NA_real_, model.vars)
  }
  else {
    user.tols <- tols |> setNames(formula.vars)
    internal.tols <- tols[attr(model.covs, "assign")] |> setNames(model.vars)
    if (!stop) {
      if (any(attr(terms(formula), "order") > 1)) {
        #message("tols look okay, but interactions were present in the formula, so make sure the order is correct.")
      }
      else {
        #message("tols look okay, but you should print them to make sure they're correct.")
      }
    }
  }

  attr(user.tols, "internal.tols") <- internal.tols
  attr(user.tols, "original.vars") <- formula.vars[attr(model.covs, "assign")] |> setNames(model.vars)
  class(user.tols) <- "optweight.tols"

  user.tols
}

#' @export
#' @rdname check_tols
check.tols <- function(...) {
  check_tols(...)
}

#' @rdname check_tols
#' @exportS3Method print optweight.tols
print.optweight.tols <- function(x, internal = FALSE, digits = 5, ...) {
  tols <- x
  internal.tols <- attr(tols, "internal.tols")
  attributes(tols) <- NULL
  names(tols) <- names(x)
  if (all(is.na(tols))) {
    cat(sprintf("- vars:\n\t%s",
                paste(names(tols), collapse = space(3L))))
  }
  else {
    cat("- tols:\n")
    print(round(tols, digits))
    if (internal) {
      cat("\n- tols used internally by optweight:\n")
      print(round(internal.tols, digits))
    }
  }
  invisible(x)
}


