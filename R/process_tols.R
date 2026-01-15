#' Construct and Check Tolerance Input
#'
#' Checks whether proposed tolerance values for `tols` are suitable in number and order for submission to [optweight()] and [optweight.svy()], and returns an object that can supplied to the `tols` argument of these functions.
#'
#' @inheritParams optweight
#' @param formula a formula with the covariates to be balanced on the right-hand side. Interactions and functions of covariates are allowed. Lists of formulas are not allowed; multiple formulas must be checked one at a time.
#' @param tols a vector of balance tolerance values in standardized mean difference units for each covariate. These should be in the order corresponding to the order of the corresponding variable in `formula`, except for interactions, which will appear after all lower-order terms. If only one value is supplied, it will be applied to all covariates.
#' @param x an `optweight.tols` object; the output of a call to `process_tols()`.
#' @param internal `logical`; whether to print the tolerance values that are to be used internally by [optweight()]. See Value section.
#' @param digits how many digits to print.
#' @param ... ignored.
#'
#' @returns
#' An `optweight.tols` object, which is a named vector of tolerance values, one for each variable specified in `formula`. This should be used as an input to the `tols` argument of [optweight()]. The `"internal.tols"` attribute contains the tolerance values to be used internally by [optweight()]. These will differ from the vector values when there are factor variables that are split up; the user only needs to submit one tolerance per factor variable, but separate tolerance values are produced for each new dummy created.
#'
#' @details
#' The purpose of `process_tols()` is to allow users to ensure that their proposed input to `tols` in [optweight()] is correct both in the number of entries and their order. This is especially important when factor variables and interactions are included in the formula because factor variables are split into several dummies and interactions are moved to the end of the variable list, both of which can cause some confusion and potential error when entering `tols` values.
#'
#' Factor variables are internally split into a dummy variable for each level, but the user only needs to specify one tolerance value per original variable; `process_tols()` automatically expands the `tols` input to match the newly created variables.
#'
#' Interactions (e.g., `a:b` or `a*b` in the `formula` input) are always sent to the end of the variable list even if they are specified elsewhere in the `formula`. It is important to run `process_tols()` to ensure the order of the proposed `tols` corresponds to the represented order of covariates used in [optweight()]. You can run `process_tols()` with no `tols` input to see the order of covariates that is required.
#'
#' Note that only one formula and vector of tolerance values can be assessed at a time; for multiple treatments, each formula and tolerance vector must be entered separately.
#'
#' @seealso [process_targets()]
#'
#' @examplesIf rlang::is_installed("cobalt")
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' # Generating tols; 0 by default
#' tols <- process_tols(treat ~ age + educ + married +
#'                        nodegree + re74,
#'                      data = lalonde)
#'
#' tols
#'
#' tols <- process_tols(treat ~ age + educ + married +
#'                        nodegree + re74,
#'                      data = lalonde,
#'                      tols = .05)
#'
#' tols
#'
#' # Checking the order of interactions; notice they go
#' # at the end even if specified at the beginning.
#' tols <- process_tols(treat ~ age:educ + married*race +
#'                        nodegree + re74,
#'                      data = lalonde,
#'                      tols = .05)
#'
#' tols
#'
#' # Internal tolerances for expanded covariates
#' print(tols, internal = TRUE)

#' @export
process_tols <- function(formula, data = NULL, tols = 0) {

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
    .err("the argument to {.arg formula} must a single formula with the covariates on the right side")
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

  .process_tols_internal(covs, tols, reported.covs,
                         if (formula.present) "formula" else "data")
}

.process_tols_internal <- function(model.covs, tols, formula.covs = NULL,
                                   tols_found_in = "formula",
                                   tols_arg = "tols") {

  chk::chk_numeric(tols)
  chk::chk_not_any_na(tols)

  model.vars <- colnames(model.covs)

  if (is_not_null(formula.covs)) {
    formula.vars <- .attr(formula.covs, "terms") |> .attr("term.labels")
    if (is_null(formula.vars)) {
      formula.vars <- "(Intercept)"
      attr(model.covs, "assign") <- 1
    }
  }
  else {
    formula.vars <- model.vars
  }

  if (is_null(names(tols))) {
    if (length(tols) == 1L) {
      tols <- setNames(rep_with(tols, formula.vars),
                       formula.vars)
    }
    else if (length(tols) == length(formula.vars)) {
      names(tols) <- formula.vars
    }
    else {
      .err("{.arg {tols_arg}} must contain {.or {unique(c(1, length(formula.vars)))}}
            {cli::qty(max(1, length(formula.vars)))} value{?s}, but {length(tols)} {?was/were} included")
    }
  }

  if (!any(names(tols) %in% formula.vars)) {
    .err("no variables named in {.arg {tols_arg}} are present in {.arg {tols_found_in}}")
  }

  overlap <- intersect(names(tols), formula.vars)

  user.tols <- setNames(rep_with(0, formula.vars),
                        formula.vars)
  user.tols[overlap] <- tols[overlap]

  if (!all(names(tols) %in% overlap)) {
    bad_tols <- add_quotes(setdiff(names(tols), overlap))
    .wrn("{bad_tols} {?was/were} named in {.arg {tols_arg}} but not present in {.arg {tols_found_in}} and so will be ignored")
  }

  if (is_not_null(formula.covs)) {
    internal.tols <- user.tols[.attr(model.covs, "assign")]
    names(internal.tols) <- model.vars

    attr(user.tols, "original.vars") <- names(user.tols)[.attr(model.covs, "assign")] |> setNames(model.vars)
  }
  else {
    internal.tols <- user.tols
  }

  attr(user.tols, "internal.tols") <- internal.tols

  class(user.tols) <- "optweight.tols"

  user.tols
}

#' @export
#' @rdname process_tols
check.tols <- function(...) {
  process_tols(...)
}

#' @exportS3Method print optweight.tols
#' @rdname process_tols
print.optweight.tols <- function(x, internal = FALSE, digits = 5, ...) {
  tols <- x
  internal.tols <- .attr(tols, "internal.tols")
  attributes(tols) <- NULL
  names(tols) <- names(x)

  if (allNA(tols)) {
    cat0(" - ", .it("variables"), ":\n\t",
         paste(names(tols), collapse = space(3L)))
  }
  else {
    cat0(" - ", .it("tols"), ":\n")
    print(round(tols, digits))

    if (internal && is_not_null(internal.tols)) {
      cat0("\n - ", .it("tols used internally by optweight"), ":\n")
      print(round(internal.tols, digits))
    }
  }

  invisible(x)
}
