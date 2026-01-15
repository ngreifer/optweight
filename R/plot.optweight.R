#' Plot Dual Variables for Covariate Constraints
#'
#' Plots the dual variables resulting from [optweight()], [optweightMV()], or [optweight.svy()] in a way similar to figure 2 of Zubizarreta (2015), which explains how to interpret these values.
#'
#' @param x an `optweight`, `optweightMV`, or `optweight.svy` object; the output of a call to [optweight()], [optweightMV()], or [optweight.svy()].
#' @param which.treat for `optweightMV` objects, an integer corresponding to which treatment to display. Only one may be displayed at a time.
#' @param type the type of plot to display; allowable options include `"variables"` (the default), which produces a row for each covariate, and `"constraints"`, which produces a row for each type of constraint (computed as the sum of the absolute dual variables for each constraint type).
#' @param \dots ignored.
#'
#' @returns
#' A `ggplot` object that can be used with other \pkg{ggplot2} functions.
#'
#' @details
#' Dual variables represent the cost of changing the constraint on the objective function minimized to estimate the weights. For covariates with large values of the dual variable, tightening the constraint will increase the variability of the weights, and relaxing the constraint will decrease the variability of the weights, both to a greater extent than would doing the same for covariate with small values of the dual variable. See [optweight()] and `vignette("optweight")` for more information on interpreting dual variables.
#'
#' @seealso
#' [optweight()], [optweightMV()], or [optweight.svy()] to estimate the weights and the dual variables.
#'
#' [plot.summary.optweight()] for plots of the distribution of weights.
#'
#' @references
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' @examplesIf rlang::is_installed("cobalt")
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' tols <- process_tols(treat ~ age + educ + married +
#'                        nodegree + re74, data = lalonde,
#'                      tols = .1)
#'
#' #Balancing covariates between treatment groups (binary)
#' ow1 <- optweight(treat ~ age + educ + married +
#'                    nodegree + re74, data = lalonde,
#'                  tols = tols,
#'                  estimand = "ATT")
#'
#' # Note the L2 divergence and effective sample
#' # size (ESS)
#' summary(ow1, weight.range = FALSE)
#'
#' # age has a low value, married is high
#' plot(ow1)
#'
#' tols["age"] <- 0
#' ow2 <- optweight(treat ~ age + educ + married +
#'                    nodegree + re74, data = lalonde,
#'                  tols = tols,
#'                  estimand = "ATT")
#'
#' # Notice that tightening the constraint on age has
#' # a negligible effect on the variability of the
#' # weights and ESS
#' summary(ow2, weight.range = FALSE)
#'
#' tols["age"] <- .1
#' tols["married"] <- 0
#' ow3 <- optweight(treat ~ age + educ + married +
#'                    nodegree + re74, data = lalonde,
#'                  tols = tols,
#'                  estimand = "ATT")
#'
#' # In contrast, tightening the constraint on married
#' # has a large effect on the variability of the
#' # weights, shrinking the ESS
#' summary(ow3, weight.range = FALSE)
#'
#' # More duals are displayed when targeting other
#' # estimands:
#' ow4 <- optweight(treat ~ age + educ + married +
#'                    nodegree + re74, data = lalonde,
#'                  estimand = "ATE")
#'
#' plot(ow4)
#'
#' # Display duals by constraint type
#' plot(ow4, type = "constraints")

#' @exportS3Method plot optweight
plot.optweight <- function(x, type = "variables", ...) {
  chk::chk_string(type)
  type <- match_arg(type, c("variables", "constraints"))

  title <- "Dual Variables for Constraints"

  .plot_optweight_internal(x$duals, title, type)
}

#' @rdname plot.optweight
#' @exportS3Method plot optweightMV
plot.optweightMV <- function(x, which.treat = 1L, type = "variables", ...) {
  chk::chk_count(which.treat)

  if (which.treat %nin% seq_len(max(x$duals$component))) {
    .err("{.arg which.treat} must correspond to an available treatment")
  }

  chk::chk_string(type)
  type <- match_arg(type, c("variables", "constraints"))

  title <- sprintf("Dual Variables for Constraints for Treatment %s",
                   which.treat)

  duals <- ss(x$duals, x$duals$component %in% c(0L, which.treat))

  .plot_optweight_internal(duals, title, type)
}

#' @rdname plot.optweight
#' @exportS3Method plot optweight.svy
plot.optweight.svy <- plot.optweight

.plot_optweight_internal <- function(d, title, type) {

  if (type == "variables") {
    constraint_types <- c("target", "balance")

    d <- ss(d, d$constraint %in% constraint_types)

    d$cov <- factor(d$cov, levels = unique(rev(d$cov)))

    d$constraint <- factor(d$constraint,
                           levels = constraint_types,
                           labels = paste("Constraint:", constraint_types))

    p <- ggplot(d) +
      geom_col(aes(y = .data$cov, x = .data$dual)) +
      scale_x_continuous(expand = expansion(c(0, .05))) +
      facet_grid(rows = vars(.data$constraint),
                 scales = "free_y", space = "free") +
      labs(x = "Absolute Dual Variable",
           y = "Variable",
           title = title) +
      theme_bw()
  }
  else {
    constraint_types <- c("target", "balance", "weight range")

    d$constraint <- factor(d$constraint,
                           levels = rev(constraint_types))

    agg <- collap(d, dual ~ constraint, FUN = sum, sort = FALSE)

    p <- ggplot(agg) +
      geom_col(aes(y = .data$constraint, x = .data$dual)) +
      scale_x_continuous(expand = expansion(c(0, .05))) +
      labs(x = "Absolute Dual Variable",
           y = "Constraint",
           title = title) +
      theme_bw()
  }

  p
}
