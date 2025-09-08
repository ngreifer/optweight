#' Plot Dual Variables for Assessing Balance Constraints
#'
#' Plots the dual variables resulting from [optweight()] in a way similar to
#' figure 2 of Zubizarreta (2015), which explained how to interpret these
#' values. These represent the cost of changing the constraint on the variance
#' of the resulting weights. For covariates with large values of the dual
#' variable, tightening the constraint will increase the variability of the
#' weights, and loosening the constraint will decrease the variability of the
#' weights, both to a greater extent than would doing the same for covariate
#' with small values of the dual variable.
#'
#' @param x an `optweight`, `optweightMV`, or `optweight.svy` object; the output of a
#' call to [optweight()], [optweightMV()], or [optweight.svy()].
#' @param which.treat For `optweightMV` objects, which treatment to display.
#' Only one may be displayed at a time.
#' @param \dots Ignored.
#'
#' @returns
#' A `ggplot` object that can be used with other \pkg{ggplot2} functions.
#'
#' @seealso [optweight()], [optweightMV()], or [optweight.svy()] to estimate
#' the weights and the dual variables
#'
#' [plot.summary.optweight()] for plots of the distribution of
#' weights
#'
#' @references
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
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
#' summary(ow1) # Note the RMSE Dev and effective
#' #              sample size (ESS)
#'
#' plot(ow1) # age has a low value, married is high
#'
#' tols["age"] <- 0
#' ow2 <- optweight(treat ~ age + educ + married +
#'                    nodegree + re74, data = lalonde,
#'                  tols = tols,
#'                  estimand = "ATT")
#'
#' summary(ow2) # Notice that tightening the constraint
#' #              on age had a negligible effect on the
#' #              variability of the weights and ESS
#'
#' tols["age"] <- .1
#' tols["married"] <- 0
#' ow3 <- optweight(treat ~ age + educ + married +
#'                    nodegree + re74, data = lalonde,
#'                  tols = tols,
#'                  estimand = "ATT")
#'
#' summary(ow3) # In contrast, tightening the constraint
#' #              on married had a large effect on the
#' #              variability of the weights, shrinking
#' #              the ESS


#' @exportS3Method plot optweight
plot.optweight <- function(x, ...) {
  d <- x$duals
  title <- "Dual Variables for Constraints"

  d$cov <- factor(d$cov, levels = rev(unique(d$cov)))
  d$constraint <- factor(d$constraint, levels = unique(d$constraint, nmax = 2L),
                         labels = paste("Constraint:", unique(d$constraint, nmax = 2L)))

  ggplot(d, mapping = aes(x = .data$cov, y = .data$dual)) +
    geom_col() +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(y = "Absolute Dual Variable",
         x = "Covariate",
         title = title) +
    theme_bw() +
    scale_y_continuous(expand = expansion(c(0, .05))) +
    facet_grid(rows = vars(.data$constraint),
               scales = "free_y", space = "free")
}

#' @rdname plot.optweight
#' @exportS3Method plot optweightMV
plot.optweightMV <- function(x, which.treat = 1, ...) {
  chk::chk_count(which.treat)

  if (which.treat %nin% seq_along(x$duals)) {
    .err("`which.treat` must correspond to an available treatment")
  }

  d <- x$duals[[which.treat]]
  title <- sprintf("Dual Variables for Constraints for Treatment %s", which.treat)

  d$cov <- factor(d$cov, levels = rev(unique(d$cov)))
  d$constraint <- factor(d$constraint, levels = unique(d$constraint, nmax = 2L),
                         labels = paste("Constraint:", unique(d$constraint, nmax = 2L)))

  ggplot(d, mapping = aes(x = .data$cov, y = .data$dual)) +
    geom_col() +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(y = "Absolute Dual Variable",
         x = "Covariate",
         title = title) +
    theme_bw() +
    scale_y_continuous(expand = expansion(c(0, .05))) +
    facet_grid(rows = vars(.data$constraint),
               scales = "free_y", space = "free")
}

#' @rdname plot.optweight
#' @exportS3Method plot optweight.svy
plot.optweight.svy <- function(x, ...) {
  d <- x$duals
  title <- "Dual Variables for Target Constraints"

  d$cov <- factor(d$cov, levels = rev(unique(d$cov)))

  p <- ggplot(d, mapping = aes(x = .data$cov, y = .data$dual)) +
    geom_col() +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(y = "Absolute Dual Variable",
         x = "Covariate",
         title = title) +
    theme_bw() +
    scale_y_continuous(expand = expansion(c(0, .05)))

  p
}
