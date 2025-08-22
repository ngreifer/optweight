#' Summarize, print, and plot information about estimated weights
#'
#' These functions summarize the weights resulting from a call to
#' [optweight()] or [optweight.svy()]. `summary()` produces summary statistics on the distribution of weights, including their
#' range and variability, and the effective sample size of the weighted sample
#' (computing using the formula in McCaffrey, Rudgeway, & Morral, 2004). `plot()` creates a histogram of the weights.
#'
#' @param object An `optweight`, `optweightMSM`, or `optweight.svy` object; the output of a call to [optweight()] or [optweight.svy()].
#' @param top How many of the largest and smallest weights to display. Default
#' is 5.
#' @param ignore.s.weights Whether or not to ignore sampling weights when
#' computing the weight summary. If `FALSE`, the default, the estimated
#' weights will be multiplied by the sampling weights (if any) before values
#' are computed.
#' @param x A `summary.optweight`, `summary.optweightMSM`, or `summary.optweight.svy` object; the output of a call to `summary.optweight()`, `summary.optweightMSM()`, or ()`summary.optweight.svy`.
#' @param ...  Additional arguments. For `plot()`, additional arguments
#' passed to [graphics::hist()] to determine the number of bins,
#' though [ggplot2::geom_histogram()] from \pkg{ggplot2} is actually
#' used to create the plot.
#'
#' @returns
#' For point treatments (i.e., `optweight` objects),
#' `summary()` returns a `summary.optweight` object with the following
#' elements:
#' \item{weight.range}{The range (minimum and maximum) weight for
#' each treatment group.}
#' \item{weight.top}{The units with the greatest weights
#' in each treatment group; how many are included is determined by `top`.}
#' \item{rms.dev}{The root-mean-squared deviation of the estimated weights from the base weights (L2 norm), weighted by the sampling weights (if any).}
#' \item{mean.abs.dev}{The mean absolute deviation of the estimated weights from the base weights (L1 norm), weighted by the sampling weights (if any).}
#' \item{max.abs.dev}{The maximum absolute deviation of the estimated weights from the base weights (L\eqn{\infinity} norm).}
#' \item{num.zeros}{The number of units with a weight equal to 0.}
#' \item{effective.sample.size}{The effective sample size for each treatment group before and after weighting.}
#'
#' For longitudinal treatments (i.e., `optweightMSM` objects), a list of
#' the above elements for each treatment period.
#'
#' For `optweight.svy` objects, a list of the above elements but with no
#' treatment group divisions.
#'
#' `plot()` returns a `ggplot` object with a histogram displaying the
#' distribution of the estimated weights. If the estimand is the ATT or ATC,
#' only the weights for the non-focal group(s) will be displayed (since the
#' weights for the focal group are all 1). A dotted line is displayed at the
#' mean of the weights (usually 1).
#'
#' @seealso
#' [plot.optweight()] for plotting the values of the dual variables.
#'
#' @references
#' McCaffrey, D. F., Ridgeway, G., & Morral, A. R. (2004).
#' Propensity Score Estimation With Boosted Regression for Evaluating Causal
#' Effects in Observational Studies. *Psychological Methods*, 9(4), 403–425. \doi{10.1037/1082-989X.9.4.403}
#'
#' @examplesIf requireNamespace("cobalt", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (ow1 <- optweight(treat ~ age + educ + married +
#'                 nodegree + re74, data = lalonde,
#'                 tols = .001,
#'                 estimand = "ATT"))
#'
#' (s <- summary(ow1))
#'
#' plot(s, breaks = 12)
#'

#' @exportS3Method summary optweight
summary.optweight <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  sw <- {
    if (ignore.s.weights || is_null(object$s.weights)) rep_with(1, object$weights)
    else object$s.weights
  }

  bw <- {
    if (is_null(object$b.weights)) rep_with(1, object$weights)
    else object$b.weights
  }

  t <- object$treat
  treat.type <- attr(t, "treat.type")

  out <- .summary_internal(t, treat.type, object$weights, sw, bw, top)

  w <- object$weights * sw
  if (is_not_null(object$focal)) {
    w <- w[t != object$focal]
    attr(w, "focal") <- object$focal
  }

  attr(out, "weights") <- w

  class(out) <- "summary.optweight"
  out
}

#' @exportS3Method print summary.optweight
print.summary.optweight <- function(x, ...) {
  cat("Summary of weights:\n\n")

  .print_summary_internal(x, ...)

  invisible(x)
}

#' @rdname summary.optweight
#' @exportS3Method summary optweightMSM
summary.optweightMSM <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  .wrn("optweights are currently not valid for longitudinal treatments. Use at your own risk")

  sw <- {
    if (ignore.s.weights || is_null(object$s.weights)) rep_with(1, object$weights)
    else object$s.weights
  }

  bw <- {
    if (is_null(object$b.weights)) rep_with(1, object$weights)
    else object$b.weights
  }

  treat.types <- vapply(object[["treat.list"]], attr, character(1L), "treat.type")

  out.list <- make_list(names(object$treat.list))
  for (ti in seq_along(object$treat.list)) {
    out.list[[ti]] <- .summary_internal(object$treat.list[[ti]],
                                        treat.types[ti], object$weights, sw, bw, top)
  }

  attr(out.list, "weights") <- object$weights * sw

  class(out.list) <- c("summary.optweightMSM", "summary.optweight")

  out.list
}

#' @exportS3Method print summary.optweightMSM
print.summary.optweightMSM <- function(x, ...) {
  only.one <- length(x) == 1L || all_apply(x, identical, x[[1L]])

  cat("Summary of weights:\n\n")
  for (ti in seq_along(x)) {
    if (!only.one) {
      cat(sprintf(" - - - - - - - - - - Time %s - - - - - - - - - -\n", ti))
    }

    .print_summary_internal(x[[ti]], ...)

    if (only.one) {
      break
    }
  }

  invisible(x)
}

#' @rdname summary.optweight
#' @exportS3Method summary optweight
summary.optweight.svy <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  sw <- {
    if (ignore.s.weights || is_null(object$s.weights)) rep_with(1, object$weights)
    else object$s.weights
  }

  bw <- {
    if (is_null(object$b.weights)) rep_with(1, object$weights)
    else object$b.weights
  }

  out <- .summary_internal(NULL, "svy", object$weights, sw, bw, top)

  attr(out, "weights") <- object$weights * sw

  class(out) <- c("summary.optweight.svy", "summary.optweight")

  out
}

#' @exportS3Method print summary.optweight.svy
print.summary.optweight.svy <- function(x, ...) {
  cat("Summary of weights:\n\n")

  .print_summary_internal(x, ...)

  invisible(x)
}

#' @rdname summary.optweight
#' @exportS3Method plot summary.optweight
plot.summary.optweight <- function(x, ...) {
  w <- attr(x, "weights")
  focal <- attr(w, "focal")

  subtitle <- if (is_not_null(focal)) sprintf("For Units Not in Treatment Group %s", add_quotes(focal))

  ggplot(mapping = aes(x = w)) +
    geom_histogram(breaks = graphics::hist(w, plot = FALSE, ...)$breaks,
                   color = "black",
                   fill = "gray", alpha = .8) +
    scale_y_continuous(expand = expansion(c(0, .05))) +
    geom_vline(xintercept = mean(w), linetype = "12") +
    labs(x = "Weight", y = "Count", title = "Distribution of Weights",
         subtitle = subtitle) +
    theme_bw()
}

.summary_internal <- function(t, treat.type, w, sw, bw, top) {
  outnames <- c("weight.range", "weight.top", "weight.ratio",
                "rms.dev", "mean.abs.dev", "max.abs.dev", "num.zeros",
                "effective.sample.size")

  out <- make_list(outnames)

  ww <- w * sw

  if (treat.type %in% c("continuous", "svy")) {
    out$weight.range <- list(all = c(min(ww[ww != 0]),
                                     max(ww[ww != 0])))
    top.weights <- sort(ww, decreasing = TRUE)[seq_len(top)]
    out$weight.top <- list(all = sort(setNames(top.weights, which(ww %in% top.weights)[seq_len(top)])))
    out$rms.dev <- c(all = rms_dev(w, bw = bw, sw = sw))
    out$mean.abs.dev <- c(all = mean_abs_dev(w, bw = bw, sw = sw))
    out$max.abs.dev <- c(all = max_abs_dev(w, bw = bw, sw = sw))
    out$num.zeros <- c(all = sum(ww == 0))

    nn <- make_df("Total", c("Unweighted", "Weighted"))
    nn[["Total"]] <- c(ESS(sw), ESS(ww))
  }
  else if (treat.type == "binary") {
    out$weight.range <- list(treated = c(min(ww[ww != 0 & t == 1]),
                                         max(ww[ww != 0 & t == 1])),
                             control = c(min(ww[ww != 0 & t == 0]),
                                         max(ww[ww != 0 & t == 0])))
    top.weights <- list(treated = sort(ww[t == 1], decreasing = TRUE)[seq_len(top)],
                        control = sort(ww[t == 0], decreasing = TRUE)[seq_len(top)])

    out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(ww[t == switch(x, control = 0, 1)] %in% top.weights[[x]])[seq_len(top)]))),
                               names(top.weights))

    out$rms.dev <- c(treated = rms_dev(w[t==1], bw = bw[t==1], sw = sw[t==1]),
                     control = rms_dev(w[t==0], bw = bw[t==0], sw = sw[t==0]))
    out$mean.abs.dev <- c(treated = mean_abs_dev(w[t==1], bw = bw[t==1], sw = sw[t==1]),
                          control = mean_abs_dev(w[t==0], bw = bw[t==0], sw = sw[t==0]))
    out$max.abs.dev <- c(treated = max_abs_dev(w[t==1], bw = bw[t==1], sw = sw[t==1]),
                         control = max_abs_dev(w[t==0], bw = bw[t==0], sw = sw[t==0]))
    out$num.zeros <- c(treated = sum(ww[t==1] == 0),
                       control = sum(ww[t==0] == 0))

    #dc <- weightit$discarded

    nn <- make_df(c("Control", "Treated"), c("Unweighted", "Weighted"))

    nn[["Control"]] <- c(ESS(sw[t == 0]), ESS(ww[t == 0]))
    nn[["Treated"]] <- c(ESS(sw[t == 1]), ESS(ww[t == 1]))
  }
  else if (treat.type == "multi-category") {
    out$weight.range <- setNames(lapply(levels(t), function(x) c(min(ww[ww != 0 & t == x]),
                                                                 max(ww[ww != 0 & t == x]))),
                                 levels(t))

    top.weights <- setNames(lapply(levels(t), function(x) sort(ww[t == x], decreasing = TRUE)[seq_len(top)]),
                            levels(t))

    out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(ww[t == x] %in% top.weights[[x]])[seq_len(top)]))),
                               names(top.weights))

    out$rms.dev <- vapply(levels(t), function(x) rms_dev(w[t==x], bw = bw[t==x], sw = sw[t==x]), numeric(1L))
    out$mean.abs.dev <- vapply(levels(t), function(x) mean_abs_dev(w[t==x], bw = bw[t==x], sw = sw[t==x]), numeric(1L))
    out$max.abs.dev <- vapply(levels(t), function(x) max_abs_dev(w[t==x], bw = bw[t==x], sw = sw[t==x]), numeric(1L))
    out$num.zeros <- vapply(levels(t), function(x) sum(ww[t==x] == 0), numeric(1L))

    nn <- make_df(levels(t), c("Unweighted", "Weighted"))

    for (i in levels(t)) {
      nn[[i]] <- c(ESS(sw[t == i]), ESS(ww[t == i]))
    }
  }

  out$effective.sample.size <- nn

  out
}

.print_summary_internal <- function(x, ...) {
  cat("- Weight ranges:\n")

  x$weight.range |>
    text_box_plot(width = 28L) |>
    round_df_char(digits = 4) |>
    print.data.frame(...)

  cat(sprintf("\n- Units with %s greatest weights by group:\n",
              length(x$weight.top[[1L]])))

  top <- max(lengths(x$weight.top))
  data.frame(unlist(lapply(names(x$weight.top), function(y) c(" ", y))),
             matrix(unlist(lapply(x$weight.top, function(y) c(names(y), character(top - length(y)),
                                                              round(y, 4), character(top - length(y))))),
                    byrow = TRUE, nrow = 2 * length(x$weight.top))) |>
    setNames(character(1 + top)) |>
    print.data.frame(row.names = FALSE)

  cat("\n")

  matrix(c(x$rms.dev, x$mean.abs.dev,
           x$max.abs.dev, x$num.zeros),
         nrow = length(x$rms.dev),
         byrow = FALSE,
         dimnames = list(names(x$rms.dev),
                         c("RMSE Dev", "Mean Abs Dev", "Max Abs Dev", "# Zeros"))) |>
    as.data.frame() |>
    round_df_char(digits = 4) |>
    print.data.frame(...)

  cat("\n- Effective Sample Sizes:\n")
  x$effective.sample.size |>
    round_df_char(digits = 3) |>
    print.data.frame(...)
  cat("\n")

  invisible(x)
}
