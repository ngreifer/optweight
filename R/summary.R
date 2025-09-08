#' Summarize, print, and plot information about estimated weights
#'
#' These functions summarize the weights resulting from a call to
#' [optweight()] or [optweight.svy()]. `summary()` produces summary statistics on the distribution of weights, including their
#' range and variability, and the effective sample size of the weighted sample
#' (computing using the formula in McCaffrey, Rudgeway, & Morral, 2004). `plot()` creates a histogram of the weights.
#'
#' @param object An `optweight`, `optweightMV`, or `optweight.svy` object; the output of a call to [optweight()] or [optweight.svy()].
#' @param top How many of the largest and smallest weights to display. Default
#' is 5.
#' @param ignore.s.weights Whether or not to ignore sampling weights when
#' computing the weight summary. If `FALSE`, the default, the estimated
#' weights will be multiplied by the sampling weights (if any) before values
#' are computed.
#' @param weight.range `logical`; whether display statistics about the range of weights and the highest and lowest weights for each group. Default is `TRUE`.
#' @param x A `summary.optweight`, `summary.optweightMV`, or `summary.optweight.svy` object; the output of a call to `summary.optweight()`, `summary.optweightMV()`, or ()`summary.optweight.svy`.
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
#' \item{rel.ent}{The relative entropy between the estimated weights and the base weights (entropy norm), weighted by the sampling weights (if any). Only computed if all weights are positive.}
#' \item{num.zeros}{The number of units with a weight equal to 0.}
#' \item{effective.sample.size}{The effective sample size for each treatment group before and after weighting.}
#'
#' For multivariate treatments (i.e., `optweightMV` objects), a list of
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
#'                     nodegree + re74, data = lalonde,
#'                   tols = .001,
#'                   estimand = "ATT"))
#'
#' (s <- summary(ow1))
#'
#' plot(s, breaks = 12)

#' @exportS3Method summary optweight
summary.optweight <- function(object, top = 5, ignore.s.weights = FALSE, weight.range = TRUE, ...) {

  chk::chk_flag(ignore.s.weights)
  chk::chk_flag(weight.range)

  if (weight.range) {
    chk::chk_count(top)
  }

  sw <- {
    if (ignore.s.weights || is_null(object$s.weights)) rep_with(1, object$weights)
    else object$s.weights
  }

  bw <- {
    if (is_null(object$b.weights)) rep_with(1, object$weights)
    else object$b.weights
  }

  t <- object$treat

  out <- .summary_internal(t, attr(t, "treat.type"),
                           object$weights, sw, bw,
                           top, weight.range)

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
#' @exportS3Method summary optweightMV
summary.optweightMV <- function(object, top = 5, ignore.s.weights = FALSE, weight.range = TRUE, ...) {
  chk::chk_flag(ignore.s.weights)
  chk::chk_flag(weight.range)

  if (weight.range) {
    chk::chk_count(top)
  }

  sw <- {
    if (ignore.s.weights || is_null(object$s.weights)) rep_with(1, object$weights)
    else object$s.weights
  }

  bw <- {
    if (is_null(object$b.weights)) rep_with(1, object$weights)
    else object$b.weights
  }

  out.list <- lapply(object$treat.list, function(t) {
    .summary_internal(t, attr(t, "treat.type"),
                      object$weights, sw, bw,
                      top, weight.range)
  })

  attr(out.list, "weights") <- object$weights * sw

  class(out.list) <- c("summary.optweightMV", "summary.optweight")

  out.list
}

#' @exportS3Method print summary.optweightMV
print.summary.optweightMV <- function(x, ...) {
  only.one <- length(x) == 1L || all_apply(x, identical, x[[1L]])

  cat("Summary of weights:\n\n")
  for (ti in seq_along(x)) {
    if (!only.one) {
      cat(sprintf(" - - - - - - - - - - Treatment %s - - - - - - - - - -\n", ti))
    }

    .print_summary_internal(x[[ti]], ...)

    if (only.one) {
      break
    }
  }

  invisible(x)
}

#' @rdname summary.optweight
#' @exportS3Method summary optweight.svy
summary.optweight.svy <- function(object, top = 5, ignore.s.weights = FALSE, weight.range = TRUE, ...) {
  chk::chk_flag(ignore.s.weights)
  chk::chk_flag(weight.range)

  if (weight.range) {
    chk::chk_count(top)
  }

  sw <- {
    if (ignore.s.weights || is_null(object$s.weights)) rep_with(1, object$weights)
    else object$s.weights
  }

  bw <- {
    if (is_null(object$b.weights)) rep_with(1, object$weights)
    else object$b.weights
  }

  out <- .summary_internal(NULL, "svy", object$weights, sw, bw, top, weight.range)

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

.summary_internal <- function(t, treat.type, w, sw, bw, top, weight.range) {
  outnames <- c("weight.range", "weight.top", "weight.ratio",
                "rms.dev", "mean.abs.dev", "max.abs.dev", "rel.ent", "num.zeros",
                "effective.sample.size")

  out <- make_list(outnames)

  treat.type[treat.type == "multinomial"] <- "multi-category"

  ww <- w * sw

  if (treat.type == "binary") {
    tx <- list(treated = which(t == 1),
               control = which(t == 0))
  }
  else if (treat.type == "multi-category") {
    tx <- lapply(levels(t), function(i) which(t == i)) |>
      setNames(levels(t))
  }
  else {
    tx <- list(all = seq_along(w))
  }

  if (weight.range) {
    out$weight.range <- lapply(tx, function(ti) c(min(ww[ti]), max(w[ti]))) |>
      setNames(names(tx))

    top.weights <- lapply(tx, function(ti) sort(ww[ti], decreasing = TRUE)[seq_len(top)]) |>
      setNames(names(tx))

    out$weight.top <- lapply(names(tx), function(i) {
      sort(setNames(top.weights[[i]], which(ww[tx[[i]]] %in% top.weights[[i]])[seq_len(top)]))
    }) |>
      setNames(names(tx))
  }

  out$rms.dev <- vapply(tx, function(ti) rms_dev(w[ti], bw = bw[ti], sw = sw[ti]), numeric(1L))
  out$mean.abs.dev <- vapply(tx, function(ti) mean_abs_dev(w[ti], bw = bw[ti], sw = sw[ti]), numeric(1L))
  out$max.abs.dev <- vapply(tx, function(ti) max_abs_dev(w[ti], bw = bw[ti], sw = sw[ti]), numeric(1L))
  out$rel.ent <- if (all(w > 0)) vapply(tx, function(ti) rel_ent(w[ti], bw = bw[ti], sw = sw[ti]), numeric(1L))
  out$num.zeros <- vapply(tx, function(ti) sum(ww[ti] == 0), numeric(1L))

  if (treat.type == "binary") {
    nn <- make_df(c("Control", "Treated"),
                  c("Unweighted", "Weighted"))

    nn[["Control"]] <- c(ESS(sw[tx$control]), ESS(ww[tx$control]))
    nn[["Treated"]] <- c(ESS(sw[tx$treated]), ESS(ww[tx$treated]))
  }
  else if (treat.type == "multi-category") {
    nn <- make_df(levels(t),
                  c("Unweighted", "Weighted"))

    for (i in levels(t)) {
      nn[[i]] <- c(ESS(sw[tx[[i]]]), ESS(ww[tx[[i]]]))
    }
  }
  else {
    nn <- make_df("Total",
                  c("Unweighted", "Weighted"))
    nn[["Total"]] <- c(ESS(sw), ESS(ww))
  }

  out$effective.sample.size <- nn

  out
}

.print_summary_internal <- function(x, ...) {
  if (is_not_null(x$weight.range)) {
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
      setNames(character(1L + top)) |>
      print.data.frame(row.names = FALSE)

    cat("\n")
  }

  matrix(c(x$rms.dev, x$mean.abs.dev,
           x$max.abs.dev, x$rel.ent, x$num.zeros),
         nrow = length(x$rms.dev),
         byrow = FALSE,
         dimnames = list(names(x$rms.dev),
                         c("RMSE Dev", "Mean Abs Dev", "Max Abs Dev",
                           "Rel Ent"[is_not_null(x$rel.ent)],
                           "# Zeros"))) |>
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
