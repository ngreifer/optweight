summary.optweight <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  outnames <- c("weight.range", "weight.top","weight.ratio",
                "coef.of.var",
                "effective.sample.size")
  out <- setNames(vector("list", length(outnames)), outnames)

  if (ignore.s.weights || is_null(object$s.weights)) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights
  w <- object$weights*sw
  t <- object$treat
  treat.type <- attr(object[["treat"]], "treat.type")

  if (treat.type == "continuous") {
    out$weight.range <- list(all = c(min(w[w > 0]),
                                     max(w[w > 0])))
    top.weights <- sort(w, decreasing = TRUE)[seq_len(top)]
    out$weight.top <- list(all = sort(setNames(top.weights, which(w %in% top.weights)[seq_len(top)])))
    out$coef.of.var <- c(all = coef.of.var(w))
    out$mean.abs.dev <- c(all = mean.abs.dev(w))

    nn <- as.data.frame(matrix(0, ncol = 1, nrow = 2))
    nn[1, ] <- ESS(sw)
    nn[2, ] <- ESS(w)
    dimnames(nn) <- list(c("Unweighted", "Weighted"),
                         c("Total"))

  }
  else if (treat.type == "binary") {
    top0 <- c(treated = min(top, sum(t == 1)),
              control = min(top, sum(t == 0)))
    out$weight.range <- list(treated = c(min(w[w > 0 & t == 1]),
                                         max(w[w > 0 & t == 1])),
                             control = c(min(w[w > 0 & t == 0]),
                                         max(w[w > 0 & t == 0])))
    top.weights <- list(treated = sort(w[t == 1], decreasing = TRUE)[seq_len(top0["treated"])],
                        control = sort(w[t == 0], decreasing = TRUE)[seq_len(top0["control"])])
    out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w %in% top.weights[[x]] & t == {if (x == "control") 0 else 1})[seq_len(top0[x])]))),
                               names(top.weights))
    top.weights <- list(treated = sort(w[t == 1], decreasing = TRUE)[seq_len(top0["treated"])],
                        control = sort(w[t == 0], decreasing = TRUE)[seq_len(top0["control"])])

    out$coef.of.var <- c(treated = coef.of.var(w[t==1]),
                         control = coef.of.var(w[t==0]),
                         overall = coef.of.var(w))
    out$mean.abs.dev <- c(treated = mean.abs.dev(w[t==1]),
                          control = mean.abs.dev(w[t==0]),
                          overall = mean.abs.dev(w))

    #dc <- weightit$discarded

    nn <- as.data.frame(matrix(0, nrow = 2, ncol = 2))
    nn[1, ] <- c(ESS(sw[t==0]), ESS(sw[t==1]))
    nn[2, ] <- c(ESS(w[t==0]), ESS(w[t==1]))
    # nn[3, ] <- c(sum(t==0 & dc==1), #Discarded
    #              sum(t==1 & dc==1))
    dimnames(nn) <- list(c("Unweighted", "Weighted"),
                         c("Control", "Treated"))
  }
  else if (treat.type == "multinomial") {
    out$weight.range <- setNames(lapply(levels(t), function(x) c(min(w[w > 0 & t == x]),
                                                                 max(w[w > 0 & t == x]))),
                                 levels(t))
    top.weights <- setNames(lapply(levels(t), function(x) sort(w[t == x], decreasing = TRUE)[seq_len(top)]),
                            levels(t))
    out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w %in% top.weights[[x]] & t == x)[seq_len(top)]))),
                               names(top.weights))
    out$coef.of.var <- c(sapply(levels(t), function(x) coef.of.var(w[t==x])),
                         overall = coef.of.var(w))
    out$mean.abs.dev <- c(sapply(levels(t), function(x) mean.abs.dev(w[t==x])),
                          overall = mean.abs.dev(w))

    nn <- as.data.frame(matrix(0, nrow = 2, ncol = nunique(t)))
    for (i in seq_len(nunique(t))) {
      nn[1, i] <- ESS(sw[t==levels(t)[i]])
      nn[2, i] <- ESS(w[t==levels(t)[i]])
      # nn[3, i] <- sum(t==levels(t)[i] & dc==1) #Discarded
    }
    dimnames(nn) <- list(c("Unweighted", "Weighted"),
                         levels(t))
  }

  out$effective.sample.size <- nn

  if (is_not_null(object$focal)) {
    w <- w[t != object$focal]
    attr(w, "focal") <- object$focal
  }
  attr(out, "weights") <- w

  class(out) <- "summary.optweight"
  return(out)
}
print.summary.optweight <- function(x, ...) {
  top <- max(lengths(x$weight.top))
  cat("Summary of weights:\n\n")
  cat("- Weight ranges:\n")
  print.data.frame(round_df_char(text_box_plot(x$weight.range, 28), 4), ...)
  df <- setNames(data.frame(do.call("c", lapply(names(x$weight.top), function(x) c(" ", x))),
                            matrix(do.call("c", lapply(x$weight.top, function(x) c(names(x), rep("", top - length(x)), round(x, 4), rep("", top - length(x))))),
                                   byrow = TRUE, nrow = 2*length(x$weight.top))),
                 rep("", 1 + top))
  cat(paste("\n- Units with", top, "greatest weights by group:\n"))
  print.data.frame(df, row.names = FALSE)
  cat("\n")
  print.data.frame(round_df_char(as.data.frame(matrix(c(x$coef.of.var, x$mean.abs.dev), ncol = 2,
                                                      byrow = FALSE,
                                                      dimnames = list(names(x$coef.of.var),
                                                                      c("Coef of Var", "Mean Abs Dev")))), 4))
  cat("\n- Effective Sample Sizes:\n")
  print.data.frame(round_df_char(x$effective.sample.size, 3))
  invisible(x)
}

summary.optweightMSM <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  warning("Optweights are currently not valid for longitudinal treatments. Use at your own risk.", call. = FALSE)

  outnames <- c("weight.range", "weight.top","weight.ratio",
                "coef.of.var", "weight.mean",
                "effective.sample.size")
  out.list <- setNames(vector("list", length(object$treat.list)),
                       names(object$treat.list))

  if (ignore.s.weights || is_null(object$s.weights)) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights
  w <- object$weights*sw
  treat.types <- sapply(object[["treat.list"]], function(y) attr(y, "treat.type"))

  for (ti in seq_along(object$treat.list)) {
    if (treat.types[ti] == "continuous") {
      out <- setNames(vector("list", length(outnames)), outnames)
      out$weight.range <- list(all = c(min(w[w > 0]),
                                       max(w[w > 0])))
      top.weights <- sort(w, decreasing = TRUE)[seq_len(top)]
      out$weight.top <- list(all = sort(setNames(top.weights, which(w %in% top.weights)[seq_len(top)])))
      out$coef.of.var <- c(all = coef.of.var(w))
      out$mean.abs.dev <- c(all = mean.abs.dev(w))

      nn <- as.data.frame(matrix(0, ncol = 1, nrow = 2))
      nn[1, ] <- ESS(sw)
      nn[2, ] <- ESS(w)
      dimnames(nn) <- list(c("Unweighted", "Weighted"),
                           c("Total"))
      out$effective.sample.size <- nn

      out.list[[ti]] <- out

    }
    else if (treat.types[ti] == "binary") {
      out <- setNames(vector("list", length(outnames)), outnames)
      t <- object$treat.list[[ti]]
      out$weight.range <- list(treated = c(min(w[w > 0 & t == 1]),
                                           max(w[w > 0 & t == 1])),
                               control = c(min(w[w > 0 & t == 0]),
                                           max(w[w > 0 & t == 0])))
      top.weights <- list(treated = sort(w[t == 1], decreasing = TRUE)[seq_len(top)],
                          control = sort(w[t == 0], decreasing = TRUE)[seq_len(top)])
      out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w[t == ifelse(x == "control", 0, 1)] %in% top.weights[[x]])[seq_len(top)]))),
                                 names(top.weights))
      out$coef.of.var <- c(treated = coef.of.var(w[t==1]),
                           control = coef.of.var(w[t==0]),
                           overall = coef.of.var(w))
      out$mean.abs.dev <- c(treated = mean.abs.dev(w[t==1]),
                            control = mean.abs.dev(w[t==0]),
                            overall = mean.abs.dev(w))

      #dc <- weightit$discarded

      nn <- as.data.frame(matrix(0, nrow = 2, ncol = 2))
      nn[1, ] <- c(ESS(sw[t==0]), ESS(sw[t==1]))
      nn[2, ] <- c(ESS(w[t==0]), ESS(w[t==1]))
      # nn[3, ] <- c(sum(t==0 & dc==1), #Discarded
      #              sum(t==1 & dc==1))
      dimnames(nn) <- list(c("Unweighted", "Weighted"),
                           c("Control", "Treated"))
      out$effective.sample.size <- nn
      out.list[[ti]] <- out

    }
    else if (treat.types[ti] == "multinomial") {

      out <- setNames(vector("list", length(outnames)), outnames)
      t <- object$treat.list[[ti]]
      out$weight.range <- setNames(lapply(levels(t), function(x) c(min(w[w > 0 & t == x]),
                                                                   max(w[w > 0 & t == x]))),
                                   levels(t))
      top.weights <- setNames(lapply(levels(t), function(x) sort(w[t == x], decreasing = TRUE)[seq_len(top)]),
                              levels(t))
      out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w[t == x] %in% top.weights[[x]])[seq_len(top)]))),
                                 names(top.weights))
      out$coef.of.var <- c(sapply(levels(t), function(x) coef.of.var(w[t==x])),
                           overall = coef.of.var(w))
      out$mean.abs.dev <- c(sapply(levels(t), function(x) mean.abs.dev(w[t==x])),
                            overall = mean.abs.dev(w))

      nn <- as.data.frame(matrix(0, nrow = 2, ncol = nunique(t)))
      for (i in seq_len(nunique(t))) {
        nn[1, i] <- ESS(sw[t==levels(t)[i]])
        nn[2, i] <- ESS(w[t==levels(t)[i]])
        # nn[3, i] <- sum(t==levels(t)[i] & dc==1) #Discarded
      }
      dimnames(nn) <- list(c("Unweighted", "Weighted"),
                           levels(t))
      out$effective.sample.size <- nn
      out.list[[ti]] <- out
    }

  }

  attr(out.list, "weights") <- w

  class(out.list) <- c("summary.optweightMSM", "summary.optweight")
  return(out.list)
}
print.summary.optweightMSM <- function(x, ...) {
  if (all(vapply(x, function(y) isTRUE(all.equal(x[[1]], y)), logical(1L)))) {
    only.one <- TRUE
  }
  else only.one <- FALSE

  cat("Summary of weights:\n\n")
  for (ti in seq_along(x)) {
    if (!only.one) cat(paste(" - - - - - - - - - - Time", ti, "- - - - - - - - - -\n"))
    cat("- Weight ranges:\n")
    print.data.frame(round_df_char(text_box_plot(x[[ti]]$weight.range, 28), 4))

    df <- setNames(data.frame(do.call("c", lapply(names(x[[ti]]$weight.top), function(y) c(" ", y))),
                              matrix(do.call("c", lapply(x[[ti]]$weight.top, function(y) c(names(y), round(y, 4)))),
                                     byrow = TRUE, nrow = 2*length(x[[ti]]$weight.top))),
                   rep("", 1 + length(x[[ti]]$weight.top[[1]])))
    cat(paste("\n- Units with", length(x[[ti]]$weight.top[[1]]), "greatest weights by group:\n"))
    print.data.frame(df, row.names = FALSE)
    cat("\n")
    print.data.frame(round_df_char(as.data.frame(matrix(c(x[[ti]]$coef.of.var, x[[ti]]$mean.abs.dev), ncol = 2,
                                                        byrow = FALSE,
                                                        dimnames = list(names(x[[ti]]$coef.of.var),
                                                                        c("Coef of Var", "Mean Abs Dev")))), 4))

    cat("\n- Effective Sample Sizes:\n")
    print.data.frame(round_df_char(x[[ti]]$effective.sample.size, 3))
    cat("\n")
    if (only.one) break
  }

  invisible(x)
}

summary.optweight.svy <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  outnames <- c("weight.range", "weight.top","weight.ratio",
                "coef.of.var",
                "effective.sample.size")
  out <- setNames(vector("list", length(outnames)), outnames)

  if (ignore.s.weights || is_null(object$s.weights)) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights
  w <- object$weights*sw

    out$weight.range <- list(all = c(min(w[w > 0]),
                                     max(w[w > 0])))
    top.weights <- sort(w, decreasing = TRUE)[seq_len(top)]
    out$weight.top <- list(all = sort(setNames(top.weights, which(w %in% top.weights)[seq_len(top)])))
    out$coef.of.var <- c(all = coef.of.var(w))
    out$mean.abs.dev <- c(all = mean.abs.dev(w))

    nn <- as.data.frame(matrix(0, ncol = 1, nrow = 2))
    nn[1, ] <- ESS(sw)
    nn[2, ] <- ESS(w)
    dimnames(nn) <- list(c("Unweighted", "Weighted"),
                         c("Total"))

  out$effective.sample.size <- nn

  attr(out, "weights") <- w

  class(out) <- c("summary.optweight.svy", "summary.optweight")
  return(out)
}
print.summary.optweight.svy <- function(x, ...) {
  top <- max(lengths(x$weight.top))
  cat("Summary of weights:\n\n")
  cat("- Weight ranges:\n")
  print.data.frame(round_df_char(text_box_plot(x$weight.range, 28), 4), ...)
  df <- setNames(data.frame(do.call("c", lapply(names(x$weight.top), function(x) c(" ", x))),
                            matrix(do.call("c", lapply(x$weight.top, function(x) c(names(x), rep("", top - length(x)), round(x, 4), rep("", top - length(x))))),
                                   byrow = TRUE, nrow = 2*length(x$weight.top))),
                 rep("", 1 + top))
  cat(paste("\n- Units with", top, "greatest weights:\n"))
  print.data.frame(df, row.names = FALSE)
  cat("\n")
  print.data.frame(round_df_char(as.data.frame(matrix(c(x$coef.of.var, x$mean.abs.dev), ncol = 2,
                                                      byrow = FALSE,
                                                      dimnames = list(names(x$coef.of.var),
                                                                      c("Coef of Var", "Mean Abs Dev")))), 4))
  cat("\n- Effective Sample Sizes:\n")
  print.data.frame(round_df_char(x$effective.sample.size, 3))
  invisible(x)
}

plot.summary.optweight <- function(x, ...) {
  w <- attr(x, "weights")
  focal <- attr(w, "focal")

  if (is_not_null(focal)) subtitle <- paste0("For Units Not in Treatment Group \"", focal, "\"")
  else subtitle <- NULL

  p <- ggplot(data = data.frame(w), mapping = aes(x = w)) +
    geom_histogram(breaks = hist(w, plot = FALSE,
                                 ...)$breaks,
                   color = "black",
                   fill = "gray", alpha = .8) +
    scale_y_continuous(expand = expand_scale(c(0, .05))) +
    geom_vline(xintercept = mean(w), linetype = "12") +
    labs(x = "Weight", y = "Count", title = "Distribution of Weights",
         subtitle = subtitle) +
    theme_bw()
  p
}
