plot.optweight <- function(x, which.time = 1, scale = c("var", "ess"), ...) {
  if ("optweightMSM" %in% class(x)) {
    if (which.time %nin% seq_along(x$duals)) stop("which.time must correspond to an available time point.", call. = FALSE)
    d <- x$duals[[which.time]]
    title <- paste("Dual Variables for Balance Constraints at Time", which.time)
    #t <- x$treat.list[[which.time]]
  }
  else {
    d <- x$duals
    title <- "Dual Variables for Balance Constraints"
    #t <- x$treat
  }

  #w <- x$weights
  # if (scale == "ess") {
  #
  # }
  # duals <- data.frame(reshape(d, direction = "long", ids = rownames(d),
  #                             idvar = "covs", varying = list(colnames(d)),
  #                             timevar = "treat", times = colnames(d), v.names = "dual"))
  #
  # duals$covs <- factor(duals$covs, levels = rev(unique(duals$covs)))
  # duals$treat <- factor(duals$treat); levels(duals$treat) <- paste("Treat =", levels(duals$treat))

  d$cov <- factor(d$cov, levels = rev(unique(d$cov)))

  p <- ggplot(d, aes(x = cov, y = dual)) +
    geom_col() +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(y = "Absolute Standardized Dual Variable",
         x = "Covariate",
         title = title) +
    theme_bw() +
    scale_y_continuous(expand = expand_scale(c(0, .05)))

  p <- p + facet_grid(cols = vars(treat),
                                                  rows = vars(constraint),
                                                  scales = "free_y", space = "free")
  #else p <- p + facet_grid(rows = vars(constraint), scales = "free_y", space = "free")
  p
}
