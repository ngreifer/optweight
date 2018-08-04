plot.optweight <- function(x, which.time = 1, ...) {
  if ("optweightMSM" %in% class(x)) {
    d <- x$duals[[which.time]]
    title <- paste("Dual Variables for Balance Constraints at Time", which.time)
  }
  else {
    d <- x$duals
    title <- "Dual Variables for Balance Constraints"
  }

  duals <- data.frame(reshape(d, direction = "long", ids = rownames(d),
                              idvar = "covs", varying = list(colnames(d)),
                              timevar = "treat", times = colnames(d), v.names = "dual"))

  duals$covs <- factor(duals$covs, labels = rev(unique(duals$covs)))
  duals$treat <- factor(duals$treat); levels(duals$treat) <- paste("Treat =", levels(duals$treat))

  p <- ggplot(duals, aes(x = covs, y = dual)) +
    geom_col() +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(y = "Absolute Standardized Dual Variable",
         x = "Covariate",
         title = title) +
    theme_bw() +
    scale_y_continuous(expand = expand_scale(c(0, .05)))

  if (nlevels(duals$treat) > 1) p <- p + facet_grid(cols = vars(treat))
  p
}
