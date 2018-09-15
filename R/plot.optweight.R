plot.optweight <- function(x, which.time = 1, ...) {
  if ("optweightMSM" %in% class(x)) {
    if (which.time %nin% seq_along(x$duals)) stop("which.time must correspond to an available time point.", call. = FALSE)
    d <- x$duals[[which.time]]
    title <- paste("Dual Variables for Constraints at Time", which.time)
  }
  else {
    d <- x$duals
    title <- "Dual Variables for Constraints"
  }

  d$cov <- factor(d$cov, levels = rev(unique(d$cov)))
  d$treat <- factor(d$treat, levels = unique(d$treat), labels = paste("Treat:", unique(d$treat)))
  d$constraint <- factor(d$constraint, levels = unique(d$constraint, nmax=2), labels = paste("Constraint:", unique(d$constraint, nmax=2)))

  p <- ggplot(d, aes(x = cov, y = dual)) +
    geom_col() +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(y = "Absolute Dual Variable",
         x = "Covariate",
         title = title) +
    theme_bw() +
    scale_y_continuous(expand = expand_scale(c(0, .05)))

  p <- p + facet_grid(cols = vars(treat),
                      rows = vars(constraint),
                      scales = "free_y", space = "free")
  p
}
