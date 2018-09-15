relax <- function(ow, step = .001, ...) {
  if ("optweight" %nin% class(ow)) {
    stop("The first argument must be an optweight object.", call. = FALSE)
  }
  if (ow$info$status_val == 1) {
    message("The optweight object seems to have already converged.")
    return(ow)
  }
  else {
    ow.obj <- as.list(ow$call)

    tols <- ow$tols
    form <- eval(ow.obj$formula)
    data <- eval(ow.obj$data)

    if (is.list(form)) times <- length(form)
    else times <- 1

    if (!is.list(tols)) tols.list <- list(tols)
    else tols.list <- tols

    if (!is.list(step)) step.list <- list(step)
    else step.list <- step
    if (length(step.list) == 1) step.list <- replicate(max(times), step.list[[1]], simplify = FALSE)

    for (i in times) {
      step[[i]] <- check.tols(form, data, tols = step, stop = TRUE)
    }
  }


  return(new.ow)
}
