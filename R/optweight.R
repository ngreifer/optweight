optweight <- function(formula, data = NULL, tols = 0, estimand = "ATE", s.weights = NULL, focal = NULL, verbose = FALSE, ...) {

  if (!is.list(formula)) formula.list <- list(formula)
  else formula.list <- formula
  times <- seq_along(formula.list)
  onetime <- length(formula.list) == 1

  if (!identical(tols, 0)) {
    if (class(tols) == "optweight.tols") tols.list <- list(tols[["internal.tols"]])
    else if (is.atomic(tols)) tols.list <- list(tols)
    else tols.list <- tols
    if (length(tols.list) == 1) tols.list <- replicate(max(times), tols.list[[1]], simplify = FALSE)
    exact <- FALSE
  }
  else {
    tols.list <- NULL
    exact <- TRUE
  }

  reported.covs.list <- covs.list <- treat.list <- vector("list", length(formula.list))
  n <- 0 * times
  for (i in times) {
    #Process treat and covs from formula and data
    t.c <- get.covs.and.treat.from.formula(formula.list[[i]], data)
    reported.covs.list[[i]] <- t.c[["reported.covs"]]
    covs.list[[i]] <- t.c[["model.covs"]]
    treat.list[[i]] <- t.c[["treat"]]
    #treat.name <- t.c[["treat.name"]]

    #Get treat type
    treat.list[[i]] <- get.treat.type(treat.list[[i]])
    treat.type <- attr(treat.list[[i]], "treat.type")

    if (onetime) {
      if (is_null(covs.list[[i]])) stop("No covariates were specified.", call. = FALSE)
      if (is_null(treat.list[[i]])) stop("No treatment variable was specified.", call. = FALSE)
      if (any(is.na(treat.list[[i]]))) {
        stop("No missing values are allowed in the treatment variable.", call. = FALSE)
      }

      #Process estimand and focal
      f.e.r <- process.focal.and.estimand(focal, estimand, treat.list[[i]], treat.type)
      focal <- f.e.r[["focal"]]
      estimand <- f.e.r[["estimand"]]
      reported.estimand <- f.e.r[["reported.estimand"]]
    }
    else {
      if (is_null(covs.list[[i]])) stop(paste0("No covariates were specified in formula ", i, "."), call. = FALSE)
      if (is_null(treat.list[[i]])) stop(paste0("No treatment variable was specified in formula ", i, "."), call. = FALSE)
      if (any(!is.finite(treat.list[[i]]))) {
        stop(paste0("No missing or non-finite values are allowed in the treatment variable. Missing or non-finite values were found in treatment ", i, "."), call. = FALSE)
      }
      if (toupper(estimand) != "ATE") stop("The only estimand allowed with longitudinal treatments is the ATE.", call. = FALSE)
      focal <- NULL
      estimand <- "ATE"
      reported.estimand <- "ATE"
    }

    n[i] <- length(treat.list[[i]])

    if (any(bad.covs <- !sapply(reported.covs.list[[i]], is.finite)) || nrow(reported.covs.list[[i]]) != n[i]) {
      stop(paste0("No missing or non-finite values are allowed in the covariates. Missing or non-finite values were found in the following covariates:\n", paste(names(reported.covs.list[[i]])[bad.covs], collapse = ", ")), call. = FALSE)
    }

    if (!exact) {
      tryCatch(ct <- check.tols(formula.list[[i]], reported.covs.list[[i]], tols.list[[i]], stop = TRUE),
               error = function(e) {
                 if (onetime) e. <- conditionMessage(e)
                 else e. <- paste0("For treatment ", i, ", ", conditionMessage(e))
                 stop(e., call. = FALSE)})
      tols.list[[i]] <- ct[["internal.tols"]]
    }
  }

  if (!all_the_same(n)) stop("The same number of units must be present in each time point.", call. = FALSE)

  #Process s.weights
  sw <- process.s.weights(s.weights, data)

  ###Run optweight.fit
  w <- optweight.fit(treat = treat.list, covs = covs.list,
                     tols = tols.list, estimand = estimand,
                     focal = focal,
                     s.weights = sw,
                     std.binary = FALSE,
                     exact = exact,
                     verbose = verbose,
                     ...)

  warn <- FALSE
  test.w <- if (is_null(sw)) w else w*sw
  if (any(sapply(treat.list, function(t) attr(t, "treat.type") == "continuous"))) {if (sd(test.w)/mean(test.w) > 4) warn <- TRUE}
  else if (any(sapply(treat.list, function(t) sapply(unique(t), function(x) sd(test.w[t == x])/mean(test.w[t == x]) > 4)))) warn <- TRUE
  if (warn) warning("Some extreme weights were generated. Examine them with summary() and maybe trim them with trim().", call. = FALSE)
  call <- match.call()

  if (onetime) {
    out <- list(weights = w,
                treat = treat.list[[1]],
                covs = reported.covs.list[[1]],
                s.weights = sw,
                estimand = reported.estimand,
                focal = focal,
                call = call)

    class(out) <- "optweight"
  }
  else {
    out <- list(weights = w,
                treat.list = treat.list,
                covs.list = reported.covs.list,
                s.weights = sw,
                estimand = reported.estimand,
                call = call)

    class(out) <- c("optweightMSM", "optweight")
  }
  return(out)
}

print.optweight <- function(x, ...) {
  cat("An optweight object\n")
  cat(paste0(" - number of obs.: ", length(x[["weights"]]), "\n"))
  if (is_not_null(x[["estimand"]])) cat(paste0(" - estimand: ", x[["estimand"]], ifelse(is_not_null(x[["focal"]]), paste0(" (focal: ", x[["focal"]], ")"), ""), "\n"))
  cat(paste0(" - sampling weights: ", ifelse(all_the_same(x[["s.weights"]]),"none", "present"), "\n"))
  invisible(x)
}
print.optweightMSM <- function(x, ...) {
  cat("An optweightMSM object\n")
  cat(paste0(" - number of obs.: ", length(x[["weights"]]), "\n"))
  cat(paste0(" - number of time points: ", length(x[["treat.list"]]), "\n"))
  cat(paste0(" - sampling weights: ", ifelse(all_the_same(x[["s.weights"]]),"none", "present"), "\n"))
  invisible(x)
}
