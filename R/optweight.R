optweight <- function(formula, data = NULL, tols = .0001, estimand = "ATE", s.weights = NULL, focal = NULL) {

  #Process treat and covs from formula and data
  t.c <- get.covs.and.treat.from.formula(formula, data)
  reported.covs <- t.c[["reported.covs"]]
  covs <- t.c[["model.covs"]]
  treat <- t.c[["treat"]]
  treat.name <- t.c[["treat.name"]]

  if (is_null(covs)) stop("No covariates were specified.", call. = FALSE)
  if (is_null(treat)) stop("No treatment variable was specified.", call. = FALSE)

  n <- length(treat)

  if (any(!sapply(reported.covs, is.finite)) || nrow(reported.covs) != n) {
    stop("No missing or non-finite values are allowed in the covariates.", call. = FALSE)
  }
  if (any(is.na(treat))) {
    stop("No missing values are allowed in the treatment variable.", call. = FALSE)
  }

  #Get treat type
  treat <- get.treat.type(treat)
  treat.type <- attr(treat, "treat.type")

  #Process estimand and focal
  estimand <- process.estimand(estimand, method, treat.type)
  f.e.r <- process.focal.and.estimand(focal, estimand, treat, treat.type)
  focal <- f.e.r[["focal"]]
  estimand <- f.e.r[["estimand"]]
  reported.estimand <- f.e.r[["reported.estimand"]]

  #Process s.weights
  sw <- process.s.weights(s.weights, data)

  ###Run optweight.fit
  w <- optweight.fit(treat = treat, covs = covs,
                       tols = tols, estimand = estimand,
                       focal = focal,
                       s.weights = sw, std.binary = FALSE)


  warn <- FALSE
  test.w <- if (is_null(sw)) w else w*sw
  if (treat.type == "continuous") {if (sd(test.w)/mean(test.w) > 4) warn <- TRUE}
  else if (any(sapply(unique(treat), function(x) sd(test.w[treat == x])/mean(test.w[treat == x]) > 4))) warn <- TRUE
  if (warn) warning("Some extreme weights were generated. Examine them with summary() and maybe trim them with trim().", call. = FALSE)
  call <- match.call()

  out <- list(weights = w,
              treat = treat,
              covs = reported.covs,
              s.weights = sw,
              estimand = reported.estimand,
              focal = focal,
              call = call)

  class(out) <- "optweight"
  return(out)
}

print.optweight <- function(x, ...) {
  cat("An optweight object\n")
  cat(paste0(" - number of obs.: ", length(x[["weights"]]), "\n"))
  if (is_not_null(x[["estimand"]])) cat(paste0(" - estimand: ", x[["estimand"]], ifelse(is_not_null(x[["focal"]]), paste0(" (focal: ", x[["focal"]], ")"), ""), "\n"))
  cat(paste0(" - sampling weights: ", ifelse(all_the_same(x[["s.weights"]]),"none", "present"), "\n"))
}
