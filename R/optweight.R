optweight <- function(formula, data, tols = .0001, estimand = "ATE", s.weights = NULL, focal = NULL, pairwise = TRUE, full.output = TRUE) {

  #Process treat and covs from formula and data
  t.c <- get.covs.and.treat.from.formula(formula, data)
  reported.covs <- t.c[["reported.covs"]]
  covs <- t.c[["model.covs"]]
  treat <- t.c[["treat"]]
  treat.name <- t.c[["treat.name"]]

  if (is_null(covs)) stop("No covariates were specified.", call. = FALSE)
  if (is_null(treat)) stop("No treatment variable was specified.", call. = FALSE)

  n <- length(treat)

  if (any(is.na(reported.covs)) || nrow(reported.covs) != n) {
    stop("No missing values are allowed in the covariates.", call. = FALSE)
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
  s.weights <- process.s.weights(s.weights, data)
  if (is_null(s.weights)) sw <- rep(1, n)
  else sw <- s.weights

  ###Run optweight.fit
  w <- optweight.fit(t.list = treat, covs.list = covs,
                       tols.list = tols, estimand = estimand,
                       focal = focal,
                       s.weights = sw, std.binary = FALSE)

  warn <- FALSE
  test.w <- w*sw
  if (treat.type == "continuous") {if (sd(test.w)/mean(test.w) > 4) warn <- TRUE}
  else if (any(sapply(unique(treat), function(x) sd(test.w[treat == x])/mean(test.w[treat == x]) > 4))) warn <- TRUE
  if (warn) warning("Some extreme weights were generated. Examine them with summary() and maybe trim them with trim().", call. = FALSE)

  if (full.output) {
    out <- list(treat = treat,
                covs = reported.covs,
                weights = w,
                s.weights = if (is_null(s.weights)) NULL else sw,
                estimand = reported.estimand,
                focal = focal)
  }
  else {
    out <- w
  }
  return(out)
}
