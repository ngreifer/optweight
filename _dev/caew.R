caew <- function(formula, data = NULL, model, targets = NULL, s.weights = NULL, verbose = FALSE, ...) {

}

caew.fit <- function(treat, covs, model, targets = NULL, s.weights = NULL, aux.covs = NULL, std.binary = FALSE, std.cont = TRUE, min.w = 1E-8, verbose = FALSE, ...) {

  #For corr.type, make sure duals process correctly
  args <- list(...)

  #Process args
  if (is_not_null(args[["eps"]])) {
    if (is_null(args[["eps_abs"]])) args[["eps_abs"]] <- args[["eps"]]
    if (is_null(args[["eps_rel"]])) args[["eps_rel"]] <- args[["eps"]]
  }
  args[names(args) %nin% names(formals(osqpSettings))] <- NULL
  if (is_null(args[["max_iter"]])) args[["max_iter"]] <- 2E5L
  if (is_null(args[["eps_abs"]])) args[["eps_abs"]] <- 1E-8
  if (is_null(args[["eps_rel"]])) args[["eps_rel"]] <- 1E-8
  args[["verbose"]] <- verbose

  key.args <- c("treat", "covs", "model")
  missing.args <- args.not.list <- setNames(rep(FALSE, length(key.args)), key.args)
  for (arg in key.args) {
    if (eval(substitute(missing(q), list(q = arg)))) {
      missing.args[arg] <- TRUE
    }
  }
  if (any(missing.args)) stop(paste(word_list(names(missing.args)[missing.args]), "must be supplied."), call. = FALSE)

  N <- length(treat)
  if (is_null(s.weights)) sw <- rep(1, N)
  else sw <- s.weights

  if (length(min.w) != 1 || !is.numeric(min.w) || min.w < 0 || min.w >= 1) stop("min.w must be a single number in the interval [0, 1).", call. = FALSE)


  if (is_null(targets)) targets <- rep(NA_real_, ncol(covs))
  else if (identical(targets, "means")) targets <- col.w.m(covs, sw)
  else if (!is.atomic(targets) || (!all(is.na(targets)) && !is.numeric(targets))) {
    stop("targets must be a vector of target values for each covariate.", call. = FALSE)
  }
  else if (length(targets) != ncol(covs)) {
    stop("targets must have the same number of values as there are covariates.", call. = FALSE)
  }

  bin.covs <- apply(covs, 2, is_binary)
  sds <- rep(NA_real_, ncol(covs))
  sds[!bin.covs] <- sqrt(col.w.v(covs[,!bin.covs, drop = FALSE], w = sw))
  sds[bin.covs] <- sqrt(col.w.v.bin(covs[,bin.covs, drop = FALSE], w = sw))

  targeted <- !is.na(targets)

  covs[, !check_if_zero(sds)] <- mat_div(covs[,!check_if_zero(sds), drop = FALSE], sds[!check_if_zero(sds)])
  targets[!check_if_zero(sds)] <- targets[!check_if_zero(sds)] / sds[!check_if_zero(sds)]

  X <- cbind(1, covs); rm(covs)
  A <- treat; rm(treat)

  if (is.function(model)) {

  }
  else if (model == "logit") {
    if (is.numeric(A) || is.logical(A)) {
      score <- t(X * (A - mean(A)))
    }
    else {
      score <- do.call("rbind", lapply(unique(A), function(a) {
        t(X * (as.numeric(A == a) - mean(A == a)))
      }))
    }
  }
  else if (model == "probit") {
    if (is.numeric(A) || is.logical(A)) {
      score <- t(X * (A - mean(A)))
    }
    else {
      score <- do.call("rbind", lapply(unique(A), function(a) {
        t(X * (as.numeric(A == a) - mean(A == a)))
      }))
    }
  }
  else if (model == "poisson") {
    score <- t(X * (A - mean(A)))
  }
  else if (model == "negbin") {
    mu <- mean(A)
    theta <- MASS::theta.ml(A, mu, weights = sw)[1]
    score <- rbind(
      t(X * (A - mu)/(1+theta*mu)),
      t(X/theta * (theta*(A - mu)/(1+theta*mu) + log(1 + theta*mu) - digamma(A + 1/theta) + digamma(1/theta)))
    )
  }
  else if (model == "linear") {
    score <- rbind(t(X * (A - mean(A))/var(A)),
                   -1 + ((A - mean(A))^2/var(A)))
  }
  else if (model == "linear2") {
    score <- rbind(t(X * (A - mean(A))/var(A)),
                   t(X * (-1 + ((A - mean(A))^2/var(A)))))
  }
}
