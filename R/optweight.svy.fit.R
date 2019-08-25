optweight.svy.fit <- function(covs, tols = 0, targets, s.weights = NULL, norm = "l2", std.binary = FALSE, std.cont = TRUE, min.w = 1E-8, verbose = FALSE, ...) {
  args <- list(...)

  #Process args
  args[names(args) %nin% names(formals(osqp::osqpSettings))] <- NULL
  if (is_null(args[["max_iter"]])) args[["max_iter"]] <- 2E5L
  if (is_null(args[["eps_abs"]])) args[["eps_abs"]] <- 1E-8
  if (is_null(args[["eps_rel"]])) args[["eps_rel"]] <- 1E-8
  args[["verbose"]] <- verbose

  key.args <- c("covs", "targets")
  missing.args <- setNames(rep(FALSE, length(key.args)), key.args)
  for (arg in key.args) {
    if (eval(substitute(missing(q), list(q = arg))) || is_null(get(arg))) {
      missing.args[arg] <- TRUE
    }
  }
  if (any(missing.args)) stop(paste(word_list(names(missing.args)[missing.args]), "must be supplied."), call. = FALSE)

  if (!all(apply(covs, 2, is.numeric))) stop("All covariates must be numeric.", call. = FALSE)
  covs <- as.matrix(covs)

  if (length(tols) == 1) tols <- rep(tols, ncol(covs))

  N <- nrow(covs)
  if (is_null(s.weights)) sw <- rep(1, N)
  else sw <- s.weights

  norm.options <- c("l2", "l1", "linf")
  if (length(norm) != 1 || !is.character(norm) || tolower(norm) %nin% norm.options) {
    stop(paste0("norm must be ", word_list(norm.options, and.or = "or", quotes = TRUE), "."), call. = FALSE)
  }
  else norm <- tolower(norm)

  if (length(min.w) != 1 || !is.numeric(min.w) || min.w < 0 || min.w >= 1) stop("min.w must be a single number in the interval [0, 1).", call. = FALSE)

  if (!is.atomic(targets) || (!all(is.na(targets)) && !is.numeric(targets))) stop("targets must be a vector of target values for each baseline covariate.", call. = FALSE)

  if (length(targets) != ncol(covs)) {
    stop("targets must have the same number of values as there are covariates.", call. = FALSE)
  }

  sds <- sqrt(col.w.v(covs, w = sw))

  targeted <- !is.na(targets)

  #tols
  if (std.binary && std.cont) vars.to.standardize <- rep(TRUE, length(tols))
  else if (!std.binary && std.cont) vars.to.standardize <- !apply(covs, 2, is_binary)
  else if (std.binary && !std.cont) vars.to.standardize <- apply(covs, 2, is_binary)
  else vars.to.standardize <- rep(FALSE, length(tols))

  tols <- ifelse(vars.to.standardize,
                 abs(tols*sds), #standardize
                 abs(tols))
  #Note: duals work incorrecly unless tols are > 0, so replace small tols with
  #sqrt(.Machine$double.eps).
  tols <- ifelse(tols < sqrt(.Machine$double.eps),
                 sqrt(.Machine$double.eps),
                 tols)

  if (norm == "l2") {
    #Minimizing variance of weights
    P = sparseMatrix(1:N, 1:N, x = 2*(sw^2)/N)
    q = -sw/N #ensures objective function value is variance of weights

    #Mean of weights  must equal 1
    E1 = matrix(sw/N, nrow = 1)
    F1l = 1
    F1u = F1l

    #All weights must be >= min; focal weights must be 1, weights where sw = 0 must be 0
    min <- min.w
    G1 = sparseMatrix(1:N, 1:N, x = 1)
    H1l <- rep(min, N)
    H1u <- ifelse(check_if_zero(sw), min, Inf)

    #Targeting constraints
    if (any(targeted)) {
      G2 = t(covs[, targeted, drop = FALSE] * sw / N)
      H2l = targets[targeted] - tols[targeted]
      H2u = targets[targeted] + tols[targeted]
    }
    else {
      G2 <- H2l <- H2u <- NULL
    }

    A  <- rbind(G1, E1, G2)
    lower <- c(H1l, F1l, H2l)
    upper <- c(H1u, F1u, H2u)

    out <- osqp::solve_osqp(P = P, q = q, A = A, l = lower, u = upper,
                             pars = do.call(osqp::osqpSettings, args))

    #Get dual vars for balance and target constraints
    G2.indices <- if (is_null(G2)) NULL else (NROW(G1)+NROW(E1)+1):(NROW(G1)+NROW(E1)+NROW(G2))

    w <- out$x
  }
  else if (norm == "l1") {
    #Minimizing mean absolute deviation of weights
    P = sparseMatrix(NULL, NULL, dims = c(2*N, 2*N))
    q = c(rep(0, N), 2*sw/N)

    #Mean of weights must equal 1
    E1 = matrix(sw/N, nrow = 1)
    F1l = 1
    F1u = F1l

    #All weights must be >= min; focal weights must be 1, weights where sw = 0 must be 0
    #Auxilliary vars must be >= 0
    G1 = sparseMatrix(1:(2*N), 1:(2*N), x = 1)
    H1l <- rep(min.w, N)
    H1u <- ifelse(check_if_zero(sw), min.w, Inf)
    H1lz <- c(H1l, rep(0, N))
    H1uz <- c(H1u, rep(Inf, N))

    #Targeting constraints
    if (any(targeted)) {
      G2 = t(covs[, targeted, drop = FALSE] * sw / N)
      H2l = targets[targeted] - tols[targeted]
      H2u = targets[targeted] + tols[targeted]
    }
    else {
      G2 <- H2l <- H2u <- NULL
    }

    #Conversion constraints
    Inxn = sparseMatrix(1:N, 1:N, x = 1)
    I = rbind(cbind(Inxn, -Inxn),
              cbind(-Inxn, -Inxn))
    jl = rep(-Inf, 2*N)
    ju = rep(1, 2*N)

    A  <- rbind(E1, G2)
    lower <- c(F1l, H2l)
    upper <- c(F1u, H2u)

    Au <- cbind(A, matrix(0, nrow = nrow(A), ncol = N))

    Az <- rbind(Au, G1, I)
    lowerz = c(lower, H1lz, jl)
    upperz = c(upper, H1uz, ju)

    out <- osqp::solve_osqp(P = P, q = q, A = Az, l = lowerz, u = upperz,
                             pars = do.call(osqp::osqpSettings, args))

    w <- out$x[1:N]

    #Get dual vars for constraints
    G2.indices <- if (is_null(G2)) NULL else (NROW(E1)+1):(NROW(E1)+NROW(G2))

  }
  else if (norm == "linf") {
    #Minimizing largest weight
    P = sparseMatrix(NULL, NULL, dims = c(2*N, 2*N))
    q = rep(sw/N, 2)

    #Mean of weights must equal 1
    E1 = matrix(sw/N, nrow = 1)
    F1l = 1
    F1u = F1l

    #All weights must be >= min; focal weights must be 1, weights where sw = 0 must be 0
    #Auxilliary var must be >= 0
    min <- min.w
    G1 = sparseMatrix(1:(2*N), 1:(2*N), x = 1)
    H1l = rep(min, N)
    H1u = ifelse(check_if_zero(sw), min, Inf)
    H1lz = c(H1l, rep(0, N))
    H1uz = c(H1u, rep(Inf, N))

    #Targeting constraints
    if (any(targeted)) {
      G2 = t(covs[, targeted, drop = FALSE] * sw / N)
      H2l = targets[targeted] - tols[targeted]
      H2u = targets[targeted] + tols[targeted]
    }
    else {
      G2 <- H2l <- H2u <- NULL
    }

    #Conversion constraints
    Inxn = sparseMatrix(1:N, 1:N, x = 1)
    I = rbind(cbind(Inxn, -Inxn),
              cbind(-Inxn, -Inxn))
    jl = rep(-Inf, 2*N)
    ju = rep(1, 2*N)

    I2 = cbind(sparseMatrix(NULL, NULL, dims = c(N-1, N)),
               matrix(1, ncol = 1, nrow = N-1),
               sparseMatrix(1:(N-1), 1:(N-1), x = -1))
    j2l = rep(0, N-1)
    j2u = rep(0, N-1)

    I = rbind(I, I2)
    jl = c(jl, j2l)
    ju = c(ju, j2u)

    A = rbind(E1, G2)
    lower = c(F1l, H2l)
    upper = c(F1u, H2u)

    Au = cbind(A, matrix(0, nrow = NROW(A), ncol = ncol(A)))

    Az = rbind(Au, G1, I)
    lowerz = c(lower, H1lz, jl)
    upperz = c(upper, H1uz, ju)

    out <- osqp::solve_osqp(P = P, q = q, A = Az, l = lowerz, u = upperz,
                             pars = do.call(osqp::osqpSettings, args))

    w <- out$x[1:N]

    #Get dual vars for constraints
    G2.indices <- if (is_null(G2)) NULL else (NROW(E1)+1):(NROW(E1)+NROW(G2))
  }

  w[w < min.w] <- min.w

  #Duals
  target_duals <- abs(out$y[G2.indices]) #G2

  if (is_not_null(target_duals)) {

    targeted.covs <- colnames(covs)[targeted]
    if (length(targeted.covs) > 0) {
      td <- data.frame(expand.grid(constraint = "target",
                                   cov = targeted.covs,
                                   stringsAsFactors = FALSE),
                       dual = target_duals[1:length(targeted.covs)]
      )
    }
    else td <- NULL

  }
  else td <- NULL

  opt_out <- list(w = w,
                  duals = td,
                  info = out$info)
  class(opt_out) <- "optweight.svy.fit"

  return(opt_out)
}
