optweight.fit <- function(treat, covs, tols = .001, estimand = "ATE", s.weights = NULL, focal = NULL, std.binary = FALSE, exact = FALSE, verbose = FALSE, ...) {
  #treat, covs, tols can be lists (for different times), or vec, mat, vec (respectively)
  args <- list(...)
  if (!missing(treat)) t.list <- treat
  if (!missing(covs)) covs.list <- covs
  if (is.atomic(t.list)) t.list <- list(t.list)
  if (is.matrix(covs.list) || is.data.frame(covs.list)) covs.list <- list(covs.list)
  t.list <- lapply(t.list, as.character)
  if (!all(vapply(covs.list, function(c) all(apply(c, 2, is.numeric)), logical(1L)))) stop("All covariates must be numeric.", call. = FALSE)
  covs.list <- lapply(covs.list, as.matrix)

  times <- seq_along(covs.list)

  if (!exact) {
    if (!missing(tols)) tols.list <- tols
    if (is.atomic(tols.list)) tols.list <- list(tols.list)
    if (length(tols.list) == 1) tols.list <- replicate(max(times), tols.list[[1]], simplify = FALSE)
    tols.list <- lapply(times, function(i) if (length(tols.list[[i]] == 1)) rep(tols.list[[i]], ncol(covs.list[[i]])) else tols.list[[i]])
  }

  N <- nrow(covs.list[[1]])
  if (is_null(s.weights)) sw <- rep(1, N)
  else sw <- s.weights
  w <- numeric(N)

  if (estimand == "ATE") {
    n <- lapply(times, function(i) vapply(split(sw, t.list[[i]]), sum, numeric(1)))
    unique.treats <- lapply(n, names)

    means <- lapply(covs.list, col.w.m, w = sw)
    sds <- lapply(times, function(i) {
      sqrt(rowMeans(matrix(sapply(unique.treats[[i]], function(t) col.w.v(covs.list[[i]][t.list[[i]]==t, , drop = FALSE], w = sw[t.list[[i]] == t]), simplify = "array"), ncol = length(unique.treats[[i]]))))
    })
    if (!exact) {
      tols <- lapply(times, function(i) {
        ifelse(apply(covs.list[[i]], 2, function(c) !std.binary && is_binary(c)), tols.list[[i]]/2, tols.list[[i]]*sds[[i]]/2)
      })
    }

  }
  else {
    if (length(t.list) > 1) stop("Only the ATE is compatible with longitduinal treatments.", call. = FALSE)
    treat <- t.list[[1]]
    covs <- covs.list[[1]]

    if (estimand == "ATT" && is_null(focal)) focal <- max(treat)
    else if (estimand == "ATC") focal <- min(treat)

    n <- vapply(split(sw[treat != focal], treat[treat != focal]), sum, numeric(1))
    unique.treats <- names(n)

    means <- col.w.m(covs[treat == focal, , drop = FALSE], w = sw[treat == focal])
    sds <- sqrt(col.w.v(covs[treat == focal, , drop = FALSE], w = sw[treat == focal]))

    if (!exact) {
      tols <-  ifelse(apply(covs, 2, function(c) !std.binary && is_binary(c)), tols.list[[1]], tols.list[[1]]*sds)
    }
    w[treat == focal] <- 1

    N <- sum(treat != focal)
    t.list <- list(treat[treat != focal])
    covs.list <- list(covs[treat != focal, , drop = FALSE])
    unique.treats <- list(unique.treats)
    tols <- list(tols)
    means <- list(means)
    sds <- list(sds)
    n <- list(n)
    sw <- sw[treat != focal]
  }

  #Minimizing squared distances of weights from the mean (1)
  A = if (is_null(s.weights)) diag(N) else diag(sw)
  B = rep(1, N)
  if (!exact) {
    #Sum of weights in each treat must equal size of group
    E1 = do.call("rbind", lapply(times, function(i) do.call("rbind", lapply(unique.treats[[i]], function(t) (t.list[[i]] == t) * sw))))
    F1l = do.call("c", lapply(times, function(i) n[[i]][unique.treats[[i]]]))
    F1u = F1l

    #All weights must be >= 0
    G1 = diag(N)
    H1l = rep(0, N)
    H1u = rep(Inf, N)

    #Balance constraints
    G2 = do.call("rbind", lapply(times, function(i) {
      do.call("rbind", lapply(unique.treats[[i]], function(t)
        t(covs.list[[i]] * (t.list[[i]] == t) * sw)
      ))}))
    H2l = do.call("c", lapply(times, function(i) {
      vapply(n[[i]][unique.treats[[i]]],
             function(n_) n_*(means[[i]]-tols[[i]]),
             numeric(length(means[[i]])))
    }))
    H2u = do.call("c", lapply(times, function(i) {
      vapply(n[[i]][unique.treats[[i]]],
             function(n_) n_*(means[[i]]+tols[[i]]),
             numeric(length(means[[i]])))
    }))

    E = rbind(E1)
    G = rbind(G1, G2)

    Fl = c(F1l)
    Fu = c(F1u)
    Hl = c(H1l, H2l)
    Hu = c(H1u, H2u)

  }
  else {
    #Sum of weights in each treat must equal size of group
    E1 = do.call("rbind", lapply(times, function(i) do.call("rbind", lapply(unique.treats[[i]], function(t) (t.list[[i]] == t) * sw))))
    F1l = do.call("c", lapply(times, function(i) n[[i]][unique.treats[[i]]]))
    F1u = F1l

    #Balance constraints
    E2 = do.call("rbind", lapply(times, function(i) {
      do.call("rbind", lapply(unique.treats[[i]], function(t)
        t(covs.list[[i]] * (t.list[[i]] == t) * sw)
      ))}))

    F2l = do.call("c", lapply(times, function(i) {
      vapply(n[[i]][unique.treats[[i]]],
             function(n_) n_*(means[[i]]),
             numeric(length(means[[i]])))
    }))
    F2u = F2l

    #All weights must be >= 0
    G1 = diag(N)
    H1l = rep(0, N)
    H1u = rep(Inf, N)

    E = rbind(E1, E2)
    G = rbind(G1)

    Fl = c(F1l, F2l)
    Fu = c(F1u, F2u)
    Hl = c(H1l)
    Hu = c(H1u)

  }

  #Process args
  args[names(args) %nin% names(formals(rosqp::osqpSettings))] <- NULL
  if (is_null(args[["adaptive_rho"]])) args[["adaptive_rho"]] <- FALSE
  if (is_null(args[["max_iter"]])) args[["max_iter"]] <- 2E5
  args[["verbose"]] <- verbose

  A <- A^2
  Amat  <- rbind(E,G)
  lower  <- c(Fl, Hl)
  upper <- c(Fu, Hu)

  out <- rosqp::solve_osqp(P = A, q = -B, A = Amat, l = lower, u = upper,
                           pars = do.call(rosqp::osqpSettings, args))
  w_ <- out$x

  if (is_null(focal)) w <- w_
  else w[w != 1] <- w_

  return(w)
}
