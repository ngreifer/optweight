optweight.fit <- function(treat, covs, tols = .001, estimand = "ATE", s.weights = NULL, focal = NULL, std.binary = FALSE, verbose = FALSE, ...) {
  #treat, covs, tols can be lists (for different times), or vec, mat, vec (respectively)
  args <- list(...)
  if (!missing(treat)) t.list <- treat
  if (!missing(covs)) covs.list <- covs
  if (is.atomic(t.list)) t.list <- list(t.list)
  if (is.matrix(covs.list) || is.data.frame(covs.list)) covs.list <- list(covs.list)
  treat.types <- vapply(t.list, function(x) {
    if (is.factor(x) || is.character(x) || is_binary(x)) "cat"
    else "cont"
  }, character(1L))
  t.list <- lapply(seq_along(treat.types), function(x) {
    if (treat.types[x] == "cat") as.character(t.list[[x]])
    else as.numeric(t.list[[x]])
  })
  if (!all(vapply(covs.list, function(c) all(apply(c, 2, is.numeric)), logical(1L)))) stop("All covariates must be numeric.", call. = FALSE)
  covs.list <- lapply(covs.list, as.matrix)

  times <- seq_along(covs.list)

  if (!missing(tols)) tols.list <- tols
  if (is.atomic(tols.list)) tols.list <- list(tols.list)
  if (length(tols.list) == 1) tols.list <- replicate(max(times), tols.list[[1]], simplify = FALSE)
  tols.list <- lapply(times, function(i) if (length(tols.list[[i]]) == 1) rep(tols.list[[i]], ncol(covs.list[[i]])) else tols.list[[i]])

  N <- nrow(covs.list[[1]])
  if (is_null(s.weights)) sw <- rep(1, N)
  else sw <- s.weights
  w <- numeric(N)

  if (length(times) > 1 && estimand != "ATE") stop("Only the ATE is compatible with longitduinal treatments.", call. = FALSE)

  if (estimand == "ATE" || length(times) > 1) {
    n <- lapply(times, function(i) {
      if (treat.types[i] == "cat") vapply(split(sw, t.list[[i]]), sum, numeric(1))
      else c(cont.treat = sum(sw))
    })
    unique.treats <- lapply(n, names)

    means <- lapply(covs.list, col.w.m, w = sw)
    covs.list <- lapply(times, function(i) {
      if (treat.types[i] == "cat") covs.list[[i]]
      else sweep(covs.list[[i]], 2, means[[i]], "-") #center covs
    })

    sds <- lapply(times, function(i) {
      if (treat.types[i] == "cat") sqrt(rowMeans(matrix(sapply(unique.treats[[i]], function(t) col.w.v(covs.list[[i]][t.list[[i]]==t, , drop = FALSE], w = sw[t.list[[i]] == t]), simplify = "array"), ncol = length(unique.treats[[i]]))))
      else sqrt(col.w.v(covs.list[[i]], w = sw))
    })
    treat.sds <- vapply(times, function(i) {
      if (treat.types[i] == "cat") NA_real_
      else sqrt(col.w.v(matrix(t.list[[i]], ncol = 1), w = sw))
    }, numeric(1L))
    treat.means <- vapply(times, function(i) {
      if (treat.types[i] == "cat") NA_real_
      else col.w.m(matrix(t.list[[i]], ncol = 1), w = sw)
    }, numeric(1L))
    t.list <- lapply(times, function(i) {
      if (treat.types[i] == "cat") t.list[[i]]
      else t.list[[i]] - treat.means[i] #center treat
    })

    tols <- lapply(times, function(i) {
      if (treat.types[i] == "cat") {
        if (std.binary) vars.to.standardize <- !check_if_zero(tols.list[[i]])
        else vars.to.standardize <- !check_if_zero(tols.list[[i]]) & !apply(covs.list[[i]], 2, is_binary)

        ifelse(vars.to.standardize,
               abs(tols.list[[i]]*sds[[i]]/2), #standardize
               abs(tols.list[[i]]/2))
      }
      else abs(tols.list[[i]]*sds[[i]]*treat.sds[i]) #Maybe make binary covs act differently
    })

  }
  else {
    treat <- t.list[[1]]
    covs <- covs.list[[1]]

    if (estimand == "ATT" && is_null(focal)) focal <- max(treat)
    else if (estimand == "ATC") focal <- min(treat)

    n <- vapply(split(sw[treat != focal], treat[treat != focal]), sum, numeric(1))
    unique.treats <- names(n)

    means <- col.w.m(covs[treat == focal, , drop = FALSE], w = sw[treat == focal])
    sds <- sqrt(col.w.v(covs[treat == focal, , drop = FALSE], w = sw[treat == focal]))

    if (std.binary) vars.to.standardize <- !check_if_zero(tols.list[[1]])
    else vars.to.standardize <- !check_if_zero(tols.list[[1]]) & !apply(covs, 2, is_binary)

    tols <- ifelse(vars.to.standardize,
           abs(tols.list[[1]]*sds), #standardize
           abs(tols.list[[1]]))
    #tols <- ifelse(check_if_zero(tols.list[[1]]) | apply(covs, 2, function(c) !std.binary && is_binary(c)), abs(tols.list[[1]]), abs(tols.list[[1]]*sds))

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
  P = if (is_null(s.weights)) sparseMatrix(1:N, 1:N, x = 1) else sparseMatrix(1:N, 1:N, x = sw^2)
  q = rep(-1, N)

  #Sum of weights in each treat must equal size of group
  E1 = do.call("rbind", lapply(times, function(i) {
    if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t) (t.list[[i]] == t) * sw))
    else sw
  }))
  F1l = do.call("c", lapply(times, function(i) n[[i]][unique.treats[[i]]]))
  F1u = F1l

  #All weights must be >= 0
  G1 = sparseMatrix(1:N, 1:N, x = 1)
  H1l = rep(0, N)
  H1u = rep(Inf, N)

  #Balance constraints
  G2 = do.call("rbind", lapply(times, function(i) {
    if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t)
      t(covs.list[[i]] * (t.list[[i]] == t) * sw / n[[i]][t])
    ))
    else t(covs.list[[i]] * t.list[[i]] * sw / n[[i]])
  }))
  H2l = do.call("c", lapply(times, function(i) {
    if (treat.types[i] == "cat") rep(means[[i]]-tols[[i]], length(unique.treats[[i]]))
    else -tols[[i]]
  }))
  H2u = do.call("c", lapply(times, function(i) {
    if (treat.types[i] == "cat") rep(means[[i]]+tols[[i]], length(unique.treats[[i]]))
    else tols[[i]]
  }))

  #For continuous treatments, constrain cov and treat means to be same as in sample
  G3 = do.call("rbind", lapply(times[treat.types == "cont"], function(i) {
    rbind(t(covs.list[[i]]), t.list[[i]])
  }))
  H3l = do.call("c", lapply(times[treat.types == "cont"], function(i) {
    rep(0, ncol(covs.list[[i]]) + 1) #variables are centered; +1 for treat
  }))
  H3u = H3l

  #Process args
  args[names(args) %nin% names(formals(rosqp::osqpSettings))] <- NULL
  if (is_null(args[["adaptive_rho"]])) args[["adaptive_rho"]] <- TRUE
  if (is_null(args[["max_iter"]])) args[["max_iter"]] <- 2E5
  if (is_null(args[["eps_abs"]])) args[["eps_abs"]] <- 1E-8
  if (is_null(args[["eps_rel"]])) args[["eps_rel"]] <- 1E-8
  args[["verbose"]] <- verbose

  A  <- rbind(G1, E1, G3, G2)
  lower <- c(H1l, F1l, H3l, H2l)
  upper <- c(H1u, F1u, H3u, H2u)

  out <- rosqp::solve_osqp(P = P, q = q, A = A, l = lower, u = upper,
                           pars = do.call(rosqp::osqpSettings, args))

  #Check for convergence
  if (out$info$status_val == -2) {
    warning(paste("The optimization failed to find a solution after", out$info$iter, "iterations. The problem may be infeasible or more interations may be required."), call. = FALSE)
  }
  else if (out$info$status_val != 1) {
    warning("The optimization failed to find a stable solution.", call. = FALSE)
  }

  #Get dual vars for constraints
  balance_duals <- out$y[-seq_len(length(out$y)-nrow(G2))]
  duals <- vector("list", length(times))

  k <- 1
  for (i in times) {
    ncovs <- ncol(covs.list[[i]])
    ntreats <- length(unique.treats[[i]])
    duals[[i]] <- as.data.frame(matrix(abs(balance_duals[k:(k+ncovs*ntreats-1)]/out$info$obj_val),
                         byrow = FALSE, ncol = ntreats, nrow = ncovs,
                         dimnames = list(colnames(covs.list[[i]]),
                                         unique.treats[[i]])))
    k <- k + ncovs*ntreats
  }

  w_ <- out$x
  w_[w_ < 0] <- 0

  if (is_null(focal)) w <- w_
  else w[w != 1] <- w_

  opt_out <- list(w = w,
                  duals = duals,
                  info = out$info)
  return(opt_out)
}
