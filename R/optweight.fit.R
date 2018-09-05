optweight.fit <- function(treat, covs, tols = .001, estimand = "ATE", targets = NULL, s.weights = NULL, focal = NULL, std.binary = FALSE, std.cont = TRUE, verbose = FALSE, ...) {
  #treat, covs, tols can be lists (for different times), or vec, mat, vec (respectively)
  args <- list(...)
  estimand <- toupper(estimand)
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

  if (length(times) > 1 && is_not_null(estimand) && estimand %nin% "ATE") stop("Only the ATE or specified targets are compatible with longitduinal treatments.", call. = FALSE)

  #if (is_null(estimand) || estimand %in% c("ATE", "ATT") || length(times) > 1) {
    n <- lapply(times, function(i) {
      if (treat.types[i] == "cat") vapply(split(sw, t.list[[i]]), sum, numeric(1)) #faster than tapply
      else c(cont.treat = sum(sw))
    })
    unique.treats <- lapply(n, names)

    means <- lapply(covs.list, col.w.m, w = sw)
    if (is_not_null(estimand)) {
      if (estimand %in% c("ATT", "ATC")) {
        targets <- lapply(times, function(i) {
          if (is_null(focal)) focal <- max(t.list[[i]])
          else if (estimand == "ATC") focal <- min(t.list[[i]])
          col.w.m(covs.list[[i]][t.list[[i]] == focal,], w = sw[t.list[[i]] == focal])
        })
        sds <- lapply(times, function(i) {
          if (is_null(focal)) focal <- max(t.list[[i]])
          else if (estimand == "ATC") focal <- min(t.list[[i]])
          sqrt(col.w.v(covs.list[[i]][t.list[[i]] == focal,], w = sw[t.list[[i]] == focal]))
        })
      }
      else if (estimand == "ATE") {
        targets <- means
        sds <- lapply(times, function(i) sqrt(rowMeans(matrix(sapply(unique.treats[[i]], function(t) col.w.v(covs.list[[i]][t.list[[i]]==t, , drop = FALSE], w = sw[t.list[[i]] == t]), simplify = "array"), ncol = length(unique.treats[[i]])))))
      }
    }
    else {
      if (is_null(targets)) targets <- lapply(covs.list, function(c) rep(NA_real_, ncol(c)))
      else if (is.atomic(targets)) targets <- list(targets)
      else if (!is.list(targets)) stop("targets must be a list of target values for each covariate.", call. = FALSE)

      for (i in seq_along(covs.list)) {
        if (length(targets[[i]]) != ncol(covs.list[[i]])) {
          stop("targets must have the same number of values as there are covariates.", call. = FALSE)
        }
      }
      sds <- lapply(times, function(i) sqrt(rowMeans(matrix(sapply(unique.treats[[i]], function(t) col.w.v(covs.list[[i]][t.list[[i]]==t, , drop = FALSE], w = sw[t.list[[i]] == t]), simplify = "array"), ncol = length(unique.treats[[i]])))))

    }

    targeted <- untargeted <- treat.sds <- treat.means <- tols <- vector("list", length(times))
    for (i in times) {
      if (treat.types[i] == "cat") {
        targeted[[i]] <- !is.na(targets[[i]])
        untargeted[[i]] <- !targeted[[i]]
        treat.sds[[i]] <- NA_real_
        treat.means[[i]] <- NA_real_
        #tols
        if (std.binary && std.cont) vars.to.standardize <- !check_if_zero(tols.list[[i]])
        else if (!std.binary && std.cont) vars.to.standardize <- !check_if_zero(tols.list[[i]]) & !apply(covs.list[[i]], 2, is_binary)
        else if (std.binary && !std.cont) vars.to.standardize <- !check_if_zero(tols.list[[i]]) & apply(covs.list[[i]], 2, is_binary)
        else vars.to.standardize <- rep(FALSE, length(tols.list[[i]]))

        tols[[i]] <- ifelse(vars.to.standardize,
                            abs(tols.list[[i]]*sds[[i]]), #standardize
                            abs(tols.list[[i]]))

      }
      else {
        targeted[[i]] <- !is.na(targets[[i]])
        untargeted[[i]] <- !targeted[[i]]
        #targeted.covs.list[[i]] <- sweep(covs.list[[i]][!is.na(targets[[i]])], 2, targets[[i]][!is.na(targets[[i]])], "-") #center covs at targets (which will be eventual means)
        #untargeted.covs.list[[i]] <- sweep(covs.list[[i]][is.na(targets[[i]])], 2, means[[i]][is.na(targets[[i]])], "-") #center covs at means
        covs.list[[i]][, targeted[[i]]] <- sweep(covs.list[[i]][, targeted[[i]], drop = FALSE], 2, targets[[i]][targeted[[i]]], "-") #center covs at targets (which will be eventual means)
        covs.list[[i]][, untargeted[[i]]] <- sweep(covs.list[[i]][, untargeted[[i]], drop = FALSE], 2, means[[i]][untargeted[[i]]], "-") #center covs at means
        sds[[i]] <- sqrt(col.w.v(covs.list[[i]], w = sw))
        treat.sds[[i]] <- sqrt(col.w.v(matrix(t.list[[i]], ncol = 1), w = sw))
        treat.means[[i]] <- col.w.m(matrix(t.list[[i]], ncol = 1), w = sw)
        t.list[[i]] <- t.list[[i]] - treat.means[[i]] #center treat
        tols[[i]] <- abs(tols.list[[i]]*sds[[i]]*treat.sds[[i]])
      }
    }
  # }
  # else {
  #   treat <- t.list[[1]]
  #   covs <- covs.list[[1]]
  #
  #   if (estimand == "ATT" && is_null(focal)) focal <- max(treat)
  #   else if (estimand == "ATC") focal <- min(treat)
  #
  #   n <- vapply(split(sw[treat != focal], treat[treat != focal]), sum, numeric(1))
  #   unique.treats <- names(n)
  #
  #   targets <- col.w.m(covs[treat == focal, , drop = FALSE], w = sw[treat == focal])
  #   sds <- sqrt(col.w.v(covs[treat == focal, , drop = FALSE], w = sw[treat == focal]))
  #
  #   if (std.binary && std.cont) vars.to.standardize <- !check_if_zero(tols.list[[1]])
  #   else if (!std.binary && std.cont) vars.to.standardize <- !check_if_zero(tols.list[[1]]) & !apply(covs, 2, is_binary)
  #   else if (std.binary && !std.cont) vars.to.standardize <- !check_if_zero(tols.list[[1]]) & apply(covs, 2, is_binary)
  #   else vars.to.standardize <- rep(FALSE, length(tols.list[[1]]))
  #
  #   tols <- ifelse(vars.to.standardize,
  #                  abs(tols.list[[1]]*sds), #standardize
  #                  abs(tols.list[[1]]))
  #   #tols <- ifelse(check_if_zero(tols.list[[1]]) | apply(covs, 2, function(c) !std.binary && is_binary(c)), abs(tols.list[[1]]), abs(tols.list[[1]]*sds))
  #
  #   w[treat == focal] <- 1
  #
  #   targeted <- list(!is.na(targets))
  #   untargeted <- list(is.na(targets))
  #   N <- sum(treat != focal)
  #   t.list <- list(treat[treat != focal])
  #   covs.list <- list(covs[treat != focal, , drop = FALSE])
  #   unique.treats <- list(unique.treats)
  #   tols <- list(tols)
  #   targets <- list(targets)
  #   sds <- list(sds)
  #   n <- list(n)
  #   sw <- sw[treat != focal]
  # }

  #Minimizing variance of weights
  P = sparseMatrix(1:N, 1:N, x = 2*(sw^2)/sum(sw))
  q = rep(-1/sum(sw), N) #ensures objective function value is variance of weights

  #Minimizing the sum of the variances in each treatment group
  #Note: equiv. to setting targets closer to smaller group
  # P = sparseMatrix(1:N, 1:N, x = (2*sw^2)/ifelse(t.list[[1]]==1, n[[1]]["1"], n[[1]]["0"]))
  # q = -1/ifelse(t.list[[1]]==1, n[[1]]["1"], n[[1]]["0"]) #ensures objective function value is variance of weights

  #Mean of weights in each treat must equal 1
  E1 = do.call("rbind", lapply(times, function(i) {
    if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t) (t.list[[i]] == t) * sw / n[[i]][t]))
    else sw/n[[i]]
  }))
  F1l = do.call("c", lapply(times, function(i) rep(1, length(unique.treats[[i]]))))
  F1u = F1l

  #All weights must be >= 0; focal weights must be 1, weights where sw = 0 must be 0
  G1 = sparseMatrix(1:N, 1:N, x = 1)
  if (is_not_null(focal)) {
    H1l <- ifelse(check_if_zero(sw), 0, ifelse(t.list[[1]] == focal, 1, 0))
    H1u <- ifelse(check_if_zero(sw), 0, ifelse(t.list[[1]] == focal, 1, Inf))
  }
  else {
    H1l <- rep(0, N)
    H1u <- ifelse(check_if_zero(sw), 0, Inf)
  }

  #Targeting constraints
  G2 = do.call("rbind", lapply(times, function(i) {
    if (any(targeted[[i]])) {
      if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t)
        if (is_null(focal) || (is_not_null(focal) && t != focal)) t(covs.list[[i]][, targeted[[i]], drop = FALSE] * (t.list[[i]] == t) * sw / n[[i]][t])
      ))
      else rbind(t(covs.list[[i]][, targeted[[i]], drop = FALSE] * sw), t.list[[i]] * sw) #variables are centered
    }
    else NULL
  }))
  H2l = do.call("c", lapply(times, function(i) {
    if (any(targeted[[i]])) {
      if (treat.types[i] == "cat") do.call("c", lapply(unique.treats[[i]], function(t) {
        if (is_null(focal)) targets[[i]][targeted[[i]]] - tols[[i]][targeted[[i]]]/2
        else if (is_not_null(focal) && t != focal) targets[[i]][targeted[[i]]] - tols[[i]][targeted[[i]]]
      }))
      else rep(0, sum(targeted[[i]]) + 1) #variables are centered at targets; +1 for treat
    }
    else NULL
  }))
  H2u = do.call("c", lapply(times, function(i) {
    if (any(targeted[[i]])) {
      if (treat.types[i] == "cat") do.call("c", lapply(unique.treats[[i]], function(t) {
        if (is_null(focal)) targets[[i]][targeted[[i]]] + tols[[i]][targeted[[i]]]/2
        else if (is_not_null(focal) && t != focal) targets[[i]][targeted[[i]]] + tols[[i]][targeted[[i]]]
      }))
      else rep(0, sum(targeted[[i]]) + 1) #variables are centered at targets; +1 for treat
    }
    else NULL
  }))

  #Balancing constraints for all covariates
  G3 <- do.call("rbind", lapply(times, function(i) {
    if (treat.types[i] == "cat") do.call("rbind", lapply(combn(unique.treats[[i]], 2, simplify = FALSE), function(comb) {
      t(covs.list[[i]] * (t.list[[i]] == comb[1]) * sw / n[[i]][comb[1]]) - t(covs.list[[i]] * (t.list[[i]] == comb[2]) * sw / n[[i]][comb[2]])
    }))
    else t(covs.list[[i]] * t.list[[i]] * sw / n[[i]]) #For cont, all have balancing constraints
  }))
  H3l <- do.call("c", lapply(times, function(i) {
    if (treat.types[i] == "cat") rep(-tols[[i]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
    else -tols[[i]]
  }))
  H3u <- -H3l #(positive tols)

  #Process args
  args[names(args) %nin% names(formals(rosqp::osqpSettings))] <- NULL
  if (is_null(args[["adaptive_rho"]])) args[["adaptive_rho"]] <- TRUE
  if (is_null(args[["max_iter"]])) args[["max_iter"]] <- 2E5
  if (is_null(args[["eps_abs"]])) args[["eps_abs"]] <- 1E-9
  if (is_null(args[["eps_rel"]])) args[["eps_rel"]] <- 1E-9
  args[["verbose"]] <- verbose

  A  <- rbind(G1, E1, G3, G2)
  lower <- c(H1l, F1l, H3l, H2l)
  upper <- c(H1u, F1u, H3u, H2u)

  out <- rosqp::solve_osqp(P = P, q = q, A = A, l = lower, u = upper,
                           pars = do.call(rosqp::osqpSettings, args))

  #Get dual vars for constraints
  # balance_duals <- out$y[-seq_len(length(out$y)-nrow(G2))]
  duals <- vector("list", length(times))

  # k <- 1
  # for (i in times) {
  #   ncovs <- ncol(covs.list[[i]])
  #   ntreats <- length(unique.treats[[i]])
  #   duals[[i]] <- as.data.frame(matrix(abs(balance_duals[k:(k+ncovs*ntreats-1)]/out$info$obj_val),
  #                                      byrow = FALSE, ncol = ntreats, nrow = ncovs,
  #                                      dimnames = list(colnames(covs.list[[i]]),
  #                                                      unique.treats[[i]])))
  #   k <- k + ncovs*ntreats
  # }

  w_ <- out$x
  w_[w_ < 0] <- 0

  if (is_null(focal)) w <- w_
  else w[w != 1] <- w_

  opt_out <- list(w = w,
                  duals = duals,
                  info = out$info)
  return(opt_out)
}
