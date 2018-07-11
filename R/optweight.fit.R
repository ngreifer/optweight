optweight.fit <- function(treat, covs, tols = .0001, estimand = "ATE", focal = NULL, s.weights = NULL, std.binary = FALSE) {
 #treat, covs, tols can be lists (for different times), or vec, mat, vec (respectively)
  if (!missing(treat)) t.list <- treat
  if (!missing(covs)) covs.list <- covs
  if (!missing(tols)) tols.list <- tols
  if (is.atomic(t.list)) t.list <- list(t.list)
  if (is.matrix(covs.list) || is.data.frame(covs.list)) covs.list <- list(covs.list)
  t.list <- lapply(t.list, as.character)
  covs.list <- lapply(covs.list, as.matrix)
  times <- seq_along(covs.list)
  if (is.atomic(tols.list)) tols.list <- list(tols.list)
  if (length(tols.list) == 1) tols.list <- replicate(max(times), tols.list[[1]], simplify = FALSE)
  tols.list <- lapply(times, function(i) if (length(tols.list[[i]] == 1)) rep(tols.list[[i]], ncol(covs.list[[i]])) else tols.list[[i]])

  unique.treats <- lapply(t.list, unique)
  N <- nrow(covs.list[[1]])
  if (is_null(s.weights)) sw <- rep(1, N)
  else sw <- s.weights
  n <- lapply(times, function(i) setNames(sapply(unique.treats[[i]], function(t) sum(sw[t.list[[i]] == t])),
                                          unique.treats[[i]]))

  if (estimand == "ATE") {
    means <- lapply(covs.list, col.w.m, w = sw)
    sds <- lapply(times, function(i) {
      sqrt(rowMeans(sapply(unique.treats[[i]], function(t) col.w.v(covs.list[[i]][t.list[[i]]==t,], w = sw[t.list[[i]] == t]))))
    })
    tols <- lapply(times, function(i) {
      ifelse(apply(covs.list[[i]], 2, function(c) !std.binary && is_binary(c)), tols.list[[i]]/2, tols.list[[i]]*sds[[i]]/2)
    })
  }
  else {
    if (estimand == "ATT" && is_null(focal)) focal <- sapply(t.list, max)
    else if (estimand == "ATC") focal <- sapply(t.list, min)

    means <- lapply(times, function(i) col.w.m(covs.list[[i]][t.list[[i]] == focal[i], , drop = FALSE], w = sw[t.list[[i]] == focal[i]]))
    sds <- lapply(times, function(i) {
      sqrt(col.w.v(covs.list[[i]][t.list[[i]] == focal[i], , drop = FALSE], w = sw[t.list[[i]] == focal[i]]))
    })
    tols <- lapply(times, function(i) {
      ifelse(apply(covs.list[[i]], 2, function(c) !std.binary && is_binary(c)), tols.list[[i]], tols.list[[i]]*sds[[i]])
    })
  }

  A = diag(N)%*%diag(sw)
  B = rep(1, N)
  E = do.call("rbind", lapply(times, function(i) do.call("rbind", lapply(unique.treats[[i]], function(t) t(as.numeric(t.list[[i]] == t))%*%diag(sw)))))
  F = unlist(n)
  if (is_not_null(focal) && is_not_null(s.weights)) {
    E0 <- do.call("rbind", lapply(times, function(i) diag(as.numeric(t.list[[i]] == focal[i]))))
    F0 <- do.call("c", lapply(times, function(i) as.numeric(t.list[[i]] == focal[i])))
    E <- rbind(E, E0)
    F <- c(F, F0)
  }
  G = do.call("rbind", c(lapply(times, function(i) {
    do.call("rbind", lapply(unique.treats[[i]], function(t)
      rbind(t(covs.list[[i]])%*%diag(as.numeric(t.list[[i]] == t))%*%diag(sw),
            -t(covs.list[[i]])%*%diag(as.numeric(t.list[[i]] == t))%*%diag(sw))
    ))}), list(diag(N))))
  H = do.call("c", c(lapply(times, function(i) {
    do.call("c", lapply(unique.treats[[i]], function(t) c(n[[i]][t]*(-tols[[i]]+means[[i]]),
                                                          n[[i]][t]*(-tols[[i]]-means[[i]]))))
  }), list(rep(0, N))))

  out <- tryCatch(limSolve::lsei(A = A, B = B, E = E, F = F, G = G, H = H),
                  warning = function(w) stop("There is no solution with the given constraints.", call. = FALSE))

  return(out$X)
}
