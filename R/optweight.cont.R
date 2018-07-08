optweight.cont <- function(t.list, covs.list, tols.list) {

  if (is.numeric(t.list)) t.list <- list(t.list)
  if (is.matrix(covs.list) || is.data.frame(covs.list)) covs.list <- list(covs.list)
  covs.list <- lapply(covs.list, as.matrix)
  times <- seq_along(covs.list)
  if (is.numeric(tols.list)) {
    if (length(tols.list) == length(covs.list)) tols.list <- lapply(times, function(i) rep(tols.list[[i]], ncol(covs.list[[i]])))
    else tols.list <- lapply(times, function(i) rep(tols.list[[1]], ncol(covs.list[[i]])))
  }

  N <- nrow(covs.list[[1]])


  covars <- lapply(covs.list, function(c) rep(0, ncol(c)))
  csds <- lapply(covs.list, apply, 2, sd)
  tsds <- lapply(t.list, sd)
  tols <- lapply(times, function(i) {
    tols.list[[i]]*csds[[i]]*tsds[[i]]
  })
  centered.covs.list <- lapply(covs.list, function(c) apply(c, 2, function(x) x - mean(x)))
  centered.t.list <- lapply(t.list, function(t) t - mean(t))

  A = diag(N)
  B = matrix(1, nrow = N, ncol = 1)
  E = rbind(matrix(1, nrow = 1, ncol = N),
            do.call("rbind", lapply(covs.list, t)))
  F = c(N, N*do.call("c", lapply(covs.list, colMeans)))
  G = do.call("rbind", c(lapply(times, function(i) {
    rbind(t(centered.covs.list[[i]])%*%diag(centered.t.list[[i]]),
          -t(centered.covs.list[[i]])%*%diag(centered.t.list[[i]]))
  }), list(diag(N))))
  H = do.call("c", c(lapply(times, function(i) {
    c(N*(-tols[[i]]+covars[[i]]),
      N*(-tols[[i]]-covars[[i]]))
  }), list(rep(0, N))))
  out <- lsei(A = A, B = B, E = E, F = F, G = G, H = H)
  w <- out$X

  return(w)
}
