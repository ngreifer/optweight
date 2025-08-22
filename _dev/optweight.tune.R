optweight.tune <- function(formula, data = NULL, candidate.tols = list(.1, .01, .001, .0001), nboot = 100, bal.tab.control = list(), ...) {
  A <- list(...)
  ow.list <- vector("list", length(candidate.tols))
  mean.boot.imbal.list <- vector("numeric", length(candidate.tols))

  if (is_not_null(bal.tab.control)) {
    bal.tab.control <- bal.tab.control[names(bal.tab.control) %in% c("binary", "continuous", "s.d.denom", "pairwise")]
    if (is_null(bal.tab.control)) bal.tab.control <- list()
  }

  if ((is_null(A[["std.binary"]]) || !isTRUE(A[["std.binary"]])) &&
      (is_not_null(bal.tab.control[["binary"]]) && !is.na(pmatch(bal.tab.control[["binary"]], "std")))) {
        warning("std.binary is FALSE but binary is set to \"std\" in bal.tab.control. These should be aligned.", call. = FALSE, immediate. = TRUE)
  }
  if (isTRUE(A[["std.binary"]]) &&
      (is_null(bal.tab.control[["binary"]]) || !is.na(pmatch(bal.tab.control[["binary"]], "raw")))) {
    warning("std.binary is TRUE but binary is set to \"raw\" in bal.tab.control. These should be aligned.", call. = FALSE, immediate. = TRUE)
  }

  if ((is_null(A[["std.cont"]]) || isTRUE(A[["std.cont"]])) &&
      (is_not_null(bal.tab.control[["continuous"]]) && !is.na(pmatch(bal.tab.control[["continuous"]], "raw")))) {
    warning("std.cont is TRUE but continuous is set to \"raw\" in bal.tab.control. These should be aligned.", call. = FALSE, immediate. = TRUE)
  }
  if (isFALSE(A[["std.cont"]]) &&
      (is_null(bal.tab.control[["continuous"]]) || !is.na(pmatch(bal.tab.control[["continuous"]], "std")))) {
    warning("std.cont is FALSE but continuous is set to \"std\" in bal.tab.control. These should be aligned.", call. = FALSE, immediate. = TRUE)
  }

  if (is.list(formula) && length(formula) > 1) {
    for (i in seq_along(candidate.tols)) {
      ow <- optweight(formula, data, tols = candidate.tols[[i]], ...)

      if (i == 1) {
        covs.list <- ow$covs.list
        treat.list <- ow$treat.list
        N <- nrow(covs.list[[1]])
      }

      mean.boot.imbal.list[i] <- mean(
        vapply(1:nboot, function(b) {
          boot.sample.indices <- sample(1:N, N, replace = TRUE)
          ow_b <- ow
          ow_b[["covs.list"]] <- lapply(ow_b[["covs.list"]], function(c) c[boot.sample.indices, , drop = FALSE])
          ow_b[["treat.list"]] <- lapply(ow_b[["treat.list"]], function(t) t[boot.sample.indices])
          ow_b[["weights"]] <- ow_b[["weights"]][boot.sample.indices]
          ow_b[["s.weights"]] <- ow_b[["s.weights"]][boot.sample.indices]

          B.list <- do.call(cobalt::bal.tab, c(list(ow_b,
                                                    method = "weighting",
                                                    msm.summary = TRUE,
                                                    multi.summary = TRUE,
                                                    quick = TRUE),
                                            bal.tab.control))[["Time.Balance"]]
          bal <- do.call("c", lapply(B.list, function(B_) {
            if ("Balance" %in% names(B_)) {
              B <- B_[["Balance"]]
              if ("Corr.Adj" %in% names(B)) bal_ <- B[["Corr.Adj"]]
              else bal_ <- B[["Diff.Adj"]]
            }
            else if ("Balance.Across.Pairs" %in% names(B_)) {
              B <- B_[["Balance.Across.Pairs"]]
              bal_ <- B[["Max.Diff.Adj"]]
            }
            return(bal_)
          }))

          return(sqrt(sum(bal^2)))
        }, numeric(1L))
      )

      ow[["covs.list"]] <- NULL
      ow[["treat.list"]] <- NULL
      ow.list[[i]] <- ow
      rm(ow)

    }
    best.ow <- ow.list[[which.min(mean.boot.imbal.list)]]
    best.ow[["covs.list"]] <- covs.list
    best.ow[["treat.list"]] <- treat.list
  }
  else {
    for (i in seq_along(candidate.tols)) {
      ow <- optweight(formula, data, tols = candidate.tols[[i]], ...)

      if (i == 1) {
        covs <- ow$covs
        N <- nrow(covs)
      }

      mean.boot.imbal.list[i] <- mean(
        vapply(1:nboot, function(b) {
          boot.sample.indices <- sample(1:N, N, replace = TRUE)
          ow_b <- ow
          ow_b[["covs"]] <- ow_b[["covs"]][boot.sample.indices, , drop = FALSE]
          ow_b[["treat"]] <- ow_b[["treat"]][boot.sample.indices]
          ow_b[["weights"]] <- ow_b[["weights"]][boot.sample.indices]
          ow_b[["s.weights"]] <- ow_b[["s.weights"]][boot.sample.indices]

          B_ <- do.call(cobalt::bal.tab, c(list(ow_b, method = "weighting",
                                               multi.summary = TRUE,
                                               quick = TRUE),
                                            bal.tab.control))

          if ("Balance" %in% names(B_)) {
            B <- B_[["Balance"]]
            if ("Corr.Adj" %in% names(B)) bal <- B[["Corr.Adj"]]
            else bal <- B[["Diff.Adj"]]
          }
          else if ("Balance.Across.Pairs" %in% names(B_)) {
            B <- B_[["Balance.Across.Pairs"]]
            bal <- B[["Max.Diff.Un"]]
          }

          return(sqrt(mean(bal^2)))
        }, numeric(1L))
      )

      ow[["covs"]] <- NULL
      ow.list[[i]] <- ow
      rm(ow)

    }
    best.ow <- ow.list[[which.min(mean.boot.imbal.list)]]
    best.ow[["covs"]] <- covs
  }
  return(best.ow)
}
