constraint_mean_w_cat <- function(treat, unique.treats, sw, n) {
  #Mean of weights in each treat must equal 1
  out <- make_list(c("A", "L", "U"))

  out$A <- do.call("rbind", lapply(unique.treats, function(i) (treat == i) * sw / n[i]))
  out$L <- rep_with(1, unique.treats)
  out$U <- rep_with(1, unique.treats)

  out$treat <- unique.treats

  out
}
constraint_mean_w_cont <- function(sw, n) {
  #Mean of weights must equal 1
  out <- make_list(c("A", "L", "U"))

  out$A <- rbind(sw / n)
  out$L <- 1
  out$U <- 1

  out
}
constraint_mean_w_svy <- function(sw, n) {
  #Mean of weights must equal 1
  out <- make_list(c("A", "L", "U"))

  out$A <- rbind(sw / n)
  out$L <- 1
  out$U <- 1

  out
}

constraint_range_w <- function(sw, min.w, focal = NULL, treat = NULL) {
  out <- make_list(c("A", "L", "U"))

  N <- length(sw)

  out$A <- Matrix::sparseMatrix(seq_len(N), seq_len(N), x = 1)
  out$L <- rep.int(min.w, N)
  out$U <- rep.int(Inf, N)

  zero_sw <- check_if_zero(sw)
  if (any(zero_sw)) {
    out$L[zero_sw] <- 0
    out$U[zero_sw] <- 0
  }

  if (is_not_null(focal)) {
    in_focal <- which(treat == focal)
    out$L[in_focal] <- 1
    out$U[in_focal] <- 1
  }

  out
}

constraint_target_cat <- function(covs, treat, sw, targets, tols, targeted, unique.treats, n, focal) {
  #Targeting constraints
  #Note: need 2 * in order to simulate tols/2 but using original tols.
  #This makes dual variables work as expected.
  out <- make_list(c("A", "L", "U"))

  if (any(targeted)) {
    covs_targeted_sw <- covs[, targeted, drop = FALSE] * sw

    out$A <- do.call("rbind", lapply(unique.treats, function(i) {
      if (is_null(focal) || i != focal) t(covs_targeted_sw * (treat == i) / n[i])
    }))

    out$L <- unlist(lapply(unique.treats, function(i) {
      if (is_null(focal)) targets[targeted] - tols[targeted] / 2
      else if (i != focal) targets[targeted] - tols[targeted]
    }))

    out$U <- unlist(lapply(unique.treats, function(i) {
      if (is_null(focal)) targets[targeted] + tols[targeted] / 2
      else if (i != focal) targets[targeted] + tols[targeted]
    }))

    out$treat <- unlist(lapply(unique.treats, function(i) {
      if (is_null(focal) || i != focal) rep.int(i, sum(targeted))
    }))

    out$covs <- unlist(lapply(unique.treats, function(i) {
      if (is_null(focal) || i != focal) names(tols[targeted])
    }))
  }

  out
}
constraint_target_cont <- function(covs, treat, sw, targeted) {
  #Targeting constraints
  out <- make_list(c("A", "L", "U"))

  out$A <- rbind(treat * sw)
  out$L <- 0
  out$U <- 0

  out$treat <- NA
  out$covs <- NA

  if (any(targeted)) {
    out$A <- rbind(out$A, t(covs[, targeted, drop = FALSE] * sw)) #variables are centered
    out$L <- c(out$L, rep.int(0, sum(targeted))) #variables are centered at targets; +1 for treat
    out$U <- c(out$U, rep.int(0, sum(targeted)))

    out$treat <- c(out$treat, rep.int(NA, sum(targeted)))
    out$covs <- c(out$covs, colnames(covs)[targeted])
  }

  out
}
constraint_target_svy <- function(covs, sw, targets, targeted, tols) {
  #Targeting constraints
  out <- make_list(c("A", "L", "U"))

  if (any(targeted)) {
    covs_targeted_sw <- covs[, targeted, drop = FALSE] * sw / sum(sw)

    out$A <- t(covs_targeted_sw)
    out$L <- targets[targeted] - tols[targeted]
    out$U <- targets[targeted] + tols[targeted]

    out$covs <- colnames(covs)[targeted]
  }

  out
}

constraint_balance_cat <- function(covs, treat, sw, tols, balanced, unique.treats, n) {
  #Balancing constraints for all covariates
  out <- make_list(c("A", "L", "U"))

  if (any(balanced)) {
    covs_balanced_sw <- covs[, balanced, drop = FALSE] * sw

    combos <- utils::combn(unique.treats, 2L, simplify = FALSE)

    out$A <- do.call("rbind", lapply(combos, function(comb) {
      t(covs_balanced_sw * ((treat == comb[1L]) / n[comb[1L]] - (treat == comb[2L]) / n[comb[2L]]))
    }))

    out$L <- rep_with(-tols[balanced], combos)
    out$U <- rep_with(tols[balanced], combos)

    out$treat.comb <- unlist(lapply(combos, function(comb) {
      sprintf("%s vs. %s", comb[1L], comb[2L])
    }))

    out$covs <- rep_with(colnames(covs)[balanced], combos)
  }

  out
}
constraint_balance_cont <- function(covs, treat, sw, tols, balanced, corr.type) {
  #Balancing constraints for all covariates
  out <- make_list(c("A", "L", "U"))

  if (any(balanced)) {
    covs_balanced_sw <- switch(corr.type,
                               pearson = covs[, balanced, drop = FALSE] * sw,
                               spearman = apply(covs[, balanced, drop = FALSE], 2L, rank) * sw)

    n <- length(treat)

    correct.factor <- 0 * 2 #see w.cov

    out$A <- t(covs_balanced_sw * treat / (n - correct.factor))
    out$L <- -tols[balanced]
    out$U <- tols[balanced]

    out$treat <- rep.int(NA, sum(balanced))
    out$covs <- colnames(covs)[balanced]
  }

  out
}

combine_constraints <- function(x, constraint.list) {
  do.call(switch(x, A = "rbind", "c"),
          grab(constraint.list, x))
}

## L2
objective_L2 <- function(bw, sw) {
  N <- length(bw)
  #Minimizing variance of weights
  P <- Matrix::sparseMatrix(seq_len(N), seq_len(N),
                            x = 2 * sw / N)
  # q = -sw/N #ensures objective function value is variance of weights
  q <- -2 * bw * sw / N

  list(P = P, q = q)
}

## L1
objective_L1 <- function(bw, sw) {
  N <- length(bw)
  eps <- 1e-9

  P <- Matrix::sparseMatrix(seq_len(N), seq_len(N),
                            x = eps * 2 * sw / N,
                            dims = c(2 * N, 2 * N))

  q <- c(eps * -2 * bw * sw / N,
         sw / N)

  list(P = P, q = q)
}

modify_constraints_L1 <- function(constraint_df, bw) {
  N <- length(bw)

  for (i in seq_row(constraint_df)) {
    if (is_not_null(constraint_df[["constraint"]][[i]][["A"]])) {
      constraint_df[["constraint"]][[i]][["A"]] <- cbind(
        constraint_df[["constraint"]][[i]][["A"]],
        Matrix::sparseMatrix(integer(), integer(),
                             dims = c(nrow(constraint_df[["constraint"]][[i]][["A"]]), N))
      )
    }
  }

  constraint_df
}

constraint_conversion_L1 <- function(bw, sw) {
  N <- length(bw)

  out <- list(
    # Constraining slack variables to be positive
    list(
      A = Matrix::sparseMatrix(seq_len(N), N + seq_len(N),
                               x = 1,
                               dims = c(N, 2 * N)),
      L = rep.int(0, N),
      U = rep.int(Inf, N),
      treat = rep.int(NA, N),
      covs = rep.int(NA, N)
    ),

    # Bounding weights by slack variables
    list(
      A = Matrix::sparseMatrix(c(seq_len(N), seq_len(N)),
                               c(seq_len(N), N + seq_len(N)),
                               x = 1,
                               dims = c(N, 2 * N)),
      L = bw,
      U = rep.int(Inf, N),
      treat = rep.int(NA, N),
      covs = rep.int(NA, N)
    ),

    list(
      A = Matrix::sparseMatrix(c(seq_len(N), seq_len(N)),
                               c(seq_len(N), N + seq_len(N)),
                               x = c(rep.int(1, N), rep.int(-1, N)),
                               dims = c(N, 2 * N)),
      L = rep.int(-Inf, N),
      U = bw,
      treat = rep.int(NA, N),
      covs = rep.int(NA, N)
    )
  )

  lapply(names(out[[1L]]), combine_constraints, out) |>
    setNames(names(out[[1L]]))
}

## Linf
objective_Linf <- function(bw, sw) {
  N <- length(bw)
  eps <- 1e-5

  #Minimizing Linf norm of w - bw
  P <- Matrix::sparseMatrix(seq_len(N), seq_len(N),
                            x = eps * 2 * sw / N,
                            dims = c(N + 1, N + 1))

  q <- c(eps * -2 * bw * sw / N,
         1)

  list(P = P, q = q)
}

modify_constraints_Linf <- function(constraint_df, bw) {
  N <- length(bw)

  for (i in seq_row(constraint_df)) {
    if (is_not_null(constraint_df[["constraint"]][[i]][["A"]])) {
      constraint_df[["constraint"]][[i]][["A"]] <- cbind(
        constraint_df[["constraint"]][[i]][["A"]],
        Matrix::sparseMatrix(integer(), integer(),
                             dims = c(nrow(constraint_df[["constraint"]][[i]][["A"]]), 1))
      )
    }
  }

  constraint_df
}

constraint_conversion_Linf <- function(bw, sw) {
  N <- length(bw)

  out <- list(
    list(
      A = Matrix::sparseMatrix(1, N + 1,
                               x = 1,
                               dims = c(1, N + 1)),
      L = 0,
      U = Inf,
      treat = NA,
      covs = NA
    ),

    # Bounding weights by slack variables
    list(
      A = Matrix::sparseMatrix(c(seq_len(N), seq_len(N)),
                               c(seq_len(N), rep.int(N + 1L, N)),
                               x = c(sw, rep.int(1L, N)),
                               dims = c(N, N + 1L)),
      L = sw * bw,
      U = rep.int(Inf, N),
      treat = rep.int(NA, N),
      covs = rep.int(NA, N)
    ),
    list(
      A = Matrix::sparseMatrix(c(seq_len(N), seq_len(N)),
                               c(seq_len(N), rep.int(N + 1L, N)),
                               x = c(sw, rep.int(-1L, N)),
                               dims = c(N, N + 1L)),
      L = rep.int(-Inf, N),
      U = sw * bw,
      treat = rep.int(NA, N),
      covs = rep.int(NA, N)
    )
  )

  lapply(names(out[[1L]]), combine_constraints, out) |>
    setNames(names(out[[1L]]))
}
