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
constraint_balance_cont <- function(covs, treat, sw, tols, balanced, n) {
  #Balancing constraints for all covariates
  out <- make_list(c("A", "L", "U"))

  if (any(balanced)) {
    # Note: denom of (n - n / length(sw)) is necessary to accord with col_w_corr()
    out$A <- t(covs[, balanced, drop = FALSE] * treat * sw / (n - n / length(sw)))
    out$L <- -tols[balanced]
    out$U <- tols[balanced]

    out$treat <- rep.int(NA, sum(balanced))
    out$covs <- colnames(covs)[balanced]
  }

  out
}

constraint_target_cat <- function(covs, treat, sw, targets, tols, targeted, unique.treats, n, focal = NULL) {
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
constraint_target_cont <- function(covs, treat, sw, n, treat.name) {
  #Targeting constraints (all covs are targeted, all have mean equal to target)
  out <- make_list(c("A", "L", "U"))

  out$A <- rbind(treat * sw / n,
                 t(covs * sw / n))
  out$L <- c(0, rep.int(0, ncol(covs)))
  out$U <- c(0, rep.int(0, ncol(covs)))

  out$treat <- c(treat.name, rep.int(NA, ncol(covs)))
  out$covs <- c(NA, colnames(covs))

  out
}
constraint_target_svy <- function(covs, sw, targets, targeted, tols, n) {
  #Targeting constraints
  out <- make_list(c("A", "L", "U"))

  if (any(targeted)) {
    out$A <- t(covs[, targeted, drop = FALSE] * sw / n)
    out$L <- targets[targeted] - tols[targeted]
    out$U <- targets[targeted] + tols[targeted]

    out$covs <- colnames(covs)[targeted]
  }

  out
}

combine_constraints <- function(x, constraint.list) {
  do.call(switch(x, A = "rbind", "c"),
          grab(constraint.list, x) |> lapply(unname))
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

  # P <- Matrix::sparseMatrix(seq_len(N), seq_len(N),
  #                           x = eps * 2 * sw / N,
  #                           dims = c(2 * N, 2 * N))
  #
  # q <- c(eps * -2 * bw * sw / N,
  #        sw / N)

  P <- NULL
  q <- c(rep.int(0, N),
         sw / N)

  list(P = P, q = q)
}
modify_constraints_L1 <- function(constraint_df, bw) {
  N <- length(bw)

  for (i in which(constraint_df[["type"]] != "conversion")) {
    if (is_not_null(constraint_df[["constraint"]][[i]][["A"]])) {
      constraint_df[["constraint"]][[i]][["A"]] <- cbind(
        constraint_df[["constraint"]][[i]][["A"]],
        Matrix::Matrix(0, nrow = nrow(constraint_df[["constraint"]][[i]][["A"]]),
                       ncol = N, sparse = TRUE)
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
  # P <- Matrix::sparseMatrix(seq_len(N), seq_len(N),
  #                           x = eps * 2 * sw / N,
  #                           dims = c(N + 1, N + 1))
  #
  # q <- c(eps * -2 * bw * sw / N,
  #        1)

  P <- NULL
  q <- c(rep.int(0, N),
         1)

  list(P = P, q = q)
}
modify_constraints_Linf <- function(constraint_df, bw) {
  for (i in which(constraint_df[["type"]] != "conversion")) {
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

## Entropy
### These work for clarabel and scs
objective_entropy <- function(bw, sw) {
  N <- length(bw)

  P <- NULL

  q <- c(rep(0, N), sw / N)

  list(P = P, q = q)
}
modify_constraints_entropy <- function(constraint_df, bw) {
  N <- length(bw)

  for (i in which(constraint_df[["type"]] != "conversion")) {
    if (is_not_null(constraint_df[["constraint"]][[i]][["A"]])) {
      #Add exponential variables
      constraint_df[["constraint"]][[i]][["A"]] <- cbind(
        constraint_df[["constraint"]][[i]][["A"]],
        Matrix::sparseMatrix(integer(), integer(),
                             dims = c(nrow(constraint_df[["constraint"]][[i]][["A"]]), N))

      )
    }
  }

  constraint_df
}
constraint_conversion_entropy <- function(bw, sw) {
  N <- length(bw)
  ind <- seq(1, 3 * N, by = 3)
  hexp <- rep.int(0, 3 * N)
  hexp[ind + 2] <- bw

  out <- list(
    list(
      A = Matrix::sparseMatrix(c(ind, ind + 1),
                               c(N + seq_len(N), seq_len(N)),
                               x = c(rep.int(1, N), rep.int(-1, N)),
                               dims = c(3 * N, 2 * N)),
      L = rep.int(-Inf, 3 * N),
      U = hexp
    )
  )

  lapply(names(out[[1L]]), combine_constraints, out) |>
    setNames(names(out[[1L]]))
}

## log
objective_log <- function(bw, sw) {
  N <- length(bw)

  P <- NULL

  q <- c(rep(0, N), -sw / N)

  list(P = P, q = q)
}
modify_constraints_log <- function(constraint_df, bw) {
  N <- length(bw)

  for (i in which(constraint_df[["type"]] != "conversion")) {
    if (is_not_null(constraint_df[["constraint"]][[i]][["A"]])) {
      #Add exponential variables
      constraint_df[["constraint"]][[i]][["A"]] <- cbind(
        constraint_df[["constraint"]][[i]][["A"]],
        Matrix::sparseMatrix(integer(), integer(),
                             dims = c(nrow(constraint_df[["constraint"]][[i]][["A"]]), N))

      )
    }
  }

  constraint_df
}
constraint_conversion_log <- function(bw, sw) {
  N <- length(bw)
  ind <- seq(1, 3 * N, by = 3)
  hexp <- rep.int(0, 3 * N)
  hexp[ind + 1] <- 1

  out <- list(
    list(
      A = Matrix::sparseMatrix(c(ind, ind + 2),
                               c(N + seq_len(N), seq_len(N)),
                               x = -1,
                               dims = c(3 * N, 2 * N)),
      L = rep.int(-Inf, 3 * N),
      U = hexp
    )
  )

  lapply(names(out[[1L]]), combine_constraints, out) |>
    setNames(names(out[[1L]]))
}

prep_constraint_df <- function(constraint_df, norm, bw, sw) {
  if (norm == "l1") {
    constraint_df <- modify_constraints_L1(constraint_df, bw) |>
      rbind(expand.grid(time = 0,
                        type = "conversion",
                        constraint = list(constraint_conversion_L1(bw, sw)),
                        stringsAsFactors = FALSE,
                        KEEP.OUT.ATTRS = FALSE))
  }
  else if (norm == "linf") {
    constraint_df <- modify_constraints_Linf(constraint_df, bw) |>
      rbind(expand.grid(time = 0,
                        type = "conversion",
                        constraint = list(constraint_conversion_Linf(bw, sw)),
                        stringsAsFactors = FALSE,
                        KEEP.OUT.ATTRS = FALSE))
  }
  else if (norm == "entropy") {
    constraint_df <- modify_constraints_entropy(constraint_df, bw) |>
      rbind(expand.grid(time = 0,
                        type = "conversion",
                        constraint = list(constraint_conversion_entropy(bw, sw)),
                        stringsAsFactors = FALSE,
                        KEEP.OUT.ATTRS = FALSE))
  }
  else if (norm == "log") {
    constraint_df <- modify_constraints_log(constraint_df, bw) |>
      rbind(expand.grid(time = 0,
                        type = "conversion",
                        constraint = list(constraint_conversion_log(bw, sw)),
                        stringsAsFactors = FALSE,
                        KEEP.OUT.ATTRS = FALSE))
  }

  constraint_df
}

prep_constraint_df_for_solver <- function(constraint_df, solver = "osqp") {

  if (solver %in% c("scs", "clarabel")) {
    for (i in which(constraint_df[["type"]] != "conversion")) {
      # If no upper bound, set A = -A, U = - L, L = -Inf
      l_only <- constraint_df[["constraint"]][[i]][["L"]] > -Inf & constraint_df[["constraint"]][[i]][["U"]] == Inf

      if (any(l_only)) {
        constraint_df[["constraint"]][[i]][["A"]][l_only, ] <- -constraint_df[["constraint"]][[i]][["A"]][l_only, , drop = FALSE]
        constraint_df[["constraint"]][[i]][["U"]][l_only] <- -constraint_df[["constraint"]][[i]][["L"]][l_only]
        constraint_df[["constraint"]][[i]][["L"]][l_only][] <- -Inf
      }

      # If unequal box constraints, append -A with U = -L and L = -Inf
      box <- constraint_df[["constraint"]][[i]][["L"]] > -Inf & constraint_df[["constraint"]][[i]][["U"]] < Inf &
        constraint_df[["constraint"]][[i]][["L"]] != constraint_df[["constraint"]][[i]][["U"]]

      if (any(box)) {
        constraint_df[["constraint"]][[i]][["A"]] <- rbind(constraint_df[["constraint"]][[i]][["A"]],
                                                           -constraint_df[["constraint"]][[i]][["A"]][box, , drop = FALSE])
        constraint_df[["constraint"]][[i]][["U"]] <- c(constraint_df[["constraint"]][[i]][["U"]],
                                                       -constraint_df[["constraint"]][[i]][["L"]][box])

        constraint_df[["constraint"]][[i]][["L"]][box] <- -Inf
        constraint_df[["constraint"]][[i]][["L"]] <- c(constraint_df[["constraint"]][[i]][["L"]],
                                                       rep(-Inf, sum(box)))

        for (j in setdiff(names(constraint_df[["constraint"]][[i]]), c("A", "U", "L"))) {
          constraint_df[["constraint"]][[i]][[j]] <- c(constraint_df[["constraint"]][[i]][[j]],
                                                       constraint_df[["constraint"]][[i]][[j]][box])
        }
      }

      # Drop constraints with infinite upper and lower bounds
      no_bound <- constraint_df[["constraint"]][[i]][["L"]] == -Inf & constraint_df[["constraint"]][[i]][["U"]] == Inf

      if (any(no_bound)) {
        constraint_df[["constraint"]][[i]][["A"]] <- constraint_df[["constraint"]][[i]][["A"]][!no_bound, , drop = FALSE]

        for (j in setdiff(names(constraint_df[["constraint"]][[i]]), "A")) {
          constraint_df[["constraint"]][[i]][[j]] <- constraint_df[["constraint"]][[i]][[j]][!no_bound]
        }
      }
    }
  }
  else if (solver == "lpsolve") {
    for (i in which(constraint_df[["type"]] != "conversion")) {
      # If unequal box constraints, append -A with U = -L and L = -Inf
      box <- constraint_df[["constraint"]][[i]][["L"]] > -Inf & constraint_df[["constraint"]][[i]][["U"]] < Inf &
        constraint_df[["constraint"]][[i]][["L"]] != constraint_df[["constraint"]][[i]][["U"]]

      if (any(box)) {
        constraint_df[["constraint"]][[i]][["A"]] <- rbind(constraint_df[["constraint"]][[i]][["A"]][!box, , drop = FALSE],
                                                           constraint_df[["constraint"]][[i]][["A"]][box, , drop = FALSE],
                                                           constraint_df[["constraint"]][[i]][["A"]][box, , drop = FALSE])

        constraint_df[["constraint"]][[i]][["U"]] <- c(constraint_df[["constraint"]][[i]][["U"]][!box],
                                                       constraint_df[["constraint"]][[i]][["U"]][box],
                                                       rep.int(Inf, sum(box)))

        constraint_df[["constraint"]][[i]][["L"]] <- c(constraint_df[["constraint"]][[i]][["L"]][!box],
                                                       rep.int(-Inf, sum(box)),
                                                       constraint_df[["constraint"]][[i]][["L"]][box])

        for (j in setdiff(names(constraint_df[["constraint"]][[i]]), c("A", "U", "L"))) {
          constraint_df[["constraint"]][[i]][[j]] <- c(constraint_df[["constraint"]][[i]][[j]][!box],
                                                       constraint_df[["constraint"]][[i]][[j]][box],
                                                       constraint_df[["constraint"]][[i]][[j]][box])
        }
      }

      # Drop constraints with infinite upper and lower bounds
      no_bound <- constraint_df[["constraint"]][[i]][["L"]] == -Inf & constraint_df[["constraint"]][[i]][["U"]] == Inf

      if (any(no_bound)) {
        constraint_df[["constraint"]][[i]][["A"]] <- constraint_df[["constraint"]][[i]][["A"]][!no_bound, , drop = FALSE]

        for (j in setdiff(names(constraint_df[["constraint"]][[i]]), "A")) {
          constraint_df[["constraint"]][[i]][[j]] <- constraint_df[["constraint"]][[i]][[j]][!no_bound]
        }
      }
    }
  }

  constraint_df$nc <- lengths(grab(constraint_df[["constraint"]], "U"))
  constraint_df$nc_cum <- cumsum(constraint_df$nc)

  constraint_df
}

prep_objective <- function(norm, bw, sw) {
  objective_fun <- switch(norm,
                          l2 = objective_L2,
                          l1 = objective_L1,
                          linf = objective_Linf,
                          entropy = objective_entropy,
                          log = objective_log)
  objective_fun(bw, sw)
}

make_process_opt_args <- function(solver) {
  if (solver == "clarabel") {
    f <- function(..., verbose = FALSE) {
      rlang::check_installed("clarabel")

      chk::chk_flag(verbose)

      args <- ...mget(names(formals(clarabel::clarabel_control)))
      eps <- ...get("eps", 1e-9)

      args[["max_iter"]] <- as.integer(...get("max_iter", 2e5))

      args[["verbose"]] <- as.logical(verbose)

      args
    }
  }
  else if (solver == "scs") {
    f <- function(..., verbose = FALSE) {
      rlang::check_installed("scs")

      chk::chk_flag(verbose)

      args <- ...mget(names(formals(scs::scs_control)))
      eps <- ...get("eps", 1e-6)

      args[["eps_rel"]] <- ...get("eps_rel", eps)
      args[["eps_abs"]] <- ...get("eps_abs", eps)
      args[["eps_infeas"]] <- ...get("eps_infeas", eps)

      args[["max_iters"]] <- as.integer(...get("max_iters", ...get("max_iter", 2e5)))

      args[["verbose"]] <- as.integer(verbose)

      args
    }
  }
  else if (solver == "osqp") {
    f <- function(..., verbose = FALSE) {
      rlang::check_installed("osqp")

      chk::chk_flag(verbose)

      args <- ...mget(names(formals(osqp::osqpSettings)))
      eps <- ...get("eps", 1e-6)

      args[["eps_rel"]] <- ...get("eps_rel", ...get("reltol", eps))
      args[["eps_abs"]] <- ...get("eps_abs", ...get("abstol", eps))
      args[["max_iter"]] <- ...get("max_iter", ...get("maxit", 2e5))

      args[["adaptive_rho_interval"]] <- ...get("adaptive_rho_interval", 10L)
      args[["polish"]] <- ...get("polish", TRUE)

      args[["verbose"]] <- verbose

      args
    }
  }
  else if (solver == "highs") {
    f <- function(..., verbose = FALSE) {
      rlang::check_installed("highs")

      chk::chk_flag(verbose)

      args <- ...mget(highs::highs_available_solver_options()[["option"]])
      eps <- ...get("eps", 1e-5)

      args[["output_flag"]] <- ...get("output_flag", verbose)
      args[["log_to_console"]] <- ...get("log_to_console", verbose)

      args
    }
  }
  else if (solver == "lpsolve") {
    f <- function(..., verbose = FALSE) {
      rlang::check_installed("lpSolve")

      list()
    }
  }
  else {
    .err("bad `solver`")
  }

  f
}

opt_fit <- function(constraint_df, objective, args, N, solver = "osqp") {
  if (solver == "clarabel") {
    A <- combine_constraints("A", constraint_df[["constraint"]])
    ub <- combine_constraints("U", constraint_df[["constraint"]])
    lb <- combine_constraints("L", constraint_df[["constraint"]])

    types <- rep(constraint_df$type, constraint_df[["nc"]])

    eq <- ub == lb

    if (any(eq)) {
      ord <- c(which(eq), which(!eq))
      A <- A[ord, , drop = FALSE]
      ub <- ub[ord]
      lb <- lb[ord]
      types <- types[ord]
      eq <- eq[ord]
    }
    else {
      ord <- seq_along(eq)
    }

    out <- clarabel::clarabel(A = A,
                              b = ub,
                              q = objective$q,
                              P = objective$P,
                              cones = list(z = as.integer(sum(types[eq] != "conversion")),
                                           l = as.integer(sum(types[!eq] != "conversion")),
                                           ep = as.integer(N)),
                              control = do.call(clarabel::clarabel_control, args))

    par_out <- out$x

    dual_out <- rep_with(0, ub)
    dual_out[ord] <- out$z

    info_out <- out[-(1:3)]

    status_val <- info_out$status
    if (is_not_null(status_val) && chk::vld_number(status_val)) {
      if (status_val == 8) {
        .wrn(sprintf("the optimization failed to find a solution after %s iterations. The problem may be infeasible or more iterations may be required. Check the dual variables to see which constraints are likely causing this issue",
                     info_out$iterations))
      }
      else if (status_val != 2) {
        .wrn("the optimization failed to find a stable solution")
      }
    }
  }
  else if (solver == "scs") {
    A <- combine_constraints("A", constraint_df[["constraint"]])
    ub <- combine_constraints("U", constraint_df[["constraint"]])
    lb <- combine_constraints("L", constraint_df[["constraint"]])

    types <- rep(constraint_df$type, constraint_df[["nc"]])

    eq <- ub == lb

    if (any(eq)) {
      ord <- c(which(eq), which(!eq))
      A <- A[ord, , drop = FALSE]
      ub <- ub[ord]
      lb <- lb[ord]
      types <- types[ord]
      eq <- eq[ord]
    }
    else {
      ord <- seq_along(eq)
    }

    out <- scs::scs(A = A,
                    b = ub,
                    obj = objective$q,
                    P = objective$P,
                    cone = list(z = as.integer(sum(types[eq] != "conversion")),
                                l = as.integer(sum(types[!eq] != "conversion")),
                                ep = as.integer(N)),
                    control = do.call(scs::scs_control, args))

    par_out <- out$x

    dual_out <- rep_with(0, ub)
    dual_out[ord] <- out$y

    info_out <- out[["info"]]

    status_val <- info_out$status_val
    if (is_not_null(status_val) && chk::vld_number(status_val)) {
      if (status_val == 2) {
        .wrn(sprintf("the optimization failed to find a solution after %s iterations. The problem may be infeasible or more iterations may be required. Check the dual variables to see which constraints are likely causing this issue",
                     info_out$iter))
      }
      else if (status_val != 1) {
        .wrn("the optimization failed to find a stable solution")
      }
    }
  }
  else if (solver == "osqp") {
    out <- osqp::solve_osqp(P = objective$P,
                            q = objective$q,
                            A = combine_constraints("A", constraint_df[["constraint"]]),
                            l = combine_constraints("L", constraint_df[["constraint"]]),
                            u = combine_constraints("U", constraint_df[["constraint"]]),
                            pars = do.call(osqp::osqpSettings, args))
    par_out <- out$x
    dual_out <- out$y
    info_out <- out$info

    status_val <- info_out$status_val
    if (is_not_null(status_val) && chk::vld_number(status_val)) {
      if (status_val == -2) {
        .wrn(sprintf("the optimization failed to find a solution after %s iterations. The problem may be infeasible or more iterations may be required. Check the dual variables to see which constraints are likely causing this issue",
                     info_out$iter))
      }
      else if (status_val != 1) {
        .wrn("the optimization failed to find a stable solution")
      }
    }
  }
  else if (solver == "highs") {
    A <- combine_constraints("A", constraint_df[["constraint"]])
    L <- combine_constraints("L", constraint_df[["constraint"]])
    U <- combine_constraints("U", constraint_df[["constraint"]])

    # Use bounding instead of constraints when possible
    types <- rep(constraint_df$type, constraint_df[["nc"]])
    range_const <- types == "range_w"

    lower <- rep(-Inf, ncol(A))
    upper <- rep(Inf, ncol(A))

    if (any(range_const)) {
      lower[seq_len(N)] <- L[range_const]
      upper[seq_len(N)] <- U[range_const]
    }

    out <- highs::highs_solve(Q = objective$P,
                              L = objective$q,
                              lower = lower,
                              upper = upper,
                              A = A[!range_const, , drop = FALSE],
                              lhs = L[!range_const],
                              rhs = U[!range_const],
                              control = do.call(highs::highs_control, args))

    par_out <- out$primal_solution

    dual_out <- rep_with(0, U)
    dual_out[!range_const] <- out$solver_msg$row_dual

    info_out <- out[-1L]

    status_val <- info_out$status
    if (is_not_null(status_val) && chk::vld_number(status_val)) {
      if (status_val != 7) {
        .wrn("the optimization failed to find a stable solution")
      }
    }
  }
  else if (solver == "lpsolve") {
    A <- combine_constraints("A", constraint_df[["constraint"]])
    L <- combine_constraints("L", constraint_df[["constraint"]])
    U <- combine_constraints("U", constraint_df[["constraint"]])

    const <- U
    const.dir <- rep_with("==", U)

    l_only <- U == Inf
    const[l_only] <- L[l_only]
    const.dir[l_only] <- ">="

    u_only <- L == -Inf
    const.dir[u_only] <- "<="

    out <- lpSolve::lp(objective.in = objective$q,
                       const.mat = as.matrix(A),
                       const.dir = const.dir,
                       const.rhs = const,
                       compute.sens	= 1)

    par_out <- out$solution
    dual_out <- out$duals[seq_len(out$const.count)]
    info_out <- out[-c(5, 9)]

    status_val <- info_out$status
    if (is_not_null(status_val) && chk::vld_number(status_val)) {
      if (status_val != 0) {
        .wrn("the optimization failed to find a stable solution")
      }
    }
  }

  list(out = out,
       par_out = par_out,
       dual_out = dual_out,
       info_out = info_out)
}

extract_weights <- function(opt_out, N, min.w, range_cons) {
  w <- opt_out$par_out[seq_len(N)]

  # Shrink tiny weights to 0
  if (abs(min.w) < .Machine$double.eps) {
    w[abs(w) < .Machine$double.eps] <- 0
  }

  # Adjust for imprecision
  w[w < range_cons$L[seq_len(N)]] <- range_cons$L[w < range_cons$L[seq_len(N)]]
  w[w > range_cons$U[seq_len(N)]] <- range_cons$U[w > range_cons$U[seq_len(N)]]

  w
}
extract_duals <- function(constraint_df, dual_out, times = 1L) {
  #Get dual vars for balance and target constraints

  duals <- make_list(length(times))

  if (is_null(dual_out)) {
    return(duals)
  }

  for (i in times) {
    td <- bd <- NULL

    ti <- which(constraint_df[["time"]] == i & constraint_df[["type"]] == "target" & constraint_df[["nc"]] > 0)
    if (is_not_null(ti)) {
      td <- do.call("rbind", lapply(ti, function(tii) {
        cons <- constraint_df[["constraint"]][[tii]][-(1:3)]
        ind <- seq(constraint_df$nc_cum[tii] - constraint_df$nc[tii],
                   constraint_df$nc_cum[tii])[-1L]

        data.frame(constraint = "target",
                   cov = if_null_then(cons$covs, NA_character_),
                   treat = if_null_then(cons$treat, NA_character_),
                   dual = abs(dual_out[ind]),
                   stringsAsFactors = FALSE)
      }))
    }

    bi <- which(constraint_df[["time"]] == i & constraint_df[["type"]] == "balance" & constraint_df[["nc"]] > 0)

    if (is_not_null(bi)) {
      bd <- do.call("rbind", lapply(bi, function(bii) {
        cons <- constraint_df[["constraint"]][[bii]][-(1:3)]
        ind <- seq(constraint_df$nc_cum[bii] - constraint_df$nc[bii],
                   constraint_df$nc_cum[bii])[-1L]

        data.frame(constraint = "balance",
                   cov = if_null_then(cons$covs, NA_character_),
                   treat = if_null_then(cons$treat.comb, NA_character_),
                   dual = abs(dual_out[ind]),
                   stringsAsFactors = FALSE)
      }))
    }

    if (is_not_null(td) || is_not_null(bd)) {
      duals[[i]] <- rbind(td, bd)
    }
  }

  duals
}
