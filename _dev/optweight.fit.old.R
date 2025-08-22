optweight.fit.old <- function(treat.list, covs.list, tols, estimand = "ATE", targets = NULL,
                          s.weights = NULL, b.weights = NULL, focal = NULL, norm = "l2",
                          std.binary = FALSE, std.cont = TRUE, min.w = 1e-8, verbose = FALSE,
                          force = FALSE, ...) {
  #For corr.type, make sure duals process correctly
  args <- list(...)

  #Process args
  # corr.type <- if (is_not_null(args[["corr.type"]])) match_arg(tolower(args[["corr.type"]]), c("pearson", "spearman", "both")) else "pearson"
  corr.type <- "pearson"

  if (is_not_null(args[["eps"]])) {
    if (is_null(args[["eps_abs"]])) args[["eps_abs"]] <- args[["eps"]]
    if (is_null(args[["eps_rel"]])) args[["eps_rel"]] <- args[["eps"]]
  }
  args[names(args) %nin% names(formals(osqp::osqpSettings))] <- NULL
  if (is_null(args[["max_iter"]])) args[["max_iter"]] <- 2E5L
  if (is_null(args[["eps_abs"]])) args[["eps_abs"]] <- 1E-8
  if (is_null(args[["eps_rel"]])) args[["eps_rel"]] <- 1E-8
  if (is_null(args[["adaptive_rho_interval"]])) args[["adaptive_rho_interval"]] <- 10L
  args[["verbose"]] <- verbose

  key.args <- c("treat.list", "covs.list", "tols")
  missing.args <- args.not.list <- rep_with(FALSE, key.args)
  for (arg in key.args) {
    if (eval(substitute(missing(q), list(q = arg)))) {
      missing.args[arg] <- TRUE
    }
    else if (!is.vector(get(arg), mode = "list")) {
      args.not.list[arg] <- TRUE
    }
  }

  if (any(missing.args)) {
    .err(sprintf("%s must be supplied",
                 word_list(names(missing.args)[missing.args])))
  }

  if (any(args.not.list)) {
    .err(sprintf("%s must be %s",
                 word_list(names(args.not.list)[args.not.list]),
                 ngettext(sum(args.not.list), "a list", "lists")))
  }

  if (!force && length(covs.list) > 1L) {
    .err("optweights are currently not valid for longitudinal treatments. Set `force = TRUE` to bypass this message at your own risk")
  }

  treat.types <- vapply(treat.list, function(x) {
    if (chk::vld_character_or_factor(x) || is_binary(x)) "cat"
    else "cont"
  }, character(1L))

  treat.list <- lapply(seq_along(treat.types), function(x) {
    if (treat.types[x] == "cat") as.character(treat.list[[x]])
    else as.numeric(treat.list[[x]])
  })

  if (!all_apply(covs.list, function(c) all(apply(c, 2L, is.numeric)))) {
    .err("all covariates must be numeric")
  }

  covs.list <- lapply(covs.list, as.matrix)
  bin.covs.list <- lapply(covs.list, is_binary_col)

  times <- seq_along(covs.list)

  tols.list <- tols
  if (length(tols.list) == 1L) {
    tols.list <- replicate(length(times), tols.list[[1L]], simplify = FALSE)
  }

  for (i in which(lengths(tols.list) == 1L)) {
    tols.list[[i]] <- rep.int(tols.list[[i]], ncol(covs.list[[i]]))
  }

  chk::chk_string(norm)
  norm <- tolower(norm)
  norm.options <- c("l2", "l1", "linf")
  chk::chk_subset(norm, norm.options)

  N <- nrow(covs.list[[1L]])
  sw <- {
    if (is_null(s.weights)) rep.int(1, N)
    else s.weights
  }

  if (is_null(b.weights)) {
    bw <- rep.int(1, N)
  }
  else if (norm == "l2") {
    bw <- b.weights
  }
  else {
    .err("only the l2 norm is compatible with `b.weights`")
  }

  estimand <- toupper(estimand)

  chk::chk_number(min.w)
  chk::chk_lt(min.w, 1)

  if (length(times) > 1L && is_not_null(estimand) && !identical(estimand, "ATE")) {
    .err("only the ATE or specified targets are compatible with longitduinal treatments")
  }

  unique.treats <- lapply(times, function(i) {
    switch(treat.types[i],
           cat = sort(unique(treat.list[[i]])),
           "treat")
  })

  n <- lapply(times, function(i) {
    switch(treat.types[i],
           cat = vapply(unique.treats[[i]],
                        function(t) sum(treat.list[[i]] == t),
                        numeric(1L)),
           c(treat = N))
  })

  means <- lapply(covs.list, col.w.m, w = sw)
  if (is_not_null(estimand) && (is_null(targets) || all(is.na(targets)))) {
    targets <- sds <- make_list(length(times))
    if (estimand %in% c("ATT", "ATC")) {
      for (i in times) {
        in_focal <- which(treat.list[[i]] == focal)

        sds[[i]] <- sqrt(col.w.v(covs.list[[i]][in_focal, , drop = FALSE],
                                 w = sw[in_focal], bin.vars = bin.covs.list[[i]]))

        targets[[i]] <- {
          if (i == 1L) col.w.m(covs.list[[i]][in_focal, , drop = FALSE],
                               w = sw[in_focal])
          else rep.int(NA_real_, ncol(covs.list[[i]]))
        }
      }

      sw[in_focal] <- 1
    }
    else if (estimand == "ATE") {
      for (i in times) {
        sds[[i]] <- sqrt(colMeans(do.call("rbind", lapply(unique.treats[[i]], function(t) {
          in.treat <- switch(treat.types[i],
                             cat = which(treat.list[[i]] == t),
                             cont = seq_along(treat.list[[i]]))

          col.w.v(covs.list[[i]][in.treat, , drop = FALSE],
                  w = sw[in.treat], bin.vars = bin.covs.list[[i]])
        }))))

        targets[[i]] <- {
          if (i == 1L) means[[i]]
          else rep.int(NA_real_, ncol(covs.list[[i]]))
        }
      }
    }
  }
  else {
    if (is_null(targets)) {
      targets <- rep(NA_real_, ncol(covs.list[[1L]]))
    }
    else if (!is.atomic(targets) || (!all(is.na(targets)) && !is.numeric(targets))) {
      .err("`targets` must be a vector of target values for each baseline covariate")
    }

    if (length(targets) != ncol(covs.list[[1L]])) {
      .err("`targets` must have the same number of values as there are baseline covariates")
    }

    targets <- c(list(targets), lapply(covs.list[-1L], function(c) rep(NA_real_, ncol(c))))

    sds <- make_list(length(times))
    for (i in times) {
      sds[[i]] <- sqrt(colMeans(do.call("rbind", lapply(unique.treats[[i]], function(t) {
        in.treat <- switch(treat.types[i],
                           cat = which(treat.list[[i]] == t),
                           cont = seq_along(treat.list[[i]]))

        col.w.v(covs.list[[i]][in.treat, , drop = FALSE],
                w = sw[in.treat], bin.vars = bin.covs.list[[i]])
      }))))
    }
  }

  targeted <- balanced <- treat.sds <- treat.means <- tols <- make_list(length(times))
  for (i in times) {
    if (treat.types[i] == "cat") {
      targeted[[i]] <- !is.na(targets[[i]])
      balanced[[i]] <- !targeted[[i]]

      treat.sds[[i]] <- NA_real_
      treat.means[[i]] <- NA_real_

      #tols
      vars.to.standardize <- {
        if (std.binary && std.cont) rep.int(TRUE, length(tols.list[[i]]))
        else if (!std.binary && std.cont) !bin.covs.list[[i]]
        else if (std.binary && !std.cont) bin.covs.list[[i]]
        else rep.int(FALSE, length(tols.list[[i]]))
      }

      tols[[i]] <- abs(tols.list[[i]])

      covs.list[[i]][, vars.to.standardize & !check_if_zero(tols.list[[i]]) & !check_if_zero(sds[[i]])] <-
        mat_div(covs.list[[i]][, vars.to.standardize & !check_if_zero(tols.list[[i]]) & !check_if_zero(sds[[i]]), drop = FALSE],
                sds[[i]][vars.to.standardize & !check_if_zero(tols.list[[i]]) & !check_if_zero(sds[[i]])])
      targets[[i]][vars.to.standardize & !check_if_zero(tols.list[[i]]) & !check_if_zero(sds[[i]])] <-
        targets[[i]][vars.to.standardize & !check_if_zero(tols.list[[i]]) & !check_if_zero(sds[[i]])] /
        sds[[i]][vars.to.standardize & !check_if_zero(tols.list[[i]]) & !check_if_zero(sds[[i]])]

      #Note: duals work incorrecly unless tols are > 0, so replace small tols with
      #sqrt(.Machine$double.eps).
      # tols[[i]] <- ifelse(tols[[i]] < sqrt(.Machine$double.eps),
      #                     sqrt(.Machine$double.eps),
      #                     tols[[i]])
    }
    else {
      sds[[i]] <- sqrt(col.w.v(covs.list[[i]], w = sw, bin.vars = bin.covs.list[[i]]))
      targeted[[i]] <- !is.na(targets[[i]])
      balanced[[i]] <- rep.int(TRUE, length(targeted[[i]]))
      covs.list[[i]][, targeted[[i]]] <- center(covs.list[[i]][, targeted[[i]], drop = FALSE],
                                                at = targets[[i]][targeted[[i]]]) #center covs at targets (which will be eventual means)
      covs.list[[i]][, !targeted[[i]]] <- center(covs.list[[i]][, !targeted[[i]], drop = FALSE],
                                                 at = means[[i]][!targeted[[i]]]) #center covs at means

      treat.sds[[i]] <- sqrt(col.w.v(treat.list[[i]], w = sw))
      treat.means[[i]] <- col.w.m(matrix(treat.list[[i]], ncol = 1L), w = sw)
      treat.list[[i]] <- treat.list[[i]] - treat.means[[i]] #center treat

      # tols[[i]] <- abs(tols.list[[i]]*sds[[i]]*treat.sds[[i]])

      tols[[i]] <- abs(tols.list[[i]])
      covs.list[[i]][, !check_if_zero(sds[[i]])] <- mat_div(covs.list[[i]][,!check_if_zero(sds[[i]]), drop = FALSE], sds[[i]][!check_if_zero(sds[[i]])])
      targets[[i]][!check_if_zero(sds[[i]])] <- targets[[i]][!check_if_zero(sds[[i]])] / sds[[i]][!check_if_zero(sds[[i]])]
      treat.list[[i]] <- treat.list[[i]]/treat.sds[[i]]

      #Note: duals work incorrecly unless tols are > 0, so replace small tols with
      #sqrt(.Machine$double.eps).
      # tols[[i]] <- ifelse(tols[[i]] < sqrt(.Machine$double.eps),
      #                     sqrt(.Machine$double.eps),
      #                     tols[[i]])
    }
  }

  if (norm == "l2") {
    #Minimizing variance of weights
    P <- Matrix::sparseMatrix(1:N, 1:N, x = 2*(sw^2)/N)
    # q = -sw/N #ensures objective function value is variance of weights
    q <- (-2*bw + mean(bw^2))*sw/N

    #Minimizing the sum of the variances in each treatment group
    #Note: equiv. to setting targets closer to smaller group
    # P = sparseMatrix(1:N, 1:N, x = (2*sw^2)/ifelse(treat.list[[1]]==1, n[[1]]["1"], n[[1]]["0"]))
    # q = -sw/ifelse(treat.list[[1]]==1, n[[1]]["1"], n[[1]]["0"]) #ensures objective function value is variance of weights

    #Mean of weights in each treat must equal 1
    A_meanw = do.call("rbind", lapply(times, function(i) {
      if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t) (treat.list[[i]] == t) * sw / n[[i]][t]))
      else sw/n[[i]]
    }))
    L_meanw = do.call("c", lapply(times, function(i) rep(1, length(unique.treats[[i]]))))
    U_meanw = L_meanw

    #All weights must be >= min.w; focal weights must be 1, weights where sw = 0 must be 0
    A_wmin = Matrix::sparseMatrix(1:N, 1:N, x = 1)
    if (is_not_null(focal)) {
      L_wmin <- ifelse(check_if_zero(sw), min.w, ifelse(treat.list[[1]] == focal, 1, min.w))
      U_wmin <- ifelse(check_if_zero(sw), min.w, ifelse(treat.list[[1]] == focal, 1, Inf))
    }
    else {
      L_wmin <- rep(min.w, N)
      U_wmin <- ifelse(check_if_zero(sw), min.w, Inf)
    }

    #Targeting constraints
    #Note: need 2 * in order to simulate tols/2 but using original tols.
    #This makes dual variables work as expected.
    A_target = do.call("rbind", lapply(times, function(i) {
      if (any(targeted[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t)
          if (is_null(focal)) 2 * t(covs.list[[i]][, targeted[[i]], drop = FALSE] * (treat.list[[i]] == t) * sw / n[[i]][t])
          else if (is_not_null(focal) && t != focal) t(covs.list[[i]][, targeted[[i]], drop = FALSE] * (treat.list[[i]] == t) * sw / n[[i]][t])
        ))
        else rbind(t(covs.list[[i]][, targeted[[i]], drop = FALSE] * sw), treat.list[[i]] * sw) #variables are centered
      }
      else NULL
    }))
    L_target = do.call("c", lapply(times, function(i) {
      if (any(targeted[[i]])) {
        if (treat.types[i] == "cat") do.call("c", lapply(unique.treats[[i]], function(t) {
          if (is_null(focal)) 2 * targets[[i]][targeted[[i]]] - tols[[i]][targeted[[i]]]
          else if (is_not_null(focal) && t != focal) targets[[i]][targeted[[i]]] - tols[[i]][targeted[[i]]]
        }))
        else rep(0, sum(targeted[[i]]) + 1) #variables are centered at targets; +1 for treat
      }
      else NULL
    }))
    U_target = do.call("c", lapply(times, function(i) {
      if (any(targeted[[i]])) {
        if (treat.types[i] == "cat") do.call("c", lapply(unique.treats[[i]], function(t) {
          if (is_null(focal)) 2 * targets[[i]][targeted[[i]]] + tols[[i]][targeted[[i]]]
          else if (is_not_null(focal) && t != focal) targets[[i]][targeted[[i]]] + tols[[i]][targeted[[i]]]
        }))
        else rep(0, sum(targeted[[i]]) + 1) #variables are centered at targets; +1 for treat
      }
      else NULL
    }))

    #Balancing constraints for all covariates
    A_balance <- do.call("rbind", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(combn(unique.treats[[i]], 2, simplify = FALSE), function(comb) {
          t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[1]) * sw / n[[i]][comb[1]]) - t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[2]) * sw / n[[i]][comb[2]])
        }))
        else {
          correct.factor <- 2 #see w.cov
          if (corr.type == "pearson")  t(covs.list[[i]][, balanced[[i]], drop = FALSE] * treat.list[[i]] * sw / (n[[i]] - correct.factor)) #For cont, all have balancing constraints
          else if (corr.type == "spearman")  t(apply(covs.list[[i]][, balanced[[i]], drop = FALSE], 2, rank) * treat.list[[i]] * sw / (n[[i]] - correct.factor)) #For cont, all have balancing constraints
          else {
            rbind(t(covs.list[[i]][, balanced[[i]], drop = FALSE] * treat.list[[i]] * sw / (n[[i]] - correct.factor)),
                  t(apply(covs.list[[i]][, balanced[[i]], drop = FALSE], 2, rank) * treat.list[[i]] * sw / (n[[i]] - correct.factor)))
          }
        }
      }
      else NULL

    }))
    L_balance <- do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(-tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else {
          if (corr.type %in% c("pearson", "spearman")) -tols[[i]][balanced[[i]]]
          else rep(-tols[[i]][balanced[[i]]], 2)
        }

      }
      else NULL
    }))
    U_balance <- do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else {
          if (corr.type %in% c("pearson", "spearman")) tols[[i]][balanced[[i]]]
          else rep(tols[[i]][balanced[[i]]], 2)
        }
      }
      else NULL
    }))

    A  <- rbind(A_wmin, A_meanw, A_balance, A_target)
    L <- c(L_wmin, L_meanw, L_balance, L_target)
    U <- c(U_wmin, U_meanw, U_balance, U_target)

    out <- osqp::solve_osqp(P = P, q = q, A = A, l = L, u = U,
                            pars = do.call(osqp::osqpSettings, args))

    #Get dual vars for balance and target constraints
    A_balance.indices <- {
      if (is_null(A_balance)) NULL
      else (NROW(A_wmin) + NROW(A_meanw) + 1L):(NROW(A_wmin) + NROW(A_meanw) + NROW(A_balance))
    }
    A_target.indices <- {
      if (is_null(A_target)) NULL
      else (NROW(A_wmin) + NROW(A_meanw) + NROW(A_balance) + 1L):(NROW(A_wmin) + NROW(A_meanw) + NROW(A_balance) + NROW(A_target))
    }

    w <- out$x
  }
  else if (norm == "l1") {
    #Minimizing mean absolute deviation of weights
    P <- sparseMatrix(NULL, NULL, dims = c(2*N, 2*N))
    q <- c(rep(0, N), 2*sw/N)

    #Mean of weights in each treat must equal 1
    A_meanw <- do.call("rbind", lapply(times, function(i) {
      switch(treat.types[i],
             cat = do.call("rbind", lapply(unique.treats[[i]], function(t) (treat.list[[i]] == t) * sw / n[[i]][t])),
             sw/n[[i]])
    }))
    L_meanw <- unlist(lapply(times, function(i) rep.int(1L, length(unique.treats[[i]]))))
    U_meanw <- L_meanw

    #All weights must be >= min; focal weights must be 1, weights where sw = 0 must be 0
    #Auxilliary vars must be >= 0
    min <- min.w
    A_wmin <- sparseMatrix(1:(2*N), 1:(2*N), x = 1)
    if (is_not_null(focal)) {
      L_wmin <- ifelse(check_if_zero(sw), min, ifelse(treat.list[[1L]] == focal, 1, min))
      U_wmin <- ifelse(check_if_zero(sw), min, ifelse(treat.list[[1L]] == focal, 1, Inf))
    }
    else {
      L_wmin <- rep(min, N)
      U_wmin <- ifelse(check_if_zero(sw), min, Inf)
    }
    Lz_wmin <- c(L_wmin, rep(0, N))
    Uz_wmin <- c(U_wmin, rep(Inf, N))

    #Targeting constraints
    A_target = do.call("rbind", lapply(times, function(i) {
      if (any(targeted[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t)
          if (is_null(focal) || (is_not_null(focal) && t != focal)) t(covs.list[[i]][, targeted[[i]], drop = FALSE] * (treat.list[[i]] == t) * sw / n[[i]][t])
        ))
        else rbind(t(covs.list[[i]][, targeted[[i]], drop = FALSE] * sw), treat.list[[i]] * sw) #variables are centered
      }
      else NULL
    }))
    L_target = do.call("c", lapply(times, function(i) {
      if (any(targeted[[i]])) {
        if (treat.types[i] == "cat") do.call("c", lapply(unique.treats[[i]], function(t) {
          if (is_null(focal)) targets[[i]][targeted[[i]]] - tols[[i]][targeted[[i]]]/2
          else if (is_not_null(focal) && t != focal) targets[[i]][targeted[[i]]] - tols[[i]][targeted[[i]]]
        }))
        else rep(0, sum(targeted[[i]]) + 1) #variables are centered at targets; +1 for treat
      }
      else NULL
    }))
    U_target = do.call("c", lapply(times, function(i) {
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
    A_balance <- do.call("rbind", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(combn(unique.treats[[i]], 2, simplify = FALSE), function(comb) {
          t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[1]) * sw / n[[i]][comb[1]]) - t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[2]) * sw / n[[i]][comb[2]])
        }))
        else {
          correct.factor <- 2
          if (corr.type == "pearson")  t(covs.list[[i]][, balanced[[i]], drop = FALSE] * treat.list[[i]] * sw / (n[[i]] - correct.factor)) #For cont, all have balancing constraints
          else if (corr.type == "spearman")  t(apply(covs.list[[i]][, balanced[[i]], drop = FALSE], 2, rank) * treat.list[[i]] * sw / (n[[i]] - correct.factor)) #For cont, all have balancing constraints
          else {
            rbind(t(covs.list[[i]][, balanced[[i]], drop = FALSE] * treat.list[[i]] * sw / (n[[i]] - correct.factor)),
                  t(apply(covs.list[[i]][, balanced[[i]], drop = FALSE], 2, rank) * treat.list[[i]] * sw / (n[[i]] - correct.factor)))
          }
        }
      }
      else NULL

    }))
    L_balance <- do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(-tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else {
          if (corr.type %in% c("pearson", "spearman")) -tols[[i]][balanced[[i]]]
          else rep(-tols[[i]][balanced[[i]]], 2)
        }
      }
      else NULL
    }))
    U_balance <- do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else {
          if (corr.type %in% c("pearson", "spearman")) tols[[i]][balanced[[i]]]
          else rep(tols[[i]][balanced[[i]]], 2)
        }
      }
      else NULL
    }))

    #Conversion constraints
    Inxn = sparseMatrix(1:N, 1:N, x = 1)
    A_conversion = rbind(cbind(Inxn, -Inxn),
                         cbind(-Inxn, -Inxn))
    # A_conversion = sparseMatrix(c(1:N, 1:N, (N+1):(2*N), (N+1):(2*N)),
    #                  c(1:N, (N+1):(2*N), 1:N, (N+1):(2*N)),
    #                  x = c(rep(1, N), rep(-1, 3*N)))
    L_conversion = rep(-Inf, 2*N)
    U_conversion = rep(1, 2*N)

    A  <- rbind(A_meanw, A_balance, A_target)
    L <- c(L_meanw, L_balance, L_target)
    U <- c(U_meanw, U_balance, U_target)

    Au <- cbind(A, matrix(0, nrow = nrow(A), ncol = N))

    Az <- rbind(Au, A_wmin, A_conversion)
    Lz = c(L, Lz_wmin, L_conversion)
    Uz = c(U, Uz_wmin, U_conversion)

    out <- solve_osqp(P = P, q = q, A = Az, l = Lz, u = Uz,
                      pars = do.call(osqpSettings, args))

    w <- out$x[1:N]

    #Get dual vars for constraints
    A_balance.indices <- if (is_not_null(A_balance)) NROW(A_meanw) + 1:NROW(A_balance)
    A_target.indices <- if (is_not_null(A_target)) NROW(A_meanw) + NROW(A_balance) + 1:NROW(A_target)

  }
  else if (norm == "linf") {
    #Minimizing largest weight
    P = sparseMatrix(NULL, NULL, dims = c(2*N, 2*N))
    q = rep(sw/N, 2)

    #Mean of weights in each treat must equal 1
    A_meanw = do.call("rbind", lapply(times, function(i) {
      if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t) (treat.list[[i]] == t) * sw / n[[i]][t]))
      else sw/n[[i]]
    }))
    L_meanw = do.call("c", lapply(times, function(i) rep(1, length(unique.treats[[i]]))))
    U_meanw = L_meanw

    #All weights must be >= min; focal weights must be 1, weights where sw = 0 must be 0
    #Auxilliary var must be >= 0
    min <- min.w
    A_wmin = sparseMatrix(1:(2*N), 1:(2*N), x = 1)
    if (is_not_null(focal)) {
      L_wmin = ifelse(check_if_zero(sw), min, ifelse(treat.list[[1]] == focal, 1, min))
      U_wmin = ifelse(check_if_zero(sw), min, ifelse(treat.list[[1]] == focal, 1, Inf))
    }
    else {
      L_wmin = rep(min, N)
      U_wmin = ifelse(check_if_zero(sw), min, Inf)
    }
    Lz_wmin = c(L_wmin, rep(0, N))
    Uz_wmin = c(U_wmin, rep(Inf, N))

    #Targeting constraints
    A_target = do.call("rbind", lapply(times, function(i) {
      if (any(targeted[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t)
          if (is_null(focal) || (is_not_null(focal) && t != focal)) t(covs.list[[i]][, targeted[[i]], drop = FALSE] * (treat.list[[i]] == t) * sw / n[[i]][t])
        ))
        else rbind(t(covs.list[[i]][, targeted[[i]], drop = FALSE] * sw), treat.list[[i]] * sw) #variables are centered
      }
      else NULL
    }))
    L_target = do.call("c", lapply(times, function(i) {
      if (any(targeted[[i]])) {
        if (treat.types[i] == "cat") do.call("c", lapply(unique.treats[[i]], function(t) {
          if (is_null(focal)) targets[[i]][targeted[[i]]] - tols[[i]][targeted[[i]]]/2
          else if (is_not_null(focal) && t != focal) targets[[i]][targeted[[i]]] - tols[[i]][targeted[[i]]]
        }))
        else rep(0, sum(targeted[[i]]) + 1) #variables are centered at targets; +1 for treat
      }
      else NULL
    }))
    U_target = do.call("c", lapply(times, function(i) {
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
    A_balance = do.call("rbind", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(combn(unique.treats[[i]], 2, simplify = FALSE), function(comb) {
          t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[1]) * sw / n[[i]][comb[1]]) - t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[2]) * sw / n[[i]][comb[2]])
        }))
        else {
          correct.factor <- 2
          if (corr.type == "pearson")  t(covs.list[[i]][, balanced[[i]], drop = FALSE] * treat.list[[i]] * sw / (n[[i]] - correct.factor)) #For cont, all have balancing constraints
          else if (corr.type == "spearman")  t(apply(covs.list[[i]][, balanced[[i]], drop = FALSE], 2, rank) * treat.list[[i]] * sw / (n[[i]] - correct.factor)) #For cont, all have balancing constraints
          else {
            rbind(t(covs.list[[i]][, balanced[[i]], drop = FALSE] * treat.list[[i]] * sw / (n[[i]] - correct.factor)),
                  t(apply(covs.list[[i]][, balanced[[i]], drop = FALSE], 2, rank) * treat.list[[i]] * sw / (n[[i]] - correct.factor)))
          }
        }
      }
      else NULL

    }))
    L_balance = do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(-tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else {
          if (corr.type %in% c("pearson", "spearman")) -tols[[i]][balanced[[i]]]
          else rep(-tols[[i]][balanced[[i]]], 2)
        }
      }
      else NULL
    }))
    U_balance = do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else {
          if (corr.type %in% c("pearson", "spearman")) tols[[i]][balanced[[i]]]
          else rep(tols[[i]][balanced[[i]]], 2)
        }
      }
      else NULL
    }))


    #Conversion constraints
    Inxn = sparseMatrix(1:N, 1:N, x = 1)
    #one = matrix(1, nrow = N, ncol = 1)
    A_conversion1 = rbind(cbind(Inxn, -Inxn),
                          cbind(-Inxn, -Inxn))
    # A_conversion = sparseMatrix(c(1:N, 1:N, (N+1):(2*N), (N+1):(2*N)),
    #                  c(1:N, rep(N+1, N), 1:N, rep(N+1, N)),
    #                  x = c(rep(1, N), rep(-1, 3*N)))
    L_conversion1 = rep(-Inf, 2*N)
    U_conversion1 = rep(1, 2*N)

    A_conversion2 = cbind(sparseMatrix(NULL, NULL, dims = c(N-1, N)),
                          matrix(1, ncol = 1, nrow = N-1),
                          sparseMatrix(1:(N-1), 1:(N-1), x = -1))
    L_conversion2 = rep(0, N - 1)
    U_conversion2 = rep(0, N - 1)

    A_conversion = rbind(A_conversion1, A_conversion2)
    L_conversion = c(L_conversion1, L_conversion2)
    U_conversion = c(U_conversion1, U_conversion2)

    A = rbind(A_meanw, A_balance, A_target)
    L = c(L_meanw, L_balance, L_target)
    U = c(U_meanw, U_balance, U_target)

    Au = cbind(A, matrix(0, nrow = NROW(A), ncol = ncol(A)))

    Az = rbind(Au, A_wmin, A_conversion)
    Lz = c(L, Lz_wmin, L_conversion)
    Uz = c(U, Uz_wmin, U_conversion)

    out <- solve_osqp(P = P, q = q, A = Az, l = Lz, u = Uz,
                      pars = do.call(osqpSettings, args))

    w <- out$x[1:N]

    #Get dual vars for constraints
    A_balance.indices <- if (is_not_null(A_balance)) NROW(A_meanw) + 1:NROW(A_balance)
    A_target.indices <- if (is_not_null(A_target)) NROW(A_meanw) + NROW(A_balance) + 1:NROW(A_target)
  }

  w[w < min.w] <- min.w

  #Duals
  balance_duals <- abs(out$y[A_balance.indices]) #A_balance
  target_duals <- abs(out$y[A_target.indices]) #A_target
  duals <- make_list(length(times))

  kb <- kt <- 1L
  for (i in times) {
    td <- bd <- NULL

    if (is_not_null(target_duals)) {
      non.focal.treats <- {
        if (is_null(focal)) unique.treats[[i]]
        else setdiff(unique.treats[[i]], focal)
      }

      targeted.covs <- colnames(covs.list[[i]])[targeted[[i]]]
      if (is_not_null(non.focal.treats) && is_not_null(targeted.covs)) {
        td <- data.frame(expand.grid(constraint = "target",
                                     cov = targeted.covs,
                                     treat = non.focal.treats,
                                     stringsAsFactors = FALSE),
                         dual = target_duals[kt:(kt + length(non.focal.treats) * length(targeted.covs) - 1L)]
        )
      }

      kt <- kt + length(non.focal.treats) * length(targeted.covs)
    }

    if (is_not_null(balance_duals)) {
      treat.combs <- switch(treat.types[i],
                            cat = vapply(combn(unique.treats[[i]], 2L, simplify = FALSE),
                                         paste, character(1L), collapse = " vs. "),
                            unique.treats[[i]])

      balanced.covs <- colnames(covs.list[[i]])[balanced[[i]]]
      if (is_not_null(treat.combs) && is_not_null(balanced.covs)) {
        bd <- data.frame(expand.grid(constraint = "balance",
                                     cov = balanced.covs,
                                     treat = treat.combs,
                                     stringsAsFactors = FALSE),
                         dual = balance_duals[kb:(kb + length(treat.combs) * length(balanced.covs) - 1)]
        )
      }

      kb <- kb + length(treat.combs) * length(balanced.covs)
    }

    if (is_not_null(td) || is_not_null(bd)) {
      duals[[i]] <- rbind(td, bd)
    }
  }

  opt_out <- list(w = w,
                  duals = duals,
                  info = out$info,
                  out = out,
                  A = A)
  class(opt_out) <- "optweight.fit"

  opt_out
}
