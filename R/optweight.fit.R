optweight.fit <- function(treat.list, covs.list, tols, estimand = "ATE", targets = NULL, s.weights = NULL, focal = NULL, norm = "l2", std.binary = FALSE, std.cont = TRUE, min.w = 1E-8, verbose = FALSE, force = FALSE, ...) {

  args <- list(...)

  #Process args
  args[names(args) %nin% names(formals(rosqp::osqpSettings))] <- NULL
  if (is_null(args[["max_iter"]])) args[["max_iter"]] <- 2E5L
  if (is_null(args[["eps_abs"]])) args[["eps_abs"]] <- 1E-8
  if (is_null(args[["eps_rel"]])) args[["eps_rel"]] <- 1E-8
  args[["verbose"]] <- verbose

  key.args <- c("treat.list", "covs.list", "tols")
  missing.args <- args.not.list <- setNames(rep(FALSE, length(key.args)), key.args)
  for (arg in key.args) {
    if (eval(substitute(missing(q), list(q = arg)))) {
      missing.args[arg] <- TRUE
    }
    else if (!is.vector(get(arg), mode = "list")) {
      args.not.list[arg] <- TRUE
    }
  }
  if (any(missing.args)) stop(paste(word.list(names(missing.args)[missing.args]), "must be supplied."), call. = FALSE)
  if (any(args.not.list)) stop(paste(word.list(names(args.not.list)[args.not.list]), "must be", ifelse(sum(args.not.list) > 1, "lists.", "a list.")), call. = FALSE)

  if (length(covs.list) > 1 && !force) stop("Optweights are currently not valid for longitudinal treatments. Set force = TRUE to bypass this message at your own risk.", call. = FALSE)

  treat.types <- vapply(treat.list, function(x) {
    if (is.factor(x) || is.character(x) || is_binary(x)) "cat"
    else "cont"
  }, character(1L))
  treat.list <- lapply(seq_along(treat.types), function(x) {
    if (treat.types[x] == "cat") as.character(treat.list[[x]])
    else as.numeric(treat.list[[x]])
  })
  if (!all(vapply(covs.list, function(c) all(apply(c, 2, is.numeric)), logical(1L)))) stop("All covariates must be numeric.", call. = FALSE)
  covs.list <- lapply(covs.list, as.matrix)

  times <- seq_along(covs.list)

  tols.list <- tols
  if (length(tols.list) == 1) tols.list <- replicate(max(times), tols.list[[1]], simplify = FALSE)
  tols.list <- lapply(times, function(i) if (length(tols.list[[i]]) == 1) rep(tols.list[[i]], ncol(covs.list[[i]])) else tols.list[[i]])

  N <- nrow(covs.list[[1]])
  if (is_null(s.weights)) sw <- rep(1, N)
  else sw <- s.weights

  norm.options <- c("l2", "l1", "linf")
  if (length(norm) != 1 || !is.character(norm) || tolower(norm) %nin% norm.options) {
    stop(paste0("norm must be ", word.list(norm.options, and.or = "or", quotes = TRUE), "."), call. = FALSE)
  }
  else norm <- tolower(norm)

  estimand <- toupper(estimand)

  if (length(min.w) != 1 || !is.numeric(min.w) || min.w < 0 || min.w >= 1) stop("min.w must be a single number in the interval [0, 1).", call. = FALSE)

  if (length(times) > 1 && is_not_null(estimand) && estimand %nin% "ATE") stop("Only the ATE or specified targets are compatible with longitduinal treatments.", call. = FALSE)

  unique.treats <- lapply(times, function(i) {
    if (treat.types[i] == "cat") sort(unique(treat.list[[i]]))
    else "treat"
  })
  n <- lapply(times, function(i) {
    if (treat.types[i] == "cat") vapply(unique.treats[[i]],
                                        function(t) sum(treat.list[[i]] == t),
                                        numeric(1)) #faster than tapply
    else c(treat = N)
  })

  means <- lapply(covs.list, col.w.m, w = sw)
  if (is_not_null(estimand) && (is_null(targets) || all(is.na(targets)))) {
    if (estimand %in% c("ATT", "ATC")) {
      targets <- lapply(times, function(i) {
        if (i == 1) {
          if (is_null(focal)) focal <- max(treat.list[[i]])
          else if (estimand == "ATC") focal <- min(treat.list[[i]])
          col.w.m(covs.list[[i]][treat.list[[i]] == focal, , drop = FALSE], w = sw[treat.list[[i]] == focal])
        }
        else rep(NA_real_, ncol(covs.list[[i]]))
      })
      sds <- lapply(times, function(i) {
        if (is_null(focal)) focal <- max(treat.list[[i]])
        else if (estimand == "ATC") focal <- min(treat.list[[i]])
        sqrt(col.w.v(covs.list[[i]][treat.list[[i]] == focal, , drop = FALSE], w = sw[treat.list[[i]] == focal]))
      })
      sw[treat.list[[1]]==focal] <- 1
    }
    else if (estimand == "ATE") {
      targets <- c(list(means[[1]]), lapply(covs.list[-1], function(c) rep(NA_real_, ncol(c))))
      sds <- lapply(times, function(i) sqrt(rowMeans(matrix(sapply(unique.treats[[i]], function(t) col.w.v(covs.list[[i]][treat.list[[i]]==t, , drop = FALSE], w = sw[treat.list[[i]] == t]), simplify = "array"), ncol = length(unique.treats[[i]])))))
    }
  }
  else {
    if (is_null(targets)) targets <- rep(NA_real_, ncol(covs.list[[1]]))
    else if (!is.atomic(targets) || (!all(is.na(targets)) && !is.numeric(targets))) stop("targets must be a vector of target values for each baseline covariate.", call. = FALSE)

    targets <- c(list(targets), lapply(covs.list[-1], function(c) rep(NA_real_, ncol(c))))

    if (length(targets[[1]]) != ncol(covs.list[[1]])) {
      stop("targets must have the same number of values as there are baseline covariates.", call. = FALSE)
    }
    sds <- lapply(times, function(i) sqrt(rowMeans(matrix(sapply(unique.treats[[i]], function(t) col.w.v(covs.list[[i]][treat.list[[i]]==t, , drop = FALSE], w = sw[treat.list[[i]] == t]), simplify = "array"), ncol = length(unique.treats[[i]])))))

  }

  targeted <- balanced <- treat.sds <- treat.means <- tols <- vector("list", length(times))
  for (i in times) {
    if (treat.types[i] == "cat") {
      targeted[[i]] <- !is.na(targets[[i]])
      balanced[[i]] <- !targeted[[i]]
      #balanced[[i]] <- rep(TRUE, length(targeted[[i]]))
      treat.sds[[i]] <- NA_real_
      treat.means[[i]] <- NA_real_

      #tols
      if (std.binary && std.cont) vars.to.standardize <- rep(TRUE, length(tols.list[[i]]))
      else if (!std.binary && std.cont) vars.to.standardize <- !apply(covs.list[[i]], 2, is_binary)
      else if (std.binary && !std.cont) vars.to.standardize <- apply(covs.list[[i]], 2, is_binary)
      else vars.to.standardize <- rep(FALSE, length(tols.list[[i]]))

      tols[[i]] <- ifelse(vars.to.standardize & !check_if_zero(tols.list[[i]]),
                          abs(tols.list[[i]]*sds[[i]]), #standardize
                          abs(tols.list[[i]]))

      #Note: duals work incorrecly unless tols are > 0, so replace small tols with
      #sqrt(.Machine$double.eps).
      # tols[[i]] <- ifelse(tols[[i]] < sqrt(.Machine$double.eps),
      #                     sqrt(.Machine$double.eps),
      #                     tols[[i]])
    }
    else {
      targeted[[i]] <- !is.na(targets[[i]])
      balanced[[i]] <- rep(TRUE, length(targeted[[i]]))
      covs.list[[i]][, targeted[[i]]] <- sweep(covs.list[[i]][, targeted[[i]], drop = FALSE], 2, targets[[i]][targeted[[i]]], "-") #center covs at targets (which will be eventual means)
      covs.list[[i]][, !targeted[[i]]] <- sweep(covs.list[[i]][, !targeted[[i]], drop = FALSE], 2, means[[i]][!targeted[[i]]], "-") #center covs at means
      sds[[i]] <- sqrt(col.w.v(covs.list[[i]], w = sw))
      treat.sds[[i]] <- sqrt(col.w.v(matrix(treat.list[[i]], ncol = 1), w = sw))
      treat.means[[i]] <- col.w.m(matrix(treat.list[[i]], ncol = 1), w = sw)
      treat.list[[i]] <- treat.list[[i]] - treat.means[[i]] #center treat

      tols[[i]] <- abs(tols.list[[i]]*sds[[i]]*treat.sds[[i]])
      #Note: duals work incorrecly unless tols are > 0, so replace small tols with
      #sqrt(.Machine$double.eps).
      # tols[[i]] <- ifelse(tols[[i]] < sqrt(.Machine$double.eps),
      #                     sqrt(.Machine$double.eps),
      #                     tols[[i]])
    }
  }

  if (norm == "l2") {
    #Minimizing variance of weights
    P = sparseMatrix(1:N, 1:N, x = 2*(sw^2)/N)
    q = -sw/N #ensures objective function value is variance of weights

    #Minimizing the sum of the variances in each treatment group
    #Note: equiv. to setting targets closer to smaller group
    # P = sparseMatrix(1:N, 1:N, x = (2*sw^2)/ifelse(treat.list[[1]]==1, n[[1]]["1"], n[[1]]["0"]))
    # q = -sw/ifelse(treat.list[[1]]==1, n[[1]]["1"], n[[1]]["0"]) #ensures objective function value is variance of weights

    #Mean of weights in each treat must equal 1
    E1 = do.call("rbind", lapply(times, function(i) {
      if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t) (treat.list[[i]] == t) * sw / n[[i]][t]))
      else sw/n[[i]]
    }))
    F1l = do.call("c", lapply(times, function(i) rep(1, length(unique.treats[[i]]))))
    F1u = F1l

    #All weights must be >= min; focal weights must be 1, weights where sw = 0 must be 0
    min <- min.w
    G1 = sparseMatrix(1:N, 1:N, x = 1)
    if (is_not_null(focal)) {
      H1l <- ifelse(check_if_zero(sw), min, ifelse(treat.list[[1]] == focal, 1, min))
      H1u <- ifelse(check_if_zero(sw), min, ifelse(treat.list[[1]] == focal, 1, Inf))
    }
    else {
      H1l <- rep(min, N)
      H1u <- ifelse(check_if_zero(sw), min, Inf)
    }

    #Targeting constraints
    #Note: need 2 * in order to simulate tols/2 but using original tols.
    #This makes dual variables work as expected.
    G2 = do.call("rbind", lapply(times, function(i) {
      if (any(targeted[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t)
          if (is_null(focal)) 2 * t(covs.list[[i]][, targeted[[i]], drop = FALSE] * (treat.list[[i]] == t) * sw / n[[i]][t])
          else if (is_not_null(focal) && t != focal) t(covs.list[[i]][, targeted[[i]], drop = FALSE] * (treat.list[[i]] == t) * sw / n[[i]][t])
        ))
        else rbind(t(covs.list[[i]][, targeted[[i]], drop = FALSE] * sw), treat.list[[i]] * sw) #variables are centered
      }
      else NULL
    }))
    H2l = do.call("c", lapply(times, function(i) {
      if (any(targeted[[i]])) {
        if (treat.types[i] == "cat") do.call("c", lapply(unique.treats[[i]], function(t) {
          if (is_null(focal)) 2 * targets[[i]][targeted[[i]]] - tols[[i]][targeted[[i]]]
          else if (is_not_null(focal) && t != focal) targets[[i]][targeted[[i]]] - tols[[i]][targeted[[i]]]
        }))
        else rep(0, sum(targeted[[i]]) + 1) #variables are centered at targets; +1 for treat
      }
      else NULL
    }))
    H2u = do.call("c", lapply(times, function(i) {
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
    G3 <- do.call("rbind", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(combn(unique.treats[[i]], 2, simplify = FALSE), function(comb) {
          t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[1]) * sw / n[[i]][comb[1]]) - t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[2]) * sw / n[[i]][comb[2]])
        }))
        else t(covs.list[[i]][, balanced[[i]], drop = FALSE] * treat.list[[i]] * sw / n[[i]]) #For cont, all have balancing constraints
      }
      else NULL

    }))
    H3l <- do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(-tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else -tols[[i]][balanced[[i]]]
      }
      else NULL
    }))
    H3u <- do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else tols[[i]][balanced[[i]]]
      }
      else NULL
    }))

    A  <- rbind(G1, E1, G3, G2)
    lower <- c(H1l, F1l, H3l, H2l)
    upper <- c(H1u, F1u, H3u, H2u)

    out <- rosqp::solve_osqp(P = P, q = q, A = A, l = lower, u = upper,
                             pars = do.call(rosqp::osqpSettings, args))

    #Get dual vars for balance and target constraints
    G3.indices <- if (is_null(G3)) NULL else (NROW(G1)+NROW(E1)+1):(NROW(G1)+NROW(E1)+NROW(G3))
    G2.indices <- if (is_null(G2)) NULL else (NROW(G1)+NROW(E1)+NROW(G3)+1):(NROW(G1)+NROW(E1)+NROW(G3)+NROW(G2))

    w <- out$x
  }
  else if (norm == "l1") {
    #Minimizing mean absolute deviation of weights
    P = sparseMatrix(NULL, NULL, dims = c(2*N, 2*N))
    q = c(rep(0, N), 2*sw/N)

    #Mean of weights in each treat must equal 1
    E1 = do.call("rbind", lapply(times, function(i) {
      if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t) (treat.list[[i]] == t) * sw / n[[i]][t]))
      else sw/n[[i]]
    }))
    F1l = do.call("c", lapply(times, function(i) rep(1, length(unique.treats[[i]]))))
    F1u = F1l

    #All weights must be >= min; focal weights must be 1, weights where sw = 0 must be 0
    #Auxilliary vars must be >= 0
    min <- min.w
    G1 = sparseMatrix(1:(2*N), 1:(2*N), x = 1)
    if (is_not_null(focal)) {
      H1l <- ifelse(check_if_zero(sw), min, ifelse(treat.list[[1]] == focal, 1, min))
      H1u <- ifelse(check_if_zero(sw), min, ifelse(treat.list[[1]] == focal, 1, Inf))
    }
    else {
      H1l <- rep(min, N)
      H1u <- ifelse(check_if_zero(sw), min, Inf)
    }
    H1lz <- c(H1l, rep(0, N))
    H1uz <- c(H1u, rep(Inf, N))

    #Targeting constraints
    G2 = do.call("rbind", lapply(times, function(i) {
      if (any(targeted[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t)
          if (is_null(focal) || (is_not_null(focal) && t != focal)) t(covs.list[[i]][, targeted[[i]], drop = FALSE] * (treat.list[[i]] == t) * sw / n[[i]][t])
        ))
        else rbind(t(covs.list[[i]][, targeted[[i]], drop = FALSE] * sw), treat.list[[i]] * sw) #variables are centered
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
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(combn(unique.treats[[i]], 2, simplify = FALSE), function(comb) {
          t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[1]) * sw / n[[i]][comb[1]]) - t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[2]) * sw / n[[i]][comb[2]])
        }))
        else t(covs.list[[i]][, balanced[[i]], drop = FALSE] * treat.list[[i]] * sw / n[[i]]) #For cont, all have balancing constraints
      }
      else NULL

    }))
    H3l <- do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(-tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else -tols[[i]][balanced[[i]]]
      }
      else NULL
    }))
    H3u <- do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else tols[[i]][balanced[[i]]]
      }
      else NULL
    }))

    #Conversion constraints
    Inxn = sparseMatrix(1:N, 1:N, x = 1)
    I = rbind(cbind(Inxn, -Inxn),
              cbind(-Inxn, -Inxn))
    # I = sparseMatrix(c(1:N, 1:N, (N+1):(2*N), (N+1):(2*N)),
    #                  c(1:N, (N+1):(2*N), 1:N, (N+1):(2*N)),
    #                  x = c(rep(1, N), rep(-1, 3*N)))
    jl = rep(-Inf, 2*N)
    ju = rep(1, 2*N)

    A  <- rbind(E1, G3, G2)
    lower <- c(F1l, H3l, H2l)
    upper <- c(F1u, H3u, H2u)

    Au <- cbind(A, matrix(0, nrow = nrow(A), ncol = N))

    Az <- rbind(Au, G1, I)
    lowerz = c(lower, H1lz, jl)
    upperz = c(upper, H1uz, ju)

    out <- rosqp::solve_osqp(P = P, q = q, A = Az, l = lowerz, u = upperz,
                             pars = do.call(rosqp::osqpSettings, args))

    w <- out$x[1:N]

    #Get dual vars for constraints
    G3.indices <- if (is_null(G3)) NULL else (NROW(E1)+1):(NROW(E1)+NROW(G3))
    G2.indices <- if (is_null(G2)) NULL else (NROW(E1)+NROW(G3)+1):(NROW(E1)+NROW(G3)+NROW(G2))

  }
  else if (norm == "linf") {
    #Minimizing largest weight
    P = sparseMatrix(NULL, NULL, dims = c(2*N, 2*N))
    q = rep(sw/N, 2)

    #Mean of weights in each treat must equal 1
    E1 = do.call("rbind", lapply(times, function(i) {
      if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t) (treat.list[[i]] == t) * sw / n[[i]][t]))
      else sw/n[[i]]
    }))
    F1l = do.call("c", lapply(times, function(i) rep(1, length(unique.treats[[i]]))))
    F1u = F1l

    #All weights must be >= min; focal weights must be 1, weights where sw = 0 must be 0
    #Auxilliary var must be >= 0
    min <- min.w
    G1 = sparseMatrix(1:(2*N), 1:(2*N), x = 1)
    if (is_not_null(focal)) {
      H1l = ifelse(check_if_zero(sw), min, ifelse(treat.list[[1]] == focal, 1, min))
      H1u = ifelse(check_if_zero(sw), min, ifelse(treat.list[[1]] == focal, 1, Inf))
    }
    else {
      H1l = rep(min, N)
      H1u = ifelse(check_if_zero(sw), min, Inf)
    }
    H1lz = c(H1l, rep(0, N))
    H1uz = c(H1u, rep(Inf, N))

    #Targeting constraints
    G2 = do.call("rbind", lapply(times, function(i) {
      if (any(targeted[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(unique.treats[[i]], function(t)
          if (is_null(focal) || (is_not_null(focal) && t != focal)) t(covs.list[[i]][, targeted[[i]], drop = FALSE] * (treat.list[[i]] == t) * sw / n[[i]][t])
        ))
        else rbind(t(covs.list[[i]][, targeted[[i]], drop = FALSE] * sw), treat.list[[i]] * sw) #variables are centered
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
    G3 = do.call("rbind", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") do.call("rbind", lapply(combn(unique.treats[[i]], 2, simplify = FALSE), function(comb) {
          t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[1]) * sw / n[[i]][comb[1]]) - t(covs.list[[i]][, balanced[[i]], drop = FALSE] * (treat.list[[i]] == comb[2]) * sw / n[[i]][comb[2]])
        }))
        else t(covs.list[[i]][, balanced[[i]], drop = FALSE] * treat.list[[i]] * sw / n[[i]]) #For cont, all have balancing constraints
      }
      else NULL

    }))
    H3l = do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(-tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else -tols[[i]][balanced[[i]]]
      }
      else NULL
    }))
    H3u = do.call("c", lapply(times, function(i) {
      if (any(balanced[[i]])) {
        if (treat.types[i] == "cat") rep(tols[[i]][balanced[[i]]], length(combn(unique.treats[[i]], 2, simplify = FALSE)))
        else tols[[i]][balanced[[i]]]
      }
      else NULL
    }))


    #Conversion constraints
    Inxn = sparseMatrix(1:N, 1:N, x = 1)
    #one = matrix(1, nrow = N, ncol = 1)
    I = rbind(cbind(Inxn, -Inxn),
              cbind(-Inxn, -Inxn))
    # I = sparseMatrix(c(1:N, 1:N, (N+1):(2*N), (N+1):(2*N)),
    #                  c(1:N, rep(N+1, N), 1:N, rep(N+1, N)),
    #                  x = c(rep(1, N), rep(-1, 3*N)))
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

    A = rbind(E1, G3, G2)
    lower = c(F1l, H3l, H2l)
    upper = c(F1u, H3u, H2u)

    Au = cbind(A, matrix(0, nrow = NROW(A), ncol = ncol(A)))

    Az = rbind(Au, G1, I)
    lowerz = c(lower, H1lz, jl)
    upperz = c(upper, H1uz, ju)

    out <- rosqp::solve_osqp(P = P, q = q, A = Az, l = lowerz, u = upperz,
                             pars = do.call(rosqp::osqpSettings, args))

    w <- out$x[1:N]

    #Get dual vars for constraints
    G3.indices <- if (is_null(G3)) NULL else (NROW(E1)+1):(NROW(E1)+NROW(G3))
    G2.indices <- if (is_null(G2)) NULL else (NROW(E1)+NROW(G3)+1):(NROW(E1)+NROW(G3)+NROW(G2))
  }

  w[w < min.w] <- min.w

  #Duals
  balance_duals <- abs(out$y[G3.indices]) #G3
  target_duals <- abs(out$y[G2.indices]) #G2
  duals <- vector("list", length(times))

  kb <- kt <- 1
  for (i in times) {
    if (is_not_null(target_duals)) {
      non.focal.treats <- if (is_null(focal)) unique.treats[[i]] else unique.treats[[i]][unique.treats[[i]] != focal]

      targeted.covs <- colnames(covs.list[[i]])[targeted[[i]]]
      if (length(non.focal.treats) > 0 && length(targeted.covs) > 0) {
        td <- data.frame(expand.grid(constraint = "target",
                                     cov = targeted.covs,
                                     treat = non.focal.treats,
                                     stringsAsFactors = FALSE),
                         dual = target_duals[kt:(kt + length(non.focal.treats) * length(targeted.covs) - 1)]
        )
      }
      else td <- NULL

      kt <- kt + length(non.focal.treats) * length(targeted.covs)
    }
    else td <- NULL
    if (is_not_null(balance_duals)) {
      if (treat.types[i] == "cat") treat.combs <- vapply(combn(unique.treats[[i]], 2, simplify = FALSE),
                                                         paste, character(1L), collapse = " vs. ")
      else treat.combs <- unique.treats[[i]]
      balanced.covs <- colnames(covs.list[[i]])[balanced[[i]]]
      if (length(treat.combs) > 0 && length(balanced.covs) > 0) {
        bd <- data.frame(expand.grid(constraint = "balance",
                                     cov = balanced.covs,
                                     treat = treat.combs,
                                     stringsAsFactors = FALSE),
                         dual = balance_duals[kb:(kb + length(treat.combs) * length(balanced.covs) - 1)]
        )
      }
      else bd <- NULL
      kb <- kb + length(treat.combs) * length(balanced.covs)
    }
    else bd <- NULL

    duals[[i]] <- rbind(td, bd)
  }

  opt_out <- list(w = w,
                  duals = duals,
                  info = out$info,
                  out = out,
                  A = A)
  class(opt_out) <- "optweight.fit"

  return(opt_out)
}
