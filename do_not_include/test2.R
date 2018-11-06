opt <- function(A, B, E, F, G, H, method = "quadprog") {
  #library(Dykstra)

  Neq  <- nrow(E)

  dvec  <- B # was: t(A)%*%B
  Dmat  <- A # t(A) %*% A
  diag(Dmat) <- diag(Dmat)+1e-8
  Amat  <- t(rbind(E,G))
  bvec  <- c(F,H)
  upper <- c(F, rep(Inf, length(H)))

  if (method == "quadprog") {
    sol   <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq=Neq)
  }
  else if (method == "dykstra") {
    sol <- Dykstra::dykstra(Dmat, dvec, Amat, bvec, meq = Neq)
  }
  else if (method == "lsei") {
    sol <- list(solution = lsei::qp(q = Dmat, p = -dvec, c = E, d = F, e = G, f = H, lower = 0,
                                    tol = sqrt(.Machine$double.eps)))
  }
  else if (method == "osqp") {

    out <- rosqp::solve_osqp(P = Dmat, q = -dvec, A = t(Amat), l = bvec, u = upper,
                             pars = rosqp::osqpSettings(adaptive_rho = FALSE,
                                                        #rho = 5e-3,
                                                        max_iter = 200000L,
                                                        verbose = FALSE))
    sol <- list(solution = out$x)
  }

  sol$X <- sol$solution
  return(sol)
}

#696.466 104.081

library(NlcOptim)
covs <- replicate(5, rnorm(1E3))
t <- as.numeric((covs %*% runif(ncol(covs), -1.5, 1.5) + rnorm(nrow(covs), sd = 2)) > .3)

objfun <- function(x) {
  bal <- apply(covs, 2, function(c) (weighted.mean(c[t==1], x[t==1]) - weighted.mean(c[t==0], x[t==0])) /
                 sqrt(mean(var(c[t==1]) + var(c[t==0]))))
  return(mean(bal^2))
}
confun <- function(x) {
  ESS <- 700
  v <- length(x)/ESS - 1
  #c <- length(x) / (1 + mean((x - 1) ^ 2))
  c <- sum(t==1)/(1 + mean((x[t==1]-1)^2)) + sum(t==0)/(1 + mean((x[t==0]-1)^2))
  #c <- mean((x-1)^2)
  return(list(ceq = NULL, c = -c + ESS))
}
lb <- rep(0, nrow(covs))
ub <- rep(Inf, nrow(covs))
Aeq <- rbind(matrix(1, nrow = 1, ncol = nrow(covs)),
             t(covs/nrow(covs)))
Beq <- c(nrow(covs),
         colMeans(covs))

out <- NlcOptim::solnl(rep(1, nrow(covs)), objfun = objfun,
             confun = confun,
             Aeq = Aeq,
             Beq = Beq,
             lb = lb, ub = ub,
             tolX = 1E-11,
             tolFun = 1E-11,
             tolCon = 1E-11,
             maxIter = 2E6)

objfun <- function(x) mean((x-1)^2)
confun <- function(x) {
  c <- abs(apply(covs, 2, function(c) (weighted.mean(c[t==1], x[t==1]) - weighted.mean(c[t==0], x[t==0])) /
                   sqrt(mean(var(c[t==1]) + var(c[t==0])))))
  tols = rep(0, ncol(covs))
  return(list(ceq = NULL, c = max(c - tols)))
}

covs <- lalonde[-c(1,4,9)]
var <- "married"
tols <- check.tols(treat ~ covs + I(educ^2) + I(educ^3) + I(re74^2), data = lalonde, tols = 0)
tols1 <- tols2 <- tols
tols1[var] <- .0
tols2[var] <- tols1[var] + .001
ow1 <- optweight(treat ~ covs + I(educ^2) + I(educ^3) + I(re74^2), data = lalonde, estimand = "ATT", tols = tols1, focal = 1)
ow2 <- optweight(treat ~ covs + I(educ^2) + I(educ^3) + I(re74^2), data = lalonde, estimand = "ATT", tols = tols2, focal = 1)
wv1 <- ow1$info$obj_val
wv2 <- ow2$info$obj_val
dual1 <- ow1$duals[ow1$duals$cov==var, "dual"]
dual2 <- ow2$duals[ow2$duals$cov==var, "dual"]

print(c(dual1 = dual1,
        dual2 = dual2,
        wv1 = wv1,
        wv2 = wv2,
        diff = wv1-wv2,
        exp = (tols2[var]-tols1[var])*mean(c(dual1, dual2))))

tols <- check.tols(treat ~ covs[-c(3)] + I(educ^2) + I(educ^3) + I(re74^2) + I(re75^2), data = lalonde, tols = 0)
d <- data.frame(tol = seq(.000, .1, length=26), var = NA, dual = NA)
for (i in 1:nrow(d)) {
  cat(i, "...")
  tols["married"] <- d$tol[i]
  ow_ <- optweight(treat ~ covs[-c(3)] + I(educ^2) + I(educ^3) + I(re74^2) + I(re75^2), data = lalonde, estimand = "ATT", tols = tols)
  if (ow_$info$status_val == 1 ) {
    d$var[i] <- summary(ow_)$coef[2]^2
    d$stddual[i] <- ow_$duals["married", 1]
    d$obj[i] <- ow_$info$obj_val
  }
}
d$ess <- 429/(1+d$var)
#########

times <- 1
covs <- lalonde[-c(1, 4, 9)]
N <- nrow(covs)
treat <- lalonde$treat
sw <- rep(1, N)
tols <- rep(0, ncol(covs))
args <- list()

P = sparseMatrix(1:N, 1:N, x = 2*sw^2)
q = rep(-1, N)

#Mean of weights in each treat must equal 1
E1 = rbind(treat/sum(treat), (1 - treat)/(sum(1-treat)))
F1l = c(1, 1)
F1u = F1l

#All weights must be >= 0
G1 = sparseMatrix(1:N, 1:N, x = 1)
H1l = rep(0, N)
H1u = rep(Inf, N)

#Balance constraints
G2 = t(covs) %*% (diag(treat/sum(treat)) - diag((1-treat)/sum(1-treat)))
H2l = -tols
H2u = tols

#Target specific value
G3 = rbind(t(covs) %*% diag(treat/sum(treat)),
           t(covs) %*% diag((1-treat)/sum(1-treat)))
H3l = rep(targets, 2)
H3u = H3l

#Process args
args[names(args) %nin% names(formals(rosqp::osqpSettings))] <- NULL
if (is_null(args[["adaptive_rho"]])) args[["adaptive_rho"]] <- TRUE
if (is_null(args[["max_iter"]])) args[["max_iter"]] <- 2E5
if (is_null(args[["eps_abs"]])) args[["eps_abs"]] <- 1E-9
if (is_null(args[["eps_rel"]])) args[["eps_rel"]] <- 1E-9
#args[["verbose"]] <- verbose

A  <- rbind(G1, E1, G3, G2)
lower <- c(H1l, F1l, H3l, H2l)
upper <- c(H1u, F1u, H3u, H2u)

out <- rosqp::solve_osqp(P = P, q = q, A = A, l = lower, u = upper,
                         pars = do.call(rosqp::osqpSettings, args))
