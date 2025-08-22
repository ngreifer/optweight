# Log liklihood score equations
# Score equation for cannonical link is always t(X * (A - mean(A)))

n <- 1e3

X1 <- rnorm(n)
X2 <- rnorm(n, 0, 2)
X <- cbind(1, X1, X2)

# Logistic regression

A <- rbinom(n, 1, (1 + exp(-(X1 + X2)))^-1)

Score <- t(X * (A - mean(A)))

# Linear regression, homo

A <- .3 + X1 - .5*X2 + rnorm(n)

Score <- rbind(t(X * (A - mean(A))/var(A)),
               -1 + ((A - mean(A))^2/var(A)))

# Linear regression, hetero

A <- .3 + X1 - .5*X2 + rnorm(n, 0, sqrt(exp(.5*X1 + .7*X2)))

Score <- rbind(t(X * (A - mean(A))/var(A)),
               t(X * (-1 + ((A - mean(A))^2/var(A)))))

# Negative binomial

A <- rpois(n, exp(.3 + X1 - .5*X2))

fit <- MASS::glm.nb(A ~ 1)
mu = exp(coef(fit)[1])
theta = 1

Score <- rbind(
  t(X * (A - mu)/(1+theta*mu)),
  t(X/theta * (theta*(A - mu)/(1+theta*mu) + log(1 + theta*mu) - digamma(A + 1/theta) + digamma(1/theta)))
)

# Poisson
A <- rpois(n, exp(.3 + X1 - .5*X2))

# fit <- glm(A ~ 1, family = poisson)
# mu = exp(coef(fit)[1])
mu = mean(A)

Score <- t(X * (A - mu))

#-----

max_iter <- 1e6
min.w <- 1e-8

#Constraints

tols <- rep(1e-6, nrow(Score))

#Sum of w must equal n
A_sumw <- Matrix::Matrix(rep(1, n), ncol = n)
l_sumw <- n
u_sumw <- n

#All weights must be positive
A_posw <- Matrix::sparseMatrix(1:n, 1:n, x = rep(1,n))
l_posw <- rep(min.w, n)
u_posw <- rep(Inf, n)

#|Cov(A,X)| must be less than AXtols
A_score <- Score
l_score <- -tols
u_score <- tols

A_ <- rbind(A_sumw, A_posw, A_score)
L <- c(l_sumw, l_posw, l_score)
U <- c(u_sumw, u_posw, u_score)

#Objective function
#Minimize the variance of the weights
#Note: need 2 in P matrix because osqp multiplies by .5

P = Matrix::sparseMatrix(1:n, 1:n, x = rep(2/n, n))

q = rep(-1/n, n) #ensures objective function value is variance of weights

#Optimize
options.list <- list(eps_abs = 1e-8, eps_rel = 1e-8,
                     verbose = FALSE, max_iter = as.integer(max_iter))
opt.out <- osqp::solve_osqp(P, q, A_, L, U,
                            pars = do.call(osqp::osqpSettings, options.list))

w <- opt.out$x
w[w < min.w] <- min.w #Numerical imprecision yields negative weights






con <- function(w, X, A) {
  crossprod(w, X) %*% (A - mean(A))
}

lb <- rep(0, n)
# ub <- rep(Inf, n)
Aeq <- rbind(matrix(as.numeric(A == 1), ncol = n)
             ,matrix(as.numeric(A == 0), ncol = n)
             , t(X[,-1] * as.numeric(A == 1)) / sum(A==1)
             , t(X[,-1] * as.numeric(A == 0)) / sum(A==0)
)
Beq <- c(sum(A == 1),
         , sum(A == 0)
         , col_w_mean(X[,-1])
         , col_w_mean(X[,-1])
)

objfun <- function(x) mean((x-1)^2)
confun <- function(x) {
  ceq <- c(0
    # ,col_w_mean(X[,-1], x, subset = A == 1) - col_w_mean(X[,-1])
    # ,col_w_mean(X[,-1], x, subset = A == 0) - col_w_mean(X[,-1])
  )
  cin <- c(0
         # , col_w_ks(X, A, x) - rep(.1, ncol(X))
  )
  return(list(ceq = ceq, c = cin))
}

out <- NlcOptim::solnl(rep(1, n),
                       objfun = objfun,
                       # confun = confun,
                       Aeq = Aeq,
                       Beq = Beq,
                       lb = lb, ub = ub,
                       tolX = 1E-3,
                       tolFun = 1E-3,
                       tolCon = 1E-3,
                       maxIter = 2E6)

col_w_mean(X[,-1], out$par, subset = A == 1) - col_w_mean(X[,-1])
col_w_mean(X[,-1], out$par, subset = A == 0) - col_w_mean(X[,-1])

col_w_ks(X, A, out$par)
