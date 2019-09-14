for (i in dir("R/")) source(paste0("R/", i))
library(Matrix)
library(ggplot2)
stop("Done sourcing.", call. = FALSE)

library("cobalt")
data("lalonde", package = "cobalt")
covs <- lalonde[-c(1, 9)]

system.time(ow <- optweight(lalonde$treat ~ covs, estimand = "ATE"))
bal.tab(ow)
summary(ow)
plot(summary(ow))
plot(ow)

(ow <- optweight(treat ~ covs, data = lalonde, estimand = "ATT"))
bal.tab(ow)
summary(ow)
plot(summary(ow))
plot(ow)

(ow <- optweight(lalonde$treat ~ covs, estimand = "ATT",
                 tols = (1:7)/100))
bal.tab(ow)

tols <- check.tols(treat ~ covs + I(re74==0) + I(re75==0), data = lalonde,
                   tols = rep(.002, 8), stop = T)

(ow <- optweight(treat ~ covs + I(re74==0) + I(re75==0), data = lalonde, estimand = "ATE",
                 tols = tols))
bal.tab(ow)

(ow <- optweight(race ~ covs[-3], data = lalonde, estimand = "ATE"))
bal.tab(ow)

(ow <- optweight(race ~ covs[-3], data = lalonde, estimand = "ATT", focal = "hispan",
                 tols = .01))
bal.tab(ow)



#Continuous Treatment
(ow <- optweight(re78 ~ covs, data = lalonde, estimand = "ATE",
                 tols = .001))
bal.tab(ow)
plot(ow)

library(twang)
data("iptwExWide", package = "twang")

owmsm <- optweight(list(tx1 ~ age + gender + use0,
                        tx2 ~ age + gender + use0 + use1 + tx1),
                   data = iptwExWide, tols = list(c(.01, .02, .03),
                                                  c(.001, .002, .003, .1, .1)))
bal.tab(owmsm)
summary(owmsm)
plot(summary(owmsm))
plot(owmsm, which.time = 2)


system.time(ow <- optweight(f.build("treat", covs), data = do.call('rbind', replicate(20, lalonde, F)),
                            estimand = "ATE", tols = .01))
bal.tab(ow)

target.covs <- c("age", "educ", "race")
targets <- check.targets(f.build("", target.covs), data = lalonde,
                         c(29, 9, .3, .2, .6))
ow.s <- optweight.svy(f.build("", target.covs), data = lalonde, targets = targets)
ow.s
summary(ow.s)
apply(splitfactor(lalonde[target.covs], drop.first = FALSE), 2, mean)
apply(splitfactor(lalonde[target.covs], drop.first = FALSE), 2, weighted.mean, w = ow.s$weights)

ow1 <- optweight(f.build("treat", target.covs), data = lalonde,
                 estimand = "ATE", s.weights = ow.s$weights, tols = 0)
bal.tab(ow1, disp.means = TRUE)

ow2 <- optweight(f.build("treat", target.covs), data = lalonde,
                 estimand = NULL, tols = 0,
                 targets = targets)
bal.tab(ow2, disp.means = TRUE)

#Generate data
N_ <- 1e3
npred_ <- 10
cor_ <- .2

#Generate predictors
var_X <- 1
X <- MASS::mvrnorm(N_, mu = rep(0, npred_),
                   Sigma = var_X*(cor_ + (1-cor_)*diag(npred_)))
colnames(X) <- paste0("X", seq_len(npred_))

#Generate treatment
#Strong effects account for 2/3 of explained variance
#Weak effects account for 1/3 of explained variance

A_r2 = .5

A_strong_covs <- colnames(X)[rep(c(TRUE, FALSE), each = npred_/2)]
A_weak_covs <- setdiff(colnames(X), A_strong_covs)
A_r2_strong <- (2/3) * A_r2
A_r2_weak <- (1/3) * A_r2

A_var_e <- 9
A_coef_strong <- sqrt(A_var_e*A_r2_strong/(length(A_strong_covs)*var_X*(1-A_r2)))
A_coef_weak <- sqrt(A_var_e*A_r2_weak/(length(A_weak_covs)*var_X*(1-A_r2)))
A_coef <- c(rep(A_coef_strong, length(A_strong_covs)),
            rep(A_coef_weak, length(A_weak_covs)))

A <- drop(X[,c(A_strong_covs, A_weak_covs)] %*% A_coef + rnorm(N_, 0, sqrt(A_var_e)))

A <- as.numeric(A>0)
ow <- optweight(A ~ X, tols = .1)
