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
                 tols = (1:9)/100))
bal.tab(ow)

(ow <- optweight(treat ~ covs + I(re74==0) + I(re75==0), data = lalonde, estimand = "ATE"))
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
