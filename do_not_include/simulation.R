#Simulation study
library(cobalt); library(WeightIt); library(optweight)
library(jtools)
nrep <- 100
npred <- 10
n <- 1e3
estimand <- "ATE"
cor_ <- .1

between <- function(x, lo, hi) {
  x >= lo & x <=hi
}

types <- c("true", "ps", "ebal", "ebcw", "ow1", "ow2", "ow3")

est.list <- lapply(seq_len(nrep), function(i) {

  X <- as.data.frame(MASS::mvrnorm(n, rep(0, npred), cor_ + (1-cor_)*diag(npred)))


  #PS model
  ps_r2 <- .75
  kps <- sqrt((ps_r2 * pi^2/3) / ((1-ps_r2) * sum(((1:npred)/10)^2)))
  ps_coefs <-  kps * cumprod(rep(-1, npred)) * (1:npred)/10
  Z_eta <- -.8 + as.matrix(X) %*% ps_coefs
  P <- plogis(Z_eta)
  Z <- rbinom(n, 1, P)

  #Outcome model
  Y_r2 <- .35
  kY <- sqrt((Y_r2 * 100) / ((1-Y_r2) * sum((1:npred)^2)))
  Y0_coefs <- kY * -cumprod(rep(-1, npred)) * rev(seq_len(npred))
  Y1_coefs <- Y0_coefs + (cumprod(rep(-1, npred)) + 1) * sd(Y0_coefs)
  Y0 <- as.matrix(X) %*% Y0_coefs + rnorm(n, 0, 10)
  Y1 <- 2 + as.matrix(X) %*% Y1_coefs + rnorm(n, 0, 10)

  Y <- Z * Y1 + (1-Z) * Y0
  d <- data.frame(X, Z, Y)
  #True ps
  true_w <- Z/P + (1-Z)/(1-P)
  if (estimand == "ATT") true_w <- P*true_w

  #Logistic PS
  ps.fit <- weightit(f.build("Z", X), data = d, method = "ps", estimand = estimand)
  ps_w <- ps.fit$weights

  #ebal w
  ebal.fit <- weightit(f.build("Z", X), data = d, method = "ebal", estimand = estimand)
  ebal_w <- ebal.fit$weights

  #optweight exact
  ow.fit_1 <- optweight(f.build("Z", X), data = d, tols = 0, estimand = estimand)
  ow_1 <- ow.fit_1$weights

  #optweight relaxed
  ow.fit_2 <- optweight(f.build("Z", X), data = d, tols = .01, estimand = estimand)
  ow_2 <- ow.fit_2$weights

  #optweight with knowledge
  ow.fit_3 <- optweight(f.build("Z", X), data = d, tols = 10^pmin(-as.numeric(cut(abs(Y0_coefs), 4)), -2),
                        estimand = estimand)
  ow_3 <- ow.fit_3$weights

  # bal.tab(X, Z, weights = data.frame(true = true_w,
  #                                    ps = ps_w,
  #                                    ebal = ebal_w,
  #                                    ow1 = ow_1,
  #                                    ow2 = ow_2,
  #                                    ow3 = ow_3),
  #         estimand = estimand, method = "w")

  #effect estimATE and CI
  true.est.fit <- glm(Y ~ Z, data = d, weights = true_w)
  true.est <- summ(true.est.fit, robust = "HC1", confint = TRUE,
                 model.fit = FALSE, model.info = FALSE,
                 pval = FALSE)$coef["Z",1:3]

  ps.est.fit <- glm(Y ~ Z, data = d, weights = ps_w)
  ps.est <- summ(ps.est.fit, robust = "HC1", confint = TRUE,
                 model.fit = FALSE, model.info = FALSE,
                 pval = FALSE)$coef["Z",1:3]

  ebal.est.fit <- glm(Y ~ Z, data = d, weights = ebal_w)
  ebal.est <- summ(ebal.est.fit, robust = "HC1", confint = TRUE,
                   model.fit = FALSE, model.info = FALSE,
                   pval = FALSE)$coef["Z",1:3]

  ebcw.est.fit <- ATE(Y = Y, Ti = Z, X = X, ATT = (estimand == "ATT"))
  ebcw.est <- summary(ebcw.est.fit)$Est[3, c(1,3,4)]

  ow1.est.fit <- glm(Y ~ Z, data = d, weights = ow_1)
  ow1.est <- summ(ow1.est.fit, robust = "HC1", confint = TRUE,
                  model.fit = FALSE, model.info = FALSE,
                  pval = FALSE)$coef["Z",1:3]

  ow2.est.fit <- glm(Y ~ Z, data = d, weights = ow_2)
  ow2.est <- summ(ow2.est.fit, robust = "HC1", confint = TRUE,
                  model.fit = FALSE, model.info = FALSE,
                  pval = FALSE)$coef["Z",1:3]

  ow3.est.fit <- glm(Y ~ Z, data = d, weights = ow_3)
  ow3.est <- summ(ow3.est.fit, robust = "HC1", confint = TRUE,
                  model.fit = FALSE, model.info = FALSE,
                  pval = FALSE)$coef["Z",1:3]


  est <- setNames(data.frame(types,
                             rbind(true.est, ps.est, ebal.est, ebcw.est, ow1.est, ow2.est, ow3.est)),
                  c("type", "est", "cilo", "cihi"))

  return(est)
})

est <- do.call("rbind", est.list)

sapply(types, function(t) {
  mean(est[est$type==t, "est"] - 2)
})

sapply(types, function(t) {
  sqrt(mean((est[est$type==t, "est"] - 2)^2))
})
sapply(types, function(t) {
  with(est[est$type==t,], mean(between(2, cilo, cihi)))
})
sapply(types, function(t) {
  with(est[est$type==t,], mean(!between(0, cilo, cihi)))
})
