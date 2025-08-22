#Simulation study
library(cobalt); library(WeightIt); library(optweight)
library(jtools); library(ATE)
nrep <- 500
npred <- 12
n <- 2e3
estimand <- "ATE"
cor_ <- .15

between <- function(x, lo, hi) {
  x >= lo & x <=hi
}

types <- c("true"
           , "reg1"
           , "reg2"
           , "ps"
           , "ebal"
           , "ebcw"
           , "ow1"
           , "ow2"
           , "ow3"
)

# est.list <- lapply(seq_len(nrep),
onerun <- function(i) {

  X <- as.data.frame(MASS::mvrnorm(n, rep(0, npred), cor_ + (1-cor_)*diag(npred)))

  #PS model
  ps_r2 <- .25
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



  est.list <- setNames(vector("list", length(types)), types)

  #effect estimATE and CI
  if ("true" %in% types) {
    #True ps
    true_w <- Z/P + (1-Z)/(1-P)
    if (estimand == "ATT") true_w <- P*true_w

    true.est.fit <- glm(Y ~ Z, data = d, weights = true_w)
    est.list[["true"]] <- summ(true.est.fit, robust = "HC1", confint = TRUE,
                               model.fit = FALSE, model.info = FALSE,
                               pval = FALSE)$coef["Z",1:3]

  }

  #effect estimATE and CI
  if ("reg1" %in% types) {
    #Reg no int
    f1 <- as.formula(paste("Y ~", paste(c("Z", colnames(X)), collapse = " + ")))
    reg1.est.fit <- glm(f1, data = d)
    est.list[["reg1"]] <- summ(reg1.est.fit, robust = FALSE, confint = TRUE,
                               model.fit = FALSE, model.info = FALSE,
                               pval = FALSE)$coef["Z",1:3]

  }

  if ("reg2" %in% types) {
    #Reg w/int
    d1 <- d; d1[colnames(X)] <- scale(X)
    f2 <- as.formula(paste("Y ~", paste("Z", colnames(X), sep = "*", collapse = " + ")))
    reg2.est.fit <- glm(f2, data = d1)
    est.list[["reg2"]] <- summ(reg2.est.fit, robust = "HC1", confint = TRUE,
                               model.fit = FALSE, model.info = FALSE,
                               pval = FALSE)$coef["Z",1:3]

  }

  if ("ps" %in% types) {
    #Logistic PS
    ps.fit <- weightit(f.build("Z", X), data = d, method = "ps", estimand = estimand)
    ps_w <- ps.fit$weights

    ps.est.fit <- glm(Y ~ Z, data = d, weights = ps_w)
    est.list[["ps"]] <- summ(ps.est.fit, robust = "HC1", confint = TRUE,
                             model.fit = FALSE, model.info = FALSE,
                             pval = FALSE)$coef["Z",1:3]
  }

  if ("ebal" %in% types) {
    #ebal w
    ebal.fit <- weightit(f.build("Z", X), data = d, method = "ebal", estimand = estimand)
    ebal_w <- ebal.fit$weights

    ebal.est.fit <- glm(Y ~ Z, data = d, weights = ebal_w)
    est.list[["ebal"]] <- summ(ebal.est.fit, robust = "HC1", confint = TRUE,
                               model.fit = FALSE, model.info = FALSE,
                               pval = FALSE)$coef["Z",1:3]
  }

  if ("ebcw" %in% types) {
    ebcw.fit <- weightit(f.build("Z", X), data = d, method = "ebcw", estimand = estimand)
    ebcw_w <- ebcw.fit$weights

    ebcw.est.fit <- ATE(Y = Y, Ti = Z, X = X, ATT = (estimand == "ATT"))
    est.list[["ebcw"]] <- summary(ebcw.est.fit)$Est[3, c(1,3,4)]
  }

  if ("ow1" %in% types) {
    #optweight exact
    ow.fit_1 <- optweight(f.build("Z", X), data = d, tols = 0, estimand = estimand)
    ow_1 <- ow.fit_1$weights

    ow1.est.fit <- glm(Y ~ Z, data = d, weights = ow_1)
    est.list[["ow1"]] <- summ(ow1.est.fit, robust = "HC1", confint = TRUE,
                              model.fit = FALSE, model.info = FALSE,
                              pval = FALSE)$coef["Z",1:3]
  }

  if ("ow2" %in% types) {
    #optweight relaxed
    ow.fit_2 <- optweight(f.build("Z", X), data = d, tols = .01, estimand = estimand)
    ow_2 <- ow.fit_2$weights

    ow2.est.fit <- glm(Y ~ Z, data = d, weights = ow_2)
    est.list[["ow2"]] <- summ(ow2.est.fit, robust = "HC1", confint = TRUE,
                              model.fit = FALSE, model.info = FALSE,
                              pval = FALSE)$coef["Z",1:3]
  }

  if ("ow3" %in% types) {
    #optweight with knowledge
    ow.fit_3 <- optweight(f.build("Z", X), data = d, tols = 10^pmin(-as.numeric(cut(abs(Y0_coefs), 4)), -2),
                          estimand = estimand)
    ow_3 <- ow.fit_3$weights

    ow3.est.fit <- glm(Y ~ Z, data = d, weights = ow_3)
    est.list[["ow3"]] <- summ(ow3.est.fit, robust = "HC1", confint = TRUE,
                              model.fit = FALSE, model.info = FALSE,
                              pval = FALSE)$coef["Z",1:3]
  }

  # bal.tab(X, Z, weights = data.frame(true = true_w,
  #                                    ps = ps_w,
  #                                    ebal = ebal_w,
  #                                    ebcw = ebcw_w
  #                                    ow1 = ow_1,
  #                                    ow2 = ow_2,
  #                                    ow3 = ow_3),
  #         estimand = estimand, method = "w")

  est <- setNames(data.frame(types, do.call('rbind', est.list)),
                  c("type", "est", "cilo", "cihi"))

  return(est)
}

est.out.list <- vector("list", nrep)
pb <- txtProgressBar(0, nrep, style = 3)
for (i in seq_len(nrep)) {
  setTxtProgressBar(pb, i)
  est.out.list[[i]] <- onerun(i)
}

est <- do.call("rbind", est.out.list)

#Bias
round(sapply(types, function(t) {
  mean(est[est$type==t, "est"] - 2)
}), 3)

#RMSE
round(sapply(types, function(t) {
  sqrt(mean((est[est$type==t, "est"] - 2)^2))
}), 3)

#Coverage prob
sapply(types, function(t) {
  with(est[est$type==t,], mean(between(2, cilo, cihi)))
})

#Power
sapply(types, function(t) {
  with(est[est$type==t,], mean(!between(0, cilo, cihi)))
})
