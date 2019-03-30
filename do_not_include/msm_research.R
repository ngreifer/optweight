#Censoring simulation
library(cobalt); library(WeightIt); library(optweight)
library(jtools)
#nrep <- 100
npred <- 3
n <- 1e6
estimand <- "ATE"
cor_ <- .1

X <- as.data.frame(MASS::mvrnorm(n, rep(0, npred), cor_ + (1-cor_)*diag(npred)))

#PS model
ps_r2 <- .5
kps <- sqrt((ps_r2 * pi^2/3) / ((1-ps_r2) * sum(((1:npred)/10)^2)))
ps_coefs <-  kps * cumprod(rep(-1, npred)) * (1:npred)/10
Z_eta <- -.8 + as.matrix(X) %*% ps_coefs
P <- plogis(Z_eta)
Z <- rbinom(n, 1, P)

#Censoring model
# cps_r2 <- .5
# kcps <- sqrt((cps_r2 * pi^2/3) / ((1-cps_r2) * sum(((1:(npred+1))/10)^2)))
# cps_coefs <-  kcps * cumprod(rep(-1, npred)) * (1:npred)/10
cps_coefs <- c(-.4, -.9, -.7, .6)
C_eta <- -.8 + cbind(as.matrix(X), Z) %*% cps_coefs
cP <- plogis(C_eta)
C <- rbinom(n, 1, cP)

#Outcome model
Y_r2 <- .25
kY <- sqrt((Y_r2 * 100) / ((1-Y_r2) * sum((1:npred)^2)))
Y0_coefs <- kY * -cumprod(rep(-1, npred)) * rev(seq_len(npred))
Y1_coefs <- Y0_coefs + (cumprod(rep(-1, npred)) + 1) * sd(Y0_coefs)
Y0 <- as.matrix(X) %*% Y0_coefs + rnorm(n, 0, 10)
Y1 <- 2 + as.matrix(X) %*% Y1_coefs + rnorm(n, 0, 10)

Y <- Z * Y1 + (1-Z) * Y0
Y[C==1] <- NA

d <- data.frame(X, Z, Y, C)

ps <- glm(f.build("Z", X), data = d, family = "binomial")$fitted
tw <- Z/ps + (1-Z)/(1-ps)

cps <- glm(update(f.build("C", X), . ~ . + Z), data = d, family = "binomial")$fitted
cw <- C/cps + (1-C)/(1-cps)

w <- tw * cw

ps.est.fit <- glm(Y ~ Z, data = d, weights = w)
(ps.est <- summ(ps.est.fit, robust = "HC1", confint = TRUE,
               model.fit = FALSE, model.info = FALSE,
               pval = FALSE)$coef["Z",1:3])



#set seed to replicate results
#Note: True effects: .14, .4
set.seed(12345)
library(jtools); library(cobalt)
#define sample size
n <- 5e4
#define confounder c
d1 <- rnorm(n,0,1)
#define treatment at time 1 as function of confounder
t1 <- as.numeric(.1*d1 + rnorm(n,0, sqrt(.99)) > .2)
#define censoring at time 2 as function of confounder and treat1
c2 <- as.numeric(.2*d1 - .3*t1 + rnorm(n,0, sqrt(.8)) > .7)
#define depression at time 2 as function of confounder and treat1
d2 <- .1*d1 + .4*t1 + rnorm(n,0, sqrt(.822))
#define treatment at time 2 as function of confounder and dep1
t2 <- as.numeric(.1*d1 + .4*d2 + .4*t1 + rnorm(n,0, sqrt(.5196)) > 0)
#define censoring at time 3 as function of dep1, treat2, and dep2
c3 <- ifelse(c2 == 1, 1, as.numeric(.1*d1 - .7*t2 + .4*d2 + rnorm(n,0, sqrt(.8)) > .7))
#define outcome depression at time 3 as function of dep1, treat2, and dep2
d3 <- .1*d1 + .4*t2 + .4*d2 + rnorm(n,0, sqrt(.4582))

d3c <- ifelse(c3 == 1, NA, d3)

#add ID variable to do mixed effects models later
id <- rep(1: length(d1))
#put all in a dataframe and write data to harddrive to use later in e.g. SPSS
df1 <- data.frame(id, d1, t1, c2, d2, t2, c3, d3, d3c)

#No censoring
ps1.fit <- glm(t1 ~ d1, data = df1, family = binomial(link = "probit"))
ps1 <- ps1.fit$fitted.values
ps1.den <- t1*ps1 + (1-t1)*(1-ps1)

sps1.fit <- glm(t1 ~ 1, data = df1, family = binomial(link = "probit"))
sps1 <- sps1.fit$fitted.values
ps1.num <- t1*sps1 + (1-t1)*(1-sps1)

ps2.fit <- glm(t2 ~ d1 + d2 + t1, data = df1, family = binomial(link = "probit"))
ps2 <- ps2.fit$fitted.values
ps2.den <- t2*ps2 + (1-t2)*(1-ps2)

sps2.fit <- glm(t2 ~ t1, data = df1, family = binomial(link = "probit"))
sps2 <- sps2.fit$fitted.values
ps2.num <- t2*sps2 + (1-t2)*(1-sps2)

w1 <- 1 / (ps1.den * ps2.den)
sw1 <- ps1.num * ps2.num * w1

bal.tab(list(t1~d1, t2~d1+d2+t1), data = df1, weights = list(w1=w1, sw1=sw1), un = T, quick = T)

summ(glm(d3~t1+t2, data = df1, weights = w1), robust = "HC0", digits = 3)

ow <- optweight(list(t1 ~ d1, t2 ~ d1 + d2), data = df1, tols = 0)

ow1 <- optweight(list(t1 ~ d1), data = df1, tols = 0)
ow2 <- optweight(list(t2 ~ d1 + d2), data = df1, tols = 0)
ow_ <- ow1$weights*ow2$weights

ow2a <- optweight(list(t2 ~ t1 + d1 + d2), data = df1, tols = 0, s.weights = ow1$weights)
ow_a <- ow2a$weights * ow1$weights

ow2b <- optweight(list(t2 ~ t1 + d1 + d2), data = df1, tols = 0, targets = c(col.w.m(cbind(t1,d1), w1), NA))
ow_b <- ow2b$weights * ow1$weights

bal.tab(list(t1~d1, t2~d1+d2+t1), data = df1, weights = list(ps = w1, ow_ = ow_, ow_a = ow_a),
        un = T, disp.means = T, quick = T)

summ(glm(d3c~t1+t2, data = df1, weights = w1), robust = "HC0", digits = 3)

#With censoring
ps1 <- glm(t1 ~ d1, data = df1, family = binomial(link = "probit"))$fitted
tw1 <- t1/ps1 + (1-t1)/(1-ps1)

cps2 <- glm(c2 ~ d1 + t1, data = df1, family = binomial(link = "probit"))$fitted
cw2 <- c2/cps2 + (1-c2)/(1-cps2)

ps2 <- glm(t2 ~ d1 + d2 + t1, data = df1[c2==0,], family = binomial(link = "probit"))$fitted
tw2 <- t2[c2==0]/ps2 + (1-t2[c2==0])/(1-ps2)

cps3 <- glm(c3 ~ d1 + d2 + t2, data = df1[c2==0,], family = binomial(link = "probit"))$fitted
cw3 <- c3[c2==0]/cps3 + (1-c3[c2==0])/(1-cps3)

w2 <- tw1[c3==0]*cw2[c3==0]*tw2[c3[c2==0]==0]*cw3[c3[c2==0]==0]

summ(glm(d3c~t1+t2, data = df1[c3==0,], weights = w2), robust = "HC0", digits = 3)
