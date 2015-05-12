



# ------------------------------------------------
# 
#  Estimating model DF
# 
# ------------------------------------------------


if (Sys.info()["user"] == "millaco") {
  setwd("~/Dropbox/SarahColin/PhD/capture_prob_paper")    
  library(setwidth)
} else 
if (Sys.info()["user"] == "millarc") {
  setwd("C:/work/repos/papers/capture_prop_paper/")
} else 
if (Sys.info()["user"] == "Millarc") {
  setwd("C:/work/SMFS-report")
}



load("rData/modelData.rData")
library(CLmodel)
library(Rcpp)
library(rstan)


source("R/penalisedSmooth.R")

wk <- subset(ef, Species == "Salmon" & LifeStage == "Fry" & CATCH_ID == 1)
wk <- subset(wk, T > 0)


mod1 <- test(Z ~ s(year) + s(doy), lambda = exp(c(1,1)))
mod2 <- test(Z ~ s(year) + s(doy), lambda = exp(c(5,1)))
mod3 <- test(Z ~ s(year) + s(doy), lambda = exp(c(10,1)))

newdat <- data.frame(year = 1996:2013, doy = 200)
X <- predictX(mod1, newdata = newdat)

plot(newdat $ year, X %*% coef(mod1), type = "l", col = 1)
lines(newdat $ year, X %*% coef(mod2), col = 2)
lines(newdat $ year, X %*% coef(mod3), col = 3)

newdat <- data.frame(year = 2010, doy = 150:300)
X <- predictX(mod1, newdata = newdat)

plot(newdat $ doy, X %*% coef(mod1), type = "l", col = 1)
lines(newdat $ doy, X %*% coef(mod2), col = 2)
lines(newdat $ doy, X %*% coef(mod3), col = 3)


# lets try and optimise over lambda...







# ------------------------------------------------
# 
#  Deviance residuals
# 
# ------------------------------------------------


if (Sys.info()["user"] == "millaco") {
  setwd("~/Dropbox/SarahColin/PhD/capture_prob_paper")    
  library(setwidth)
} else 
if (Sys.info()["user"] == "millarc") {
  setwd("C:/work/repos/papers/capture_prop_paper/")
} else 
if (Sys.info()["user"] == "Millarc") {
  setwd("C:/work/SMFS-report")
}



library(Rcpp)
library(rstan)


# One way of assessing the adequacy of a model is to compare it with a more general 
# model with the maximum number of parameters that can be estimated. It is referred to
# as the saturated model. In the saturated model there is basically one parameter per
# observation. The deviance assesses the goodness of fit for the model by looking at the
# difference between the log-likelihood functions of the saturated model and the model
# under investigation, i.e. (b y) (b y) sat l , − l , . Here sat b denotes the maximum likelihood
# estimator of the parameter vector of the saturated model, sat β , and b is the maximum
# likelihood estimator of the parameters of the model under investigation, β . The
# maximum likelihood estimator is the estimator that maximises the likelihood function.
# The deviance is defined as
# D = 2{l(bsat , y)− l(b, y)}.

library(devtools)
clmodel <- as.package("c:/work/repos/CLmodel")
check(clmodel)
install(clmodel)

# load data
load("rData/modelData.rData")

# compile model
m0 <- efp(X ~ 1, data = ef, passes = "Runs")

# drop any zero observations
ef <- subset(ef, T > 0)

# fit saturated model
n <- nrow(ef)
if (FALSE) {
  pest <- rep(NA, n)
  for (i in which(ef$T>0)) {
    if (i%%10==0) cat("done", i, "of", n, "     \r"); flush.console()
    pest[i] <- efp(X ~ 1, data = ef[i,], passes = "Runs", verbose = FALSE)$fitted
  }
  save(pest, file = "pest.rData")
}
load("pest.rData")

# loglik vectors
ef $ R <- with(ef, s - 1 - Z)

# keep fitted values to use as observations
ef $ yhat <- pest
ef $ llsat <- with(ef, T * log(yhat) + T * R * log(1-yhat) -T * log(1 - (1-yhat)^s) )


# fit a big model  and check overdispersion
full <- c("Trust*fyear",
         "factor(HACode)", 
         "s(Water_W, k = 6)",
         "s(Elevation_, k = 6)", 
         "s(Distance_s, k = 6)",
         "s(sinSlope, k = 6)",
         "s(Upcatch_km, k = 6)",
         "s(Urban, k = 6)",
         "s(NCTrees, k = 6)",
         "s(CTrees, k = 6)",
         "s(Mixed, k = 6)",
         "s(Marsh, k = 6)",
         "s(Other, k = 6)",
         "s(doy, k = 6)"
         )
fullf <- formula(paste("X ~", paste(full, collapse = " + ")))

# and get the individual likelihood components
mod <- efp(fullf, data = ef, passes = "Runs")
# get p predictions for ef
ef $ p <- fitted(mod)

# get logLik components
ef $ ll <- with(ef, T * log(p) + T * R * log(1-p) -T * log(1 - (1-p)^s) )

# calculate Deviance components and residuals
ef $ devcomp <- with(ef, 2 * (llsat - ll))
ef $ devres <- with(ef, sign(yhat - p)*sqrt(devcomp))

# calculate deviance 
D <- sum(ef $ devcomp)

# and check for overdispersion
npar <- mod $ rank
print(1 - pchisq(D, n-npar), 10) # so highly significant overdispersion

# scale estimate
phi <- D / (n-npar)
phi
#[1] 3.558088




