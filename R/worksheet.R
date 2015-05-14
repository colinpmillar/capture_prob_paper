



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
  setwd("C:/work/repos/papers/capture_prop_paper/")
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
  setwd("C:/work/repos/papers/capture_prop_paper/")
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
#check(clmodel)
#install(clmodel)
load_all(clmodel)

# load data
load("rData/modelData.rData")
# drop any zero observations and non 3 pass fishings
ef <- subset(ef, T > 0 & Runs == 3 & Species == "Salmon" & LifeStage == "Fry")

# compile model
m0 <- efp(X ~ 1, data = ef, passes = "Runs")


# fit sitewise saturated model
n <- nrow(ef)
if (FALSE) {
  pest <- rep(NA, n)
  for (i in 1:n) {
    if (i%%10==0) cat("done", i, "of", n, "     \r"); flush.console()
    pest[i] <- efp(X ~ 1, data = ef[i,], passes = "Runs", verbose = FALSE)$fitted
  }
  save(pest, file = "pest.rData")
}
load("pest.rData")

# calculate deviance components per observation
ef $ p <- pest
nij <- as.matrix(dplyr::select(ef, n_R1, n_R2, n_R3))
pij <- t(sapply(ef $ p, function(p) p * (1-p)^c(0:2)))
muij <- pij * ef$T/(1-(1-ef$p)^3)

dij <- 2*(nij * log(nij/muij) - (nij - muij))
dij[nij == 0] <- muij[nij == 0]
dij[abs(dij) < 1e-9] <- 0
rij <- sign(nij - muij) * sqrt(dij)




# calculate deviance per site visit and compare to chisq 1
ef $ Di <- rowSums(dij)
ef $ pchi <- pchisq(ef $ Di, 1)
plot(ef $ Di)
hist(ef $ pchi, nclass = 50) # should be uniform i think?
mean(pchisq(ef $ Di, 1) > 0.95)
mean(pchisq(ef $ Di, 1) > 0.99)

#
dord <- rev(order(ef $ Di))
pval <- sapply(1:200, function(i) 1 - pchisq(sum(Di[-dord[1:i]]), n-i))
plot(1:200, pval)
outliers <- dord[which(pval < 0.05)]
keep <- (1:n)[-outliers]

par(mfrow = c(2,2))
plot(x <- seq(0, 1, length=1000), sapply(x, function(x) 1-mean(pchisq(ef $ Di, 1) > x)), 
     type = "l", ylab = "", xlab = "", sub = "observed Devs")
abline(a=0,b=1)

plot(x <- seq(0, 1, length=1000), sapply(x, function(x) 1-mean(pchisq(ef $ Di[keep], 1) > x)), 
     type = "l", ylab = "", xlab = "", sub = "Drop worse 110 - chisq(sumD, n)<0.95")
abline(a=0,b=1)

plot(x <- seq(0, 1, length=1000), sapply(x, function(x) 1-mean(pchisq(ef $ Di[clean], 1) > x)), 
     type = "l", ylab = "", xlab = "", sub = "drop obs with n1 <= n2")
abline(a=0,b=1)
mtext("qq plot for uniformity of chisq(dev)>x", outer = TRUE, side = 3, line = -3)

# drop first 110 as outliers - sites where constant p assumtion is not valid.
ef2 <- ef[keep,]

# now check against complex model for site to site variabily
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

# and get the individual ps
mod <- efp(fullf, data = ef2, passes = "Runs")
# get p predictions for ef
ef2 $ pbig <- fitted(mod)

# get logLik components
ef2 $ R <- with(ef2, s - 1 - Z)
ef2 $ llsat <- with(ef2, T * log(p) + T * R * log(1-p) - T * log(1 - (1-p)^s) )
ef2 $ llbig <- with(ef2, T * log(pbig) + T * R * log(1-pbig) - T * log(1 - (1-pbig)^s) )

# calculate Deviance components and residuals
ef2 $ devcomp <- with(ef2, 2 * (llsat - llbig))
ef2 $ devres <- with(ef2, sign(p - pbig)*sqrt(devcomp))
plot(ef2 $ devres)

# calculate deviance 
dord <- rev(order(abs(ef2 $ devres)))
pval <- sapply(1:1000, function(i) 1 - pchisq(sum(ef2 $ devcomp[-dord[1:i]]), n-i))
plot(1:1000, pval)
outliers <- dord[which(pval < 0.05)]
keep <- (1:n)[-outliers]

plot(ef2 $ devres)
points(keep, ef2 $ devres[keep], pch = 16)


# and check for overdispersion
npar <- mod $ rank
print(1 - pchisq(D, n-npar), 10) # so highly significant overdispersion

# scale estimate
sum(ef2 $ devcomp)/(nrow(ef2)-mod $ rank)

#[1] 3.07

# drop massive residual
drop <- which(abs(ef2 $ devres) > 10)
sum(ef2 $ devcomp[-drop])/(nrow(ef2[-drop,])-mod $ rank)
ef2[drop,]

plot(ef2 $ devres[-drop])


