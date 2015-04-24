

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


