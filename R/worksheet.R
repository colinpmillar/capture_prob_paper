

if (Sys.info()["user"] == "millaco") {
  setwd("~/Dropbox/SarahColin/PhD/capture_prob_paper")    
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


logLik(test(exp(1)))
logLik(test(exp(2)))
logLik(test(exp(3)))
logLik(test(exp(4)))

coef(test(exp(1)))

X <- test(1) $ G

plot(wk $ year, X %*% coef(test(exp(1))), pch = 16, col = 1)

points(wk $ year, X %*% coef(test(exp(2))), pch = 16, col = 2)
points(wk $ year, X %*% coef(test(exp(3))), pch = 16, col = 3)
points(wk $ year, X %*% coef(test(exp(15))), pch = 16, col = 4)

