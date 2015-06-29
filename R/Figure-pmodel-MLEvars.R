


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


# load fits and model data
load("rData/bestpmodel.rData")
load("rData/densmodelData.rData")
source("R/ModelSelectionFunctions.R")
 
library(mgcv)

getScale(ef3, best)

# predict p for data
# get a gam container for best model
g1 <- gam(G = best $ Gsetup)
qr.G <- qr(best $ G)
rank.deficient <- qr.G $ pivot[abs(diag(qr.G $ qr)) < 1e-7]
whichkeep <- -rank.deficient
if (!length(whichkeep)) whichkeep <- 1:length(best $ coefficients) 
names(g1 $ coefficients[-1 * whichkeep])
g1 $ coefficients[] <- 0
g1 $ coefficients[whichkeep] <- best $ coefficients       
g1 $ Vp[] <- 0
diag(g1 $ Vp[]) <- 1e-5
g1 $ Vp[whichkeep, whichkeep] <- best $ Vb
g1 $ family <- binomial()
var.summary <- best $ Gsetup $ var.summary

# get stuff for simulating
X <- predict(g1, type = "lpmatrix", newdata = ef3)
b <- coef(g1)
Vp <- getScale(ef3, best) * g1 $ Vp

nsim <- 999
# 1. simulate ps
bsim <- MASS::mvrnorm(nsim, b, Vp)
psim <- 1/(1 + exp(-X %*% t(bsim)))
pcatch <- 1 - (1-psim)^ef3$Runs

# 2. simulate T, given N and p, T is binomial
Nhat <- round(ef3$T/pcatch)
Tsim <- Nhat # to set the dimensions
Tsim[] <- rbinom(length(Nhat), c(Nhat), c(pcatch))


# 3. compute N given T and p
Nsim <- Tsim/pcatch
Nsd <-  apply(Nsim, 1, sd)
Nest <-  rowMeans(Nsim)
Ncv <- Nsd / Nest


hist(Ncv * 100, nclass = 100)

cv.best <- Ncv
save(cv.best, file = "cv.best.rData")


#
# now using constant p model
#

library(CLmodel)
base <- efp(X ~ 1, data = ef3, passes = "Runs", hessian = TRUE)
#base <- efp(X ~ LifeStage, data = ef3, passes = "Runs", hessian = TRUE)
getScale(ef3, base)

# predict p for data
# get a gam container for base model
g1 <- gam(G = base $ Gsetup)
qr.G <- qr(base $ G)
rank.deficient <- qr.G $ pivot[abs(diag(qr.G $ qr)) < 1e-7]
whichkeep <- -rank.deficient
if (!length(whichkeep)) whichkeep <- 1:length(base $ coefficients) 
names(g1 $ coefficients[-1 * whichkeep])
g1 $ coefficients[] <- 0
g1 $ coefficients[whichkeep] <- base $ coefficients       
g1 $ Vp[] <- 0
diag(g1 $ Vp[]) <- 1e-5
g1 $ Vp[whichkeep, whichkeep] <- base $ Vb
g1 $ family <- binomial()
var.summary <- base $ Gsetup $ var.summary

# get stuff for simulating
X <- predict(g1, type = "lpmatrix", newdata = ef3)
b <- coef(g1)
Vp <- getScale(ef3, base) * g1 $ Vp

nsim <- 5000
# 1. simulate ps
bsim <- MASS::mvrnorm(nsim, b, Vp)
psim <- 1/(1 + exp(-X %*% t(bsim)))
pcatch <- 1 - (1-psim)^ef3$Runs

# 2. simulate T, given N and p, T is binomial
Nhat <- round(ef3$T/pcatch)
Tsim <- Nhat # to set the dimensions
Tsim[] <- rbinom(length(Nhat), c(Nhat), c(pcatch))

# 3. compute N given T and p
Nsim <- Tsim/pcatch
Nsd <-  apply(Nsim, 1, sd)
Nest <-  rowMeans(Nsim)
Ncv <- Nsd / Nest

hist(Ncv * 100, nclass = 50)

cv.const <- Ncv
save(cv.const, file = "cv.const.rData")


#
# now using saturated p model
#

# fit sitewise saturated model
n <- nrow(ef3)
if (FALSE) {
  psd <- rep(NA, n)
  for (i in 1:n) {
    if (i%%10==0) cat("done", i, "of", n, "     \r"); flush.console()
    mod <- efp(X ~ 1, data = ef4[i,], passes = "Runs", verbose = FALSE, hessian = TRUE)
    psd[i] <- sqrt(mod $ Vb)
  }
  save(psd, file = "psd.rData")
}
load("psd.rData")

# remove cvs where p == 0
ef4 <- ef3
ef4 $ psd <- psd
ef4 <- subset(ef4, psat > 1e-8)
# simulate densities

nsim <- 999
# 1. simulate ps
bsim <- t(sapply(1:nrow(ef4), function(i) rnorm(nsim, log(ef4 $ psat[i]/(1-ef4 $ psat[i])), ef4 $ psd[i])))
psim <- 1/(1 + exp(bsim))
pcatch <- 1 - (1-psim)^ef4$Runs

## NOTE: restric pcatch to be between zero and 1
pcatch[pcatch < 1e-9] <- 1e-9

# 2. simulate T, given N and p, T is binomial
Nhat <- round(ef4$T/pcatch)
Tsim <- Nhat # to set the dimensions
Tsim[] <- rbinom(length(Nhat), c(Nhat), c(pcatch))

# 3. compute N given T and p
Nsim <- Tsim/pcatch
Nsd <-  apply(Nsim, 1, sd)
Nest <-  rowMeans(Nsim)
Ncv <- Nsd / Nest

hist(Ncv * 100, nclass = 50)

cv.sat <- Ncv
save(cv.sat, file = "cv.sat.rData")





# ----------------------------------
#
# combine plots
#
# ----------------------------------


load("cv.const.rData")
load("cv.sat.rData")
load("cv.best.rData")



{
png(file = "figures/pmodel_mlevars.png", width = 7, height = 5, units = "in", res = 500)

axisfun <- function()
  axis(2, line = -1.8, at = seq(0, 0.12, by = 0.02), labels = c("0.00", "", "0.04", "", "0.08", "", "0.12"), las = 1)

xlim = c(0,200)
breaks <- seq(0, 300, by = 1)

heights <- c(2.5, 2.5, 1) * 0.0005
heights2 <- c(2.5, 2.5, 1.4)
layout(matrix(c(1,2,3), ncol = 1), heights = heights2)
par(mar = c(1,0.5,2,1), oma = c(5,1,0,0))

hist(cv.const * 100, breaks = breaks[1:(sum(breaks < max(cv.const)*100)+1)], col = "lightblue", xlim = xlim, ylim = c(0,heights[1])*100, xlab = "", axes = FALSE, ylab = "", main = "", freq = FALSE)
axisfun()
mtext(side = 3, "a", adj = 0, font = 2)

hist(cv.best * 100, breaks = breaks[1:(sum(breaks < max(cv.best)*100)+1)], col = "lightblue", xlim = xlim, ylim = c(0,heights[2])*100, xlab = "", axes = FALSE, ylab = "", main = "", freq = FALSE)
axisfun()
mtext(side = 3, "b", adj = 0, font = 2)

hist(cv.sat[cv.sat < 3] * 100, breaks = breaks, col = "lightblue", xlim = xlim, ylim = c(0,heights[3])*100, xlab = "", axes = FALSE, ylab = "", main = "", freq = FALSE)
axisfun()
mtext(side = 3, "c", adj = 0, font = 2)

axis(1, line = 1)
mtext(side = 1, "Coefficient of variation", outer = TRUE, line = 2)

dev.off()
}
