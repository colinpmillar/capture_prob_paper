


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

library(mgcv)

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
Vp <- g1 $ Vp

# simulate densities
bsim <- MASS::mvrnorm(1000, b, Vp)
psim <- 1/(1 + exp(-X %*% t(bsim)))
Lsim <- 1/( (1-(1-psim)^ef3$Runs) * ef3$Area)
Lse <- apply(Lsim, 1, sd)
Lcv <- apply(Lsim, 1, sd) / apply(Lsim, 1, mean)

plot(density(Lse))
hist(Lcv, nclass = 100)

Lcv.best <- Lcv

#
# now using constant p model
#

library(CLmodel)
base <- efp(X ~ 1, data = ef3, passes = "Runs", hessian = TRUE)
#base <- efp(X ~ LifeStage, data = ef3, passes = "Runs", hessian = TRUE)

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
Vp <- g1 $ Vp

# simulate densities
bsim <- MASS::mvrnorm(1000, b, Vp)
psim <- 1/(1 + exp(-X %*% t(bsim)))
Lsim <- 1/( (1-(1-psim)^ef3$Runs) * ef3$Area)
Lse <- apply(Lsim, 1, sd)
Lcv <- apply(Lsim, 1, sd) / apply(Lsim, 1, mean)

plot(density(Lse))
hist(Lcv, nclass = 100)

Lcv.base <- Lcv


#
# now using saturated p model
#

# fit sitewise saturated model
n <- nrow(ef3)
if (FALSE) {
  psd <- rep(NA, n)
  for (i in 1:n) {
    if (i%%10==0) cat("done", i, "of", n, "     \r"); flush.console()
    mod <- efp(X ~ 1, data = ef3[i,], passes = "Runs", verbose = FALSE, hessian = TRUE)
    psd[i] <- sqrt(mod $ Vb)
  }
  save(psd, file = "psd.rData")
}
load("psd.rData")


# simulate densities
bsim <- t(sapply(1:nrow(ef3), function(i) rnorm(1000, ef3 $ psat[i], psd[i])))
psim <- 1/(1 + exp(-bsim))
Lsim <- 1/ ((1-(1-psim)^ef3$Runs) * ef3$Area)
Lse <- apply(Lsim, 1, sd)
Lcv <- apply(Lsim, 1, sd) / apply(Lsim, 1, mean)

plot(density(Lse))
hist(Lcv, nclass = 100)

Lcv.sat <- Lcv




# ----------------------------------
#
# combine plots
#
# ----------------------------------

png(file = "figures/pmodel_mlevars.png", width = 7, height = 4, units = "in", res = 500)

par(mfrow = c(1,3))
hist(Lcv.best * 100, nclass = 1000, col = "lightblue", xlim = c(0,5))
hist(Lcv.sat * 100, nclass = 1000, col = "lightblue", xlim = c(0,100))
hist(Lcv.base * 100, nclass = 1000, col = "lightblue", xlim = c(0, 0.5))
abline(v = Lcv.base[1]*100, col = "red")

dev.off()

