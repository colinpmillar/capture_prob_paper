



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
ef <- subset(ef, T > 0 & Runs == 3 & Species == "Salmon")

# compile model
m0 <- efp(X ~ 1, data = ef, passes = "Runs")


# fit sitewise saturated model
#if (FALSE) {
#  n <- nrow(ef)
#  pest <- rep(NA, n)
#  for (i in 1:n) {
#    if (i%%10==0) cat("done", i, "of", n, "     \r"); flush.console()
#    pest[i] <- efp(X ~ 1, data = ef[i,], passes = "Runs", verbose = FALSE)$fitted
#  }
#  save(pest, file = "pest.rData")
#}
#load("pest.rData")

# try an alternative to estimating deviance
# this exploits that we only have 3 pass data

#f <- function(pi, Z, s = 3) Z -  pi/(1-pi) + s*pi^s/(1-pi^s)
#qest <- rep(NA, nrow(ef))
#uZs <- sort(unique(ef[["Z"]]), decreasing = TRUE)
#for (i in 1:nrow(pairs)) {
#  if (i%%100==0) cat("done", i, "of", nrow(pairs), "     \r"); flush.console()
#  id <- which(ef$Z == uZs[i])
#  qest[id] <- uniroot(f, c(0,1-1e-9), Z = 3 - 1 - uZs[i]) $ root
#}
#pest2 <- 1 - qest

# try an alternative to estimating deviance
# this exploits that we only have 3 pass data
pest <- with(ef, ifelse(X < T*(s-1)/2, 0, (3 * X - T - sqrt(T^2 + 6*X*T - 3*X^2))/(2*X)))

# calculate deviance components per observation
ef $ psat <- pest
nij <- as.matrix(dplyr::select(ef, n_R1, n_R2, n_R3))
pij <- t(sapply(ef $ psat, function(p) p * (1-p)^c(0:2)))
muij <- pij * ef$T/(1-(1-ef$psat)^3)

dij <- 2*(nij * log(nij/muij) - (nij - muij))
dij[nij == 0] <- muij[nij == 0]
dij[abs(dij) < 1e-9] <- 0
rij <- sign(nij - muij) * sqrt(dij)

plot(c(rij))


# calculate deviance per site visit and compare to chisq 1
ef $ Di <- rowSums(dij)
ef $ pchi <- pchisq(ef $ Di, 1)
plot(ef $ Di)
hist(ef $ pchi, nclass = 50) # should be uniform i think?
mean(pchisq(ef $ Di, 1) > 0.95)
mean(pchisq(ef $ Di, 1) > 0.99)

#
dord <- rev(order(ef $ Di))
pval <- sapply(1:200, function(i) 1 - pchisq(sum(ef $ Di[-dord[1:i]]), nrow(ef)-i))
plot(1:200, pval)
outliers <- dord[which(pval < 0.05)]
keep <- (1:nrow(ef))[-outliers]
ef $ outliers <- 1:nrow(ef) %in% outliers
plot(ef $ Di, col = as.numeric(ef $ outliers) + 1)

length(outliers)


par(mfrow = c(2,1))
plot(x <- seq(0, 1, length=1000), sapply(x, function(x) 1-mean(pchisq(ef $ Di, 1) > x)), 
     type = "l", ylab = "", xlab = "", sub = "observed Devs")
abline(a=0,b=1)

plot(x <- seq(0, 1, length=1000), sapply(x, function(x) 1-mean(pchisq(ef $ Di[keep], 1) > x)), 
     type = "l", ylab = "", xlab = "", sub = "Drop worse 110 - chisq(sumD, n)<0.95")
abline(a=0,b=1)

#plot(x <- seq(0, 1, length=1000), sapply(x, function(x) 1-mean(pchisq(ef $ Di[clean], 1) > x)), 
#     type = "l", ylab = "", xlab = "", sub = "drop obs with n1 <= n2")
#abline(a=0,b=1)
mtext("qq plot for uniformity of chisq(dev)>x", outer = TRUE, side = 3, line = -3)

# drop first 124 as outliers - sites where constant p assumtion is not valid.
ef2 <- subset(ef, !outliers)
save(ef2, file = "ef2.rData")


# -------------------------------------------------------------
#
#   now check against complex model for site to site variabily
#
# -------------------------------------------------------------

#    sort out new covariates
load("ef2.rData")

# woodland
ef2 $ woodland <- with(ef2, CTrees + NCTrees + Mixed)

# numbers caught on first pass
ef2 $ EventID <- as.numeric(interaction(ef2 $ Site_OBJECTID, ef2 $ Date))
df <- data.frame(totalN = tapply(ef2 $ n_R1, ef2 $ EventID, sum))
ef2 $ totalN <- df[paste(ef2 $ EventID),]

# drop out dodgy sites / samples and covariate outliers
ef2 <- subset(ef2, Site_OBJECTID != 2180)
ef2 <- subset(ef2, !(grepl("GIR", Site.Name) & year < 2006))
ef2 <- subset(ef2, !(Date == "13/09/2002" & Site_OBJECTID == 1527))
ef2 <- subset(ef2, !(Date == "13/09/2002" & Site_OBJECTID == 4004))
ef2 <- subset(ef2, !(Date == "27/09/2004" & Site_OBJECTID == 1955))
ef2 <- subset(ef2, !(Date == "03/08/2006" & Site_OBJECTID == 2626))
ef2 <- subset(ef2, keep)
ef3 <- ef2

# fit a big model  and check overdispersion
full <- c("LifeStage*Trust*fyear",
          "s(totalN, k = 4)",
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
if (FALSE) {
  bigmod <- efp(fullf, data = ef2, passes = "Runs")
  save(bigmod, ef2, file = "bigmod.rData")
}
load("bigmod.rData")


# get p predictions for ef
ef2 $ pbig <- fitted(bigmod)

# get logLik components
ef2 $ R <- with(ef2, s - 1 - Z)
ef2 $ llsat <- with(ef2, T * log(psat) + T * R * log(1-psat) - T * log(1 - (1-psat)^s) )
ef2 $ llbig <- with(ef2, T * log(pbig) + T * R * log(1-pbig) - T * log(1 - (1-pbig)^s) )

# calculate Deviance components and residuals
ef2 $ devcomp <- with(ef2, 2 * (llsat - llbig))
ef2 $ devres <- with(ef2, sign(psat - pbig)*sqrt(devcomp))
plot(ef2 $ devres)

# calculate deviance 
dord <- rev(order(abs(ef2 $ devres)))
n <- nrow(ef2)
pval <- sapply(1:1000, function(i) 1 - pchisq(sum(ef2 $ devcomp[-dord[1:i]]), n-i))
plot(1:1000, pval)
outliers <- dord[which(pval < 0.05)]
keep <- (1:n)[-outliers]

plot(ef2 $ devres, col = 2)
points(keep, ef2 $ devres[keep])

cols <- c("Site_OBJECTID","Site.Name", "Dataset", "Width", "Date", "Runs", "Area",
          "Trust", "n_R1", "n_R2", "n_R3","psat", "pbig", "devres")
head(ef2[order(ef2 $ devres, decreasing = TRUE),cols], 20)
head(ef2[order(ef2 $ devres),cols], 20)


# and check for overdispersion
npar <- bigmod $ rank
print(1 - pchisq(sum(ef2 $ devcomp), n-npar), 10) # so highly significant overdispersion

# scale estimate
sum(ef2 $ devcomp)/(nrow(ef2)-bigmod $ rank)

#[1] 2.108


# plot residuals against time and organisation
lattice::xyplot(devres ~ I(year + doy/365) | Trust, data = ef2, pch = ".")


# ------------------------------------------------
# 
#  Model fitting
# 
# ------------------------------------------------


phiBIC <- function(object, ..., phi = 2.11) {
  ll <- if ("stats4" %in% loadedNamespaces()) stats4:::logLik
    else logLik
  Nobs <- if ("stats4" %in% loadedNamespaces()) stats4:::nobs
    else nobs
  lls <- ll(object)
  nos <- attr(lls, "nobs")
  -2 * as.numeric(lls)/phi + log(nos) * attr(lls, "df")
}

# should maybe return the best model as well summaries...
summaryMods <- function(lst, m0 = NULL, order = TRUE, fn = phiBIC, phi = 2.11) {
  #aics <- sapply(lst, AIC)
  fn <- match.fun(fn)
  aics <- sapply(lst, fn)

  tab <- 
   data.frame(
    forms = sapply(lst, function(x) paste(deparse(x$formula, width.cutoff = 500L))),
    aic = aics
    )

  if (!is.null(m0)) tab $ Daic <- tab $ aic - fn(m0, phi = phi)
  if (order) tab <- tab[order(aics),]

  unique(tab)  
}

runModels <- function(chosen, f1s, data = ef3, fn = phiBIC, phi = 2.11, ...) {
  # build a formula list of additions
  formsadd <- lapply(f1s[!chosen], function(x) as.formula(paste0("X ~ ", paste(c(f1s[chosen], x), collapse = " + "))))

  descadd <- f1s[!chosen]
  dropadd <- rep("add", sum(!chosen))

  #build a formula list dropping elements
  if (any(chosen)) {
    if (sum(chosen)==1) {
      formsdrop <- list(X ~ 1)
    } else {
      formsdrop <- lapply(seq(sum(chosen)), function(i) as.formula(paste0("X ~ ", paste(f1s[chosen][-i], collapse = " + ") )))
    }
    descdrop <- f1s[chosen]
    forms <- c(formsadd, formsdrop)
    desc <- c(descadd, descdrop)
    dropadd <- c(dropadd, rep("drop", sum(chosen)))
  } else {
    forms <- formsadd
    desc <- descadd
  }

  mods <- lapply(forms, efp, data = data, passes = "Runs", ...)

  if (all(!chosen)) {
    m0 <- efp(X ~ 1, data = data, passes = "Runs", ...)
  } else {
    m0 <- efp(as.formula(paste0("X ~ ", paste(f1s[chosen], collapse = " + "))), data = data, passes = "Runs", ...)    
  }
  
  tab <- cbind(what = desc, step = dropadd, summaryMods(mods, m0 = m0, order = FALSE, fn = fn, phi = phi)[,-1] )
  tab[order(tab $ Daic),]
}
  
runSelection <- function(forms, data = ef3, fn = phiBIC, phi = 2.11) {
  chosen <- rep(FALSE, length(forms))
  out <- list(Daic = -1)
  tol <- 0
  outlist <- list() # grow this - bad!
  i <- 0
  while(out $ Daic[1] < tol & !all(chosen)) {
    out <- runModels(chosen, forms, verbose = FALSE, data = data, fn = fn, phi = phi)
    print(head(out, 10), digits = 3)
    if ( out $ Daic[1] < 0) {
      # then model is an improvement
      which <- which(forms %in% out $ what[1])
      chosen[which] <- !chosen[which]
      cat("\t", paste(forms[chosen], collapse = " + "), "\n ----------------------------------------------- \n")
      outlist[[i <- i + 1]] <- out
    } 
  }

  if (all(!chosen)) {
    mod <- efp(X ~ 1, data = data, passes = "Runs")
  } else {
    mod <- efp(as.formula(paste0("X ~ ", paste(forms[chosen], collapse = " + "))), data = data, passes = "Runs")
  }

  list(history = outlist, chosen = chosen, mod = mod)
}



## Q model for spatial regions
library(CLdata)
require(Matrix)
{
  hma <- hma[!hma $ HACode %in% c(107:108),]
  hmaadj <- spdep::poly2nb(hma, queen = FALSE)
  hmaadj <- spdep::nb2mat(hmaadj, style = "B", zero.policy = TRUE)
  # add connections for inner and outer hebs
  hmaadj[hma $ HAName == "Inner Hebrides", hma $ HAName == "Outer Hebrides"] <- 1
  hmaadj[hma $ HAName == "Outer Hebrides", hma $ HAName == "Inner Hebrides"] <- 1
  hmaadj[hma $ HAName == "Inner Hebrides", c(21, 42, 43)] <- 1
  hmaadj[c(21, 42, 43), hma $ HAName == "Inner Hebrides"] <- 1
  hmaadj <- as(hmaadj, "dgTMatrix")

  Q <- hmaadj
  if (any(Q @ x == 0)) stop("something went wrong")
  Q @ x[] <- -1/Q @ x
  diag(Q) <- rowSums(hmaadj)
  Q <- as.matrix(Q)
  colnames(Q) <- rownames(Q) <- hma $ HACode
  Qhma <- Q
}



# now some model selection
f1s <- c("LifeStage",
         "Trust",
         #"Trust:fyear",
         #"te(NEAR_X, NEAR_Y, k = c(4,4))",
         #"ti(NEAR_X, NEAR_Y, k = c(4,4))",
         #"s(HACode, k = 8, bs = 'gmrf', xt = list(penalty = Qhma))",
         "fyear",
         "s(totalN, k = 3)",
         "poly(Water_W, 1)",
         "poly(Elevation_, 1)", 
         "poly(Distance_s, 1)",
         "poly(sinSlope, 1)",
         "poly(Upcatch_km, 1)",
         "poly(Urban, 1)",
         "poly(woodland, 1)",
         "poly(Marsh, 1)",         
         "poly(Other, 1)",
         "poly(doy, 1)",
         "s(Water_W, k = 3)",
         "s(Elevation_, k = 3)", 
         "s(Distance_s, k = 3)",
         "s(sinSlope, k = 3)",
         "s(Upcatch_km, k = 3)",
         "s(Urban, k = 3)",
         "s(woodland, k = 3)",
         "s(Marsh, k = 3)",         
         "s(Other, k = 3)",
         "s(doy, k = 3)",
         "s(totalN, k = 3, by = LifeStage)",
         "poly(Water_W, 1):LifeStage",
         "poly(Elevation_, 1):LifeStage", 
         "poly(Distance_s, 1):LifeStage",
         "poly(sinSlope, 1):LifeStage",
         "poly(Upcatch_km, 1):LifeStage",
         "poly(Urban, 1):LifeStage",
         "poly(woodland, 1):LifeStage",
         "poly(Marsh, 1):LifeStage",         
         "poly(Other, 1):LifeStage",
         "poly(doy, 1):LifeStage",
         "s(Water_W, k = 3, by = LifeStage)",
         "s(Elevation_, k = 3, by = LifeStage)", 
         "s(Distance_s, k = 3, by = LifeStage)",
         "s(sinSlope, k = 3, by = LifeStage)",
         "s(Upcatch_km, k = 3, by = LifeStage)",
         "s(Urban, k = 3, by = LifeStage)",
         "s(woodland, k = 3, by = LifeStage)",
         "s(Marsh, k = 3, by = LifeStage)",         
         "s(Other, k = 3, by = LifeStage)",
         "s(doy, k = 3, by = LifeStage)"
         )


#   Loop over adding in covariates
ef3 $ LifeStage <- factor(ef3 $ LifeStage)
ef3 $ totalN <- ef3 $ totalN^.25
runSelection(f1s[29:30], data = ef3, fn = phiBIC, phi = 2.1085)

if (TRUE) {
  sel1 <- runSelection(f1s, data = ef3, fn = phiBIC, phi = 2.108524)
  save(sel1, ef3, file = "sel1.rData")
}
load("sel1.rData")

# re-estimate phi

wk <- ef3
mod <- sel1 $ mod


getScale <- function(wk, mod) {
  # get p predictions for ef
  wk $ pbig <- fitted(mod)

  # get logLik components
  wk $ R <- with(wk, s - 1 - Z)
  wk $ llsat <- with(wk, T * log(psat) + T * R * log(1-psat) - T * log(1 - (1-psat)^s) )
  wk $ llbig <- with(wk, T * log(pbig) + T * R * log(1-pbig) - T * log(1 - (1-pbig)^s) )

  # calculate Deviance components and residuals
  wk $ devcomp <- with(wk, 2 * (llsat - llbig))
  wk $ devres <- with(wk, sign(psat - pbig)*sqrt(devcomp))

  # scale estimate
  sum(wk $ devcomp)/(nrow(wk)-mod $ rank)
}

















sel1 $ mod

#   Loop over adding in covariates seperately for lifestages
#sel2 <- runSelection(f1s[-1], data = subset(ef3, LifeStage == "Fry"), fn = phiBIC)
#sel3 <- runSelection(f1s[-1], data = subset(ef3, LifeStage == "Parr"), fn = phiBIC)
#sel2[[3]]
#sel3[[3]]


#  best model on all data is:
as.formula(paste0("X ~ ", paste(f1s[sel1 $ chosen], collapse = " + ")))
#X ~ LifeStage + Trust + fyear + poly(Water_W, 1) + poly(Elevation_,1) + 
#    poly(Distance_s, 1) + poly(doy, 1) + 
#    s(totalN, k = 3, by = LifeStage) + s(doy, k = 3, by = LifeStage)


#   Drop terms for importance
finalfs <- f1s[sel1 $ chosen]
forms <- lapply(seq_along(finalfs), function(i) as.formula(paste0("X ~ ", paste(finalfs[-i], collapse = " + ") )))
mods <- lapply(forms, efp, data = ef3, passes = "Runs", verbose = FALSE)

m0 <- efp(as.formula(paste0("X ~ ", paste(finalfs, collapse = " + "))), data = ef3, passes = "Runs")    
  
tab <- cbind(dropped = finalfs, summaryMods(mods, m0 = m0, order = FALSE)[,-1] ) 
tab[order(tab $ Daic, decreasing = TRUE),]
#                           dropped      aic         Daic
#2                            Trust 220956.1 562.01555060
#1                        LifeStage 220708.0 313.94311487
#8 s(totalN, k = 3, by = LifeStage) 220554.9 160.84287716
#9    s(doy, k = 3, by = LifeStage) 220493.0  98.95337701
#3                            fyear 220468.8  74.74429456
#6              poly(Distance_s, 1) 220441.4  47.32345352
#5              poly(Elevation_, 1) 220422.7  28.67079813
#4                 poly(Water_W, 1) 220401.5   7.45802116
#7                     poly(doy, 1) 220394.1   0.02171448


# ------------------------------------------------
# 
#  Model prep
# 
# ------------------------------------------------

best <- efp(X ~ LifeStage + Trust + fyear + 
                poly(Water_W, 1) + poly(Elevation_,1) + poly(Distance_s, 1) + 
            s(totalN, k = 3, by = LifeStage) + s(doy, k = 3, by = LifeStage), 
            data = ef3, passes = "Runs", hessian = TRUE)    


#summary(best)
best

# constant p model
base <- efp(X ~ 1, data = ef3, passes = "Runs", hessian = TRUE)
#base <- efp(X ~ LifeStage, data = ef3, passes = "Runs", hessian = TRUE)


# predict p for data
# get a gam container
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


## check factors
setdiff(levels(ef3 $ Trust), levels(g1 $ model $ Trust))

X <- predict(g1, type = "lpmatrix", newdata = ef3)

ef3 $ p <- 1/(1 + exp(-X %*% coef(g1)))
ef3 $ logitp.se <- sqrt(diag(X %*% g1 $ Vp %*% t(X)))

ef3 $ offset <- with(ef3, log( (1-(1-p)^Runs) * Area) )

psim <- matrix(rnorm(1000 * nrow(ef3)), ncol = 1000, nrow = nrow(ef3))
psim <- psim * ef3 $ logitp.se + c(X %*% coef(g1))
psim <- 1/(1 + exp(-psim))
offsetsim <- with(ef3, (1-(1-psim)^Runs) * Area )
ef3 $ offset.se <- apply(offsetsim, 1, sd)
ef3 $ weights <- 1/ef3 $ offset.se^2

# add in constant p

ef3 $ constantp <- 1/(1 + exp(-coef(base)))


save(best, g1, file = "rData/bestpmodel.rData")
save(ef3, file = "rData/densmodelData.rData")


# ------------------------------------------------
# 
#  Model summaries
# 
# ------------------------------------------------

# residuals
ef4 <- ef3

# get p predictions for ef
ef4 $ pbig <- fitted(best)

# get logLik components
ef4 $ R <- with(ef4, s - 1 - Z)
ef4 $ llsat <- with(ef4, T * log(psat) + T * R * log(1-psat) - T * log(1 - (1-psat)^s) )
ef4 $ llbig <- with(ef4, T * log(pbig) + T * R * log(1-pbig) - T * log(1 - (1-pbig)^s) )

# calculate Deviance components and residuals
ef4 $ devcomp <- with(ef4, 2 * (llsat - llbig))
ef4 $ devcomp[abs(ef4 $ devcomp) < 1e-9] <- 0
ef4 $ devres <- with(ef4, sign(psat - pbig)*sqrt(devcomp))
ef4 $ devres.scaled <- ef4 $ devres / 2.11
plot(ef4 $ devres)
hist(ef4 $ devres)
qqnorm(ef4 $ devres)

# plot
library(lattice)
xyplot(devres.scaled ~ I(year + doy/365) | Trust,
       data = ef4, type = c("p", "g"), pch = 16, cex = 0.5)

xyplot(devres.scaled ~ pbig | Trust, group = grepl("GIR",Site.Name),
       data = ef4, type = c("p", "g"), pch = 16, cex = 0.5)


xyplot(devres.scaled ~ I(year + doy/365) | Site.Name,
       data = subset(ef4, Trust == "MSS"), type = c("p", "g"), pch = 16, cex = 0.5)

xyplot(devres.scaled ~ I(year + doy/365) | Site.Name, groups = year,
       data = subset(ef4, grepl("GIR", Site.Name)), type = c("p", "g"), pch = 16, cex = 0.5)


xyplot(pbig ~ psat | cut(llsat, breaks = 50),
       data = ef4, type = c("p", "g"), pch = 16, cex = 0.5)



cols <- c("Site_OBJECTID","Site.Name", "Dataset", "Width", "Date", "Runs", "Area",
          "Trust", "n_R1", "n_R2", "n_R3","psat", "pbig", "devres", "devres.scaled")

head(ef4[order(ef4 $ devres.scaled, decreasing = TRUE),cols], 20)
head(ef4[order(ef4 $ devres.scaled),cols], 20)

subset(ef4,  totalN > 700)

subset(ef4,  Trust == "Nith")

# ------------------------------------------------
# 
#  Plots of effects
# 
# ------------------------------------------------

