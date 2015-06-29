




# ------------------------------------------------
# 
#  Full analysis
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
# keep salmon, 3 pass fishings within the covariate bounds
ef <- subset(ef, Runs == 3 & Species == "Salmon" & keep)
#ef <- subset(ef, Species == "Salmon" & keep)

# summaries full data
sampledf <- data.frame(
   siteid = ef $ Site_OBJECTID,
   visitid = paste0(ef $ Site_OBJECTID, ":", ef $ Date),
   sampleid = 1:nrow(ef))

sapply(sampledf, function(x) length(unique(x)))

unique(table(sampledf $ visitid))
table(table(sampledf $ visitid))

# double sampling sometimes - different trusts sampling same site on same day
#oddvisits <- names(which(table(sampledf $ visitid) == 4))
#ef[sampledf $ visitid %in% oddvisits[1],]

# drop any zero observations
ef <- ef[!is.na(ef $ T), ]
ef <- ef[ef $ T > 0, ]



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
ef $ psat <- replace(pest, pest == 0, 1e-9)
ef $ psat <- replace(ef $ psat, pest == 1, 1-1e-9)

# calculate deviance components per observation
# need to adjust p estimates that are zero to get a deviance for them...
# the problem is that you can have an estimate of zero p, while there are still counts 
nij <- as.matrix(dplyr::select(ef, n_R1, n_R2, n_R3))
pij <- t(sapply(replace(pest, pest == 0, 1e-9), function(p) p * (1-p)^c(0:2)))
muij <- pij * ef$T/(1-(1-replace(pest, pest == 0, 1e-9))^3)

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
# 111

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

# some covariate summaries
length(unique(ef3 $ Trust))
length(unique(ef3 $ HACode))

range(apply(table(ef $ CATCH_ID, ef $ HACode), 2, function(x) sum(x>0)))

length(unique(ef3 $ year))

# fit a big model  and check overdispersion
full <- c("LifeStage*Trust*fyear",
          "factor(HACode)", 
         "s(Water_W, k = 3)",
         "s(Elevation_, k = 3)", 
         "s(Distance_s, k = 3)",
         "s(sinSlope, k = 3)",
         "s(Upcatch_km, k = 3)",
         "s(Urban, k = 3)",
         "s(woodland, k = 3)",
         "s(Marsh, k = 3)",         
         "s(Other, k = 3)",
         "s(doy, k = 3)"
         )
fullf <- formula(paste("X ~", paste(full, collapse = " + ")))

# and get the individual ps
if (FALSE) {
  bigmod <- efp(fullf, data = ef3, passes = "Runs")
  save(bigmod, ef3, file = "bigmod.rData")
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

#[1] 2.193646
#[1] 2.198613? newer estimate...

# plot residuals against time and organisation
lattice::xyplot(devres ~ I(year + doy/365) | Trust, data = ef2, pch = 16, cex = 0.6, grid = TRUE)

# is there correlation accross lifestage

ef2 $ sitevisit <- as.numeric(factor(paste0(ef2 $ Site_OBJECTID, ef2$Date)))
wk <- subset(ef2, sitevisit %in% which(table(ef2 $ sitevisit) == 2))
wk2 <- by(wk, wk$sitevisit, function(x) c(x$devres[x$LifeStage == "Fry"], x$devres[x$LifeStage == "Parr"]))
wk3 <- t(matrix(unlist(wk2), nrow = 2))
head(wk3)
plot(wk3, xlab = "Fry residuals", ylab = "Parr residuals", xlim = range(wk3), ylim = range(wk3))
abline(h = 0, v = 0, col = grey(0.5))
cov(wk3)

# boostrap distribution of covariance
cor.sim <- function(...) {
  ids <- sample(1:nrow(wk3), replace = TRUE)
  cor(wk3[ids,])[1,2]
}

scor <- sapply(1:1000, cor.sim)
hist(scor)
cis <- 2*mean(scor) - quantile(scor, c(0.975, 0.025))
# 0.1910985 0.2741948 
# so there is some correlation between lifestage...

# ------------------------------------------------
# ------------------------------------------------
# 
#  Model fitting
# 
# ------------------------------------------------
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

library(devtools)
clmodel <- as.package("c:/work/repos/CLmodel")
#check(clmodel)
#install(clmodel)
load_all(clmodel)



load("bigmod.rData")
source("R/ModelSelectionFunctions.R")
 
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
{
f1s <- c("LifeStage",
         "Trust",
         "fyear",
         "Trust:fyear",
         "Trust:LifeStage",
         "fyear:LifeStage",
         "Trust:fyear:LifeStage",
         "factor(HACode)",      
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
         "s(NEAR_X, k = 4)",
         "s(NEAR_Y, k = 4)",
         "te(NEAR_X, NEAR_Y, k = c(4,4))",
         "ti(NEAR_X, NEAR_Y, k = c(4,4))",
         "s(HACode, k = 12, bs = 'gmrf', xt = list(penalty = Qhma))",
         "s(HACode, k = 16, bs = 'gmrf', xt = list(penalty = Qhma))",
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
         "s(doy, k = 3, by = LifeStage)",
         "s(totalN, k = 3)",
         "s(totalN, k = 3, by = LifeStage)"
         )[1:53]
}

#   Loop over adding in covariates
ef3 $ LifeStage <- factor(ef3 $ LifeStage)
#ef3 $ totalN <- ef3 $ totalN^.25
runSelection(f1s[29:30], data = ef3, fn = phiBIC, phi = 2.1085)

start <- rep(FALSE, 53)
start[1:19] <- TRUE
fullf <- formula(paste("X ~", paste(f1s[start], collapse = " + ")))
bigmodtest <- efp(fullf, data = ef3, passes = "Runs")

logLik(bigmod)
logLik(bigmodtest)

# cheat!
phiBIC(bigmod)
start[7] <- FALSE
fullf <- formula(paste("X ~", paste(f1s[start], collapse = " + ")))
modnext <- efp(fullf, data = ef3, passes = "Runs")
phiBIC(modnext) - phiBIC(bigmod)

start[4] <- FALSE
fullf <- formula(paste("X ~", paste(f1s[start], collapse = " + ")))
modnext1 <- efp(fullf, data = ef3, passes = "Runs")
phiBIC(modnext1) - phiBIC(modnext)

start[5] <- FALSE
fullf <- formula(paste("X ~", paste(f1s[start], collapse = " + ")))
modnext2 <- efp(fullf, data = ef3, passes = "Runs")
phiBIC(modnext2) - phiBIC(modnext1)

start[6] <- FALSE
fullf <- formula(paste("X ~", paste(f1s[start], collapse = " + ")))
modnext3 <- efp(fullf, data = ef3, passes = "Runs")
phiBIC(modnext3) - phiBIC(modnext2)


{
f1s <- c("LifeStage",
         "Trust",
         "fyear",
         "factor(HACode)",      
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
         "te(NEAR_X, NEAR_Y, k = c(4,4))",
         "s(HACode, k = 12, bs = 'gmrf', xt = list(penalty = Qhma))",
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
}

start <- rep(FALSE, length(f1s))
start[1:24] <- TRUE
logLik(bigmod)
getScale(ef3, bigmod)
if (FALSE) {
  sel1 <- runSelection(f1s, data = ef3, fn = phiBIC, phi = 2.20, start = start)
  save(sel1, ef3, file = "sel1.rData")
}
load("sel1.rData")

if (FALSE) {
  getScale(ef3, bigmod)
  getScale(ef3, sel1 $ mod)
  start <- sel1 $ chosen
  start[c(14,23)] <- FALSE
  sel2 <- runSelection(f1s, data = ef3, fn = phiBIC, phi = 2.20, start = start)


  phiBIC(sel2 $ mod, phi = 2.19)
  getScale(ef3, sel2 $ mod)
  #[1] 2.430373
  start <- sel2 $ chosen
  start[28] <- FALSE
  start[16] <- TRUE
  sel3 <- runSelection(f1s, data = ef3, fn = phiBIC, phi = 2.43, start = start)
  getScale(ef3, sel3 $ mod)
  # [1] 2.433985

  save(sel3, ef3, file = "sel3.rData")
}
load("sel3.rData")



#  best model on all data is:
as.formula(paste0("X ~ ", paste(f1s[sel3 $ chosen], collapse = " + ")))
#X ~ LifeStage + Trust + fyear + poly(Water_W, 1) + poly(Elevation_, 1) +
#    poly(Distance_s, 1) + s(doy, k = 3, by = LifeStage)

#   Drop terms for importance
#  choose terms to drop
finalfs <- f1s[sel3 $ chosen]
forms <- lapply(seq_along(finalfs), function(i) as.formula(paste0("X ~ ", paste(finalfs[-i], collapse = " + ") )))
dropped <- finalfs

# remove Lifestage completely
forms[[1]] <- as.formula(paste0("X ~ ", paste(c(finalfs[-c(1,7)], f1s[14]), collapse = " + ")))
dropped[1] <- "LifeStage"

# drop only lifestage interactions
forms[[8]] <- as.formula(paste0("X ~ ", paste(c(finalfs[-7], f1s[14]), collapse = " + ")))
dropped[8] <- "s(doy):LifeStage"



# fit
mods <- lapply(forms, efp, data = ef3, passes = "Runs", verbose = FALSE)

m0 <- efp(as.formula(paste0("X ~ ", paste(finalfs, collapse = " + "))), data = ef3, passes = "Runs")    
  
tab <- cbind(dropped = dropped, summaryMods(mods, m0 = m0, order = FALSE, fn = phiBIC, phi = 2.434)[,-1] ) 
tab[order(tab $ Daic, decreasing = TRUE),]
if (FALSE){
                        dropped      aic      llik df   phi       Daic           Fp
2                         Trust 192200.6 -95993.78 25 2.434 599.121108 0.000000e+00
1                     LifeStage 191868.5 -95738.26 46 2.434 267.030974 0.000000e+00
7 s(doy, k = 3, by = LifeStage) 191714.5 -95665.49 45 2.434 112.972028 0.000000e+00
3                         fyear 191637.4 -95678.10 33 2.434  35.931798 0.000000e+00
6           poly(Distance_s, 1) 191631.8 -95611.36 48 2.434  30.271694 5.099405e-10
8              s(doy):LifeStage 191612.9 -95606.17 47 2.434  11.378164 7.013889e-07
5           poly(Elevation_, 1) 191610.1 -95600.54 48 2.434   8.637276 3.495075e-05
4              poly(Water_W, 1) 191607.8 -95599.37 48 2.434   6.294991 1.199766e-04
}


# ------------------------------------------------
# 
#  Model prep
# 
# ------------------------------------------------

best <- efp(X ~ LifeStage + Trust + fyear + poly(Water_W, 1) + poly(Elevation_, 1) +
            poly(Distance_s, 1) + s(doy, k = 3, by = LifeStage), 
            data = ef3, passes = "Runs", hessian = TRUE)    


#summary(best)
best
getScale(ef3, best)
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

