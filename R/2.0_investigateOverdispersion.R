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

# One way of assessing the adequacy of a model is to compare it with a more general 
# model with the maximum number of parameters that can be estimated. It is referred to
# as the saturated model. In the saturated model there is basically one parameter per
# observation. The deviance assesses the goodness of fit for the model by looking at the
# difference between the log-likelihood functions of the saturated model and the model
# under investigation, i.e. (b y) (b y) sat l , − l , . Here sat b denotes the maximum likelihood
# estimator of the parameter vector of the saturated model, sat , and b is the maximum
# likelihood estimator of the parameters of the model under investigation . The
# maximum likelihood estimator is the estimator that maximises the likelihood function.
# The deviance is defined as
# D = 2{l(bsat, y) − l(b, y)}.


efpackage <- devtools::as.package("c:/work/repos/faskally/ef")
#devtools::check(efpackage)
#devtools::install(efpackage)
devtools::load_all(efpackage)

# -------------------------------------------------------------
#
#   now check against complex model for site to site variabily
#
# -------------------------------------------------------------

#    sort out new covariates
load("intermediate_rData/screenedData.rData")

## Q model for spatial regions
library(CLdata)
library(gmrf)
if (FALSE) {
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

  save(Qhma, file = "intermediate_rData/Qhma.rData")
}
load("intermediate_rData/Qhma.rData")

# quick test
m1 <- efp(n ~ s(HACode, k = 12, bs = 'gmrf', xt = list(penalty = Qhma)), 
          data = ef, pass = pass, verbose = TRUE)

# fit a big model  and check overdispersion
big <- c("LifeStage",
         "Trust",
         "fyear",
         "pass23",
         "Trust:pass23",
         "LifeStage:pass23",
         "Trust:fyear",
         "s(HACode, k = 12, bs = 'gmrf', xt = list(penalty = Qhma))", 
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
bigf <- formula(paste("n ~", paste(big, collapse = " + ")))

# and get the individual ps
if (FALSE) {
  bigmod <- efp(bigf, data = ef, pass = pass, verbose = TRUE, hessian = FALSE)
  save(bigmod, file = "intermediate_rData/bigmod.rData")
}
load("intermediate_rData/bigmod.rData")
ef $ bigfit <- bigmod $ fitted

# now estimate a sample effect for each sample
if (FALSE) {
  n <- nrow(ef) # data points
  N <- n / 6    # site visits (samples)
  llsample <- rep(NA, N)
  sIDs <- sort(unique(ef $ sampleID))
  for (i in 1:N) {
    if (i%%10 == 0) cat("done", i, "of", N, "     \r"); flush.console()
    samp <- subset(ef, sampleID == sIDs[i])
    offset <- logit(samp $ bigfit)  
    mod <- efp(n ~ 1, offset = offset, data = samp, pass = pass, hessian = FALSE)
    llsample[i] <- logLik(mod)
  }
  cat("Done!/n")
  save(llsample, file = "intermediate_rData/llsample.rData")
}
load("intermediate_rData/llsample.rData")

load("intermediate_rData/llsat.rData")
load("intermediate_rData/llscreen.rData")

# compare deviance
n <- nrow(ef) # data points
N <- n / 6    # site visits (samples)
m <- bigmod $ df.null - bigmod $ df.residual
## within sample overdispersion
2 * (sum(llsat) - sum(llsample)) / (3*N - m)
# 1.389

phi <- 2 * (sum(llsample) - logLik(bigmod)) / N
phi
# 2.981 - slightly larger than last time

# and check for overdispersion
npar <- bigmod $ rank
print(1 - pchisq(sum(ef2 $ devcomp), n-npar), 10) # so highly significant overdispersion

# scale estimate
sum(ef2 $ devcomp)/(nrow(ef2)-bigmod $ rank)

#[1] 2.193646
#[1] 2.198613? newer estimate...

save(phi, file = "intermediate_rData/phi.rData")

