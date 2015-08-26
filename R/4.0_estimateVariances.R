
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

# ------------------------------------------------
# ------------------------------------------------
# 
#  Variance estimation
# 
# ------------------------------------------------
# ------------------------------------------------

efpackage <- devtools::as.package("c:/work/repos/faskally/ef")
#devtools::check(efpackage)
#devtools::install(efpackage)
devtools::load_all(efpackage)

# load appropriate data
load("intermediate_rData/phi.rData") # phi
load("intermediate_rData/screenedData.rData") # ef

# load final model?
finalf <-  n ~ LifeStage + Trust + fyear + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage
                   
ef <- within(ef, {
              cDistance_s = scale(Distance_s)
              cWater_W = scale(Water_W)
              cElevation_ = scale(Elevation_)
              fyear = fyear[drop = TRUE]
              })
contrasts(ef $ fyear) <- "contr.sum"
contrasts(ef $ Trust) <- "contr.treatment"
#levels(ef $ Trust) <- c("Tweed", )
ef $ LifeStage <- factor(ef $ LifeStage)

final <- efp(finalf, data = ef, pass = pass)


# bootstrap
if (FALSE) {
  nboot <- 1000
  bootb <- matrix(NA, nboot, length(coef(final)))
  ef <- ef[order(ef$sampleID),]
  sids <- paste(unique(ef $ sampleID))
  ids <- 1:length(sids)
  names(ids) <- sids
  for (i in 1:nboot) {
    cat("                                        \rrunning set ", i, "of", nboot)
    flush.console()
    bids <- sample(sids, length(sids), replace = TRUE)
    bids <- ids[bids]
    dat <- ef[1:6 + rep((bids-1)*6, each = 6),]
    while (nlevels(dat $ Trust[drop=TRUE]) != nlevels(ef $ Trust[drop=TRUE]) |
           nlevels(dat $ fyear[drop=TRUE]) != nlevels(ef $ fyear[drop=TRUE])) {
      bids <- sample(sids, length(sids), replace = TRUE)
      bids <- ids[bids]
      dat <- ef[1:6 + rep((bids-1)*6, each = 6),]    
    }
    mod <- efp(finalf, data = dat, pass = pass)
    bootb[i,] <- coef(mod)  
  }
  save(bootb, file = "intermediate_rData/bootb.rData") # ef
}
load("intermediate_rData/bootb.rData") # ef


H2 <- var(bootb)
dim(H2)
library(Matrix)
image(Matrix(cov2cor(H2)))
H2


image(Matrix(cov2cor(H1) - cov2cor(H2)))

# simple
H1 <- final $ Vb

plot(cbind(sqrt(diag(H1) * 3.58), sqrt(diag(H2))))


simb <- MASS::mvrnorm(1000, coef(final), H1)
X <- final $ Gsetup $ X

simp <- X %*% t(simb)

tmp <- 
  apply(1-simp, 2, function(p) 
       1 - tapply(p, list(ef$LifeStage, ef$sampleID), prod))
dim(tmp)


cbind(coef(final), sqrt(diag(H1)))

Cov <- cov2cor(H1)
diag(Cov) <- 0
image(Matrix(Cov))
A <- Cov[!grepl("year|Trust", rownames(Cov)), !grepl("year|Trust", rownames(Cov))]
ids <- which(apply(A, 2, function(x) any(abs(x) > 0.5)))
A <- A[ids,ids]
tmp <- round(A, 3)
tmp[lower.tri(tmp)] <- 0
Matrix(tmp)
image(Matrix(A))



# ------------------------------------------------
# ------------------------------------------------
# 
#  Model summary
# 
# ------------------------------------------------
# ------------------------------------------------
 

#  best model on all data is:
finalf <-  n ~ LifeStage + Trust + fyear + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage

forms <- list(
  LifeStage = n ~ Trust + fyear + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   s(doy, k = 3),
  Trust = n ~ LifeStage + fyear + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage,
  fyear = n ~ LifeStage + Trust + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage,
  pass = n ~ LifeStage + Trust + fyear + cWater_W + 
                   cElevation_ + cDistance_s + 
                   s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage,
  width = n ~ LifeStage + Trust + fyear + pass23 + 
                   cElevation_ + cDistance_s + 
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage,
  altitude = n ~ LifeStage + Trust + fyear + pass23 + cWater_W + 
                   cDistance_s + 
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage),
  distance = n ~ LifeStage + Trust + fyear + pass23 + cWater_W + 
                   cElevation_ + 
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage,
  lifestage.pass = n ~ LifeStage + Trust + fyear + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage,
  "s(doy).lifestage" = n ~ LifeStage + Trust + fyear + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   LifeStage:pass23 + s(doy, k = 3) + 
                   cElevation_:LifeStage,
  altitude.lifestage = n ~ LifeStage + Trust + fyear + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage))



# fit
mods <- lapply(forms, efp, data = ef, pass = pass)

m0 <- efp(finalf, data = ef, pass = pass)    
  
# reestimate phi
load("intermediate_rData/llsample.rData")

phi <- 2 * (sum(llsample) - logLik(m0)) / 2838
phi
# 3.589

tab <- cbind(dropped = dropped, summaryMods(mods, m0 = m0, order = FALSE, 
                                            fn = phiBIC, phi = 2.981, n = 2838)[,-1] ) 
tab[order(tab $ Dic, decreasing = TRUE),]
if (FALSE){
                                 dropped       ic      llik df   phi        Dic           Fp
Trust                              Trust 174647.5 -87212.43 28 2.981 551.556672 0.000000e+00
LifeStage                      LifeStage 174296.4 -86961.33 47 2.981 200.427334 0.000000e+00
fyear                              fyear 174179.2 -86946.49 36 2.981  83.284539 0.000000e+00
pass                              pass23 174149.9 -86876.17 50 2.981  53.952215 1.110223e-15
distance                     cDistance_s 174114.8 -86854.66 51 2.981  18.881634 2.377123e-07
altitude                     cElevation_ 174108.0 -86855.21 50 2.981  12.030573 9.221397e-07
width                           cWater_W 174103.3 -86848.89 51 2.981   7.352147 9.376014e-05
altitude.lifestage cElevation_:LifeStage 174102.9 -86848.72 51 2.981   7.017830 1.118159e-04
s(doy).lifestage        s(doy):LifeStage 174100.3 -86851.37 50 2.981   4.360578 4.130748e-05
lifestage.pass          LifeStage:pass23 174097.7 -86846.10 51 2.981   1.776307 1.834205e-03
}


tab2 <- cbind(dropped = dropped, summaryMods(mods, m0 = m0, order = FALSE, 
                                            fn = phiBIC, phi = 3.5896, n = 2838)[,-1] ) 
tab2[order(tab2 $ Dic, decreasing = TRUE),]

                                 dropped       ic      llik df    phi        Dic           Fp
Trust                              Trust 145074.6 -72425.97 28 3.5896 425.690069 0.000000e+00
LifeStage                      LifeStage 144808.6 -72217.44 47 3.5896 159.705658 0.000000e+00
fyear                              fyear 144696.5 -72205.12 36 3.5896  47.595544 0.000000e+00
pass                              pass23 144691.0 -72146.72 50 3.5896  42.108806 3.408385e-13
distance                     cDistance_s 144663.2 -72128.85 51 3.5896  14.332310 2.468911e-06
altitude                     cElevation_ 144656.2 -72129.31 50 3.5896   7.294784 9.634000e-06
width                           cWater_W 144653.6 -72124.07 51 3.5896   4.757594 3.701283e-04
altitude.lifestage cElevation_:LifeStage 144653.3 -72123.93 51 3.5896   4.479959 4.291206e-04
s(doy).lifestage        s(doy):LifeStage 144649.8 -72126.12 50 3.5896   0.925201 2.275501e-04
lifestage.pass          LifeStage:pass23 144649.0 -72121.75 51 3.5896   0.127112 4.513279e-03
