


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

efpackage <- devtools::as.package("c:/work/repos/faskally/ef")
#devtools::check(efpackage)
#devtools::install(efpackage)
devtools::load_all(efpackage)

# load appropriate data
load("intermediate_rData/phi.rData") # phi
load("intermediate_rData/screenedData.rData") # ef

source("R/ModelSelectionFunctions.R")



# ------------------------------------------------
# ------------------------------------------------
# 
#  Model summary table
# 
# ------------------------------------------------
# ------------------------------------------------
 
# calculate nicer covariates
ef <- within(ef, {
              cDistance_s = c(scale(Distance_s))
              cWater_W = c(scale(Water_W))
              cElevation_ = c(scale(Elevation_))
              fyear = factor(fyear)
              Trust = factor(Trust)
              LifeStage = factor(LifeStage)
              })
contrasts(ef $ fyear) <- "contr.sum"
contrasts(ef $ Trust) <- "contr.treatment"
#levels(ef $ Trust) <- c("Tweed", )


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
  "s(doy)" = n ~ LifeStage + Trust + fyear + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   LifeStage:pass23 + 
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



# fits
mods <- lapply(forms, efp, data = ef, pass = pass)

m0 <- efp(finalf, data = ef, pass = pass)    
  
# origional phi
load("intermediate_rData/phi.rData") # phi


Nsamples <- sum(tapply(ef$n, ef$sampleID, sum) > 0)

tab <- cbind(dropped = names(forms), summaryMods(mods, m0 = m0, order = FALSE, 
                                            fn = phiBIC, phi = phi, n = Nsamples)[,-1] ) 
tab[order(tab $ Dic, decreasing = TRUE),]
if (FALSE){
                              dropped       ic      llik df      phi        Dic           Fp
Trust                           Trust 174591.7 -87184.55 28 2.981916 559.121271 0.000000e+00
LifeStage                   LifeStage 174232.9 -86933.61 46 2.981916 200.337164 0.000000e+00
s(doy)                         s(doy) 174150.2 -86888.28 47 2.981916 117.633825 0.000000e+00
fyear                           fyear 174115.8 -86918.76 35 2.981916  83.190881 0.000000e+00
pass                             pass 174086.6 -86848.49 49 2.981916  53.949242 1.110223e-15
distance                     distance 174051.5 -86826.98 50 2.981916  18.880385 2.379075e-07
altitude                     altitude 174044.6 -86827.53 49 2.981916  12.028586 9.233694e-07
width                           width 174039.9 -86821.21 50 2.981916   7.336196 9.456842e-05
altitude.lifestage altitude.lifestage 174039.6 -86821.05 50 2.981916   7.018909 1.117731e-04
s(doy).lifestage     s(doy).lifestage 174036.9 -86823.69 49 2.981916   4.348364 4.157323e-05
lifestage.pass         lifestage.pass 174034.4 -86818.43 50 2.981916   1.783920 1.826999e-03
}


# reestimate phi
load("intermediate_rData/phi_new.rData") # phi

tab2 <- cbind(dropped = names(forms), summaryMods(mods, m0 = m0, order = FALSE, 
                                            fn = phiBIC, phi = phi_new, n = Nsamples)[,-1] ) 
tab2[order(tab2 $ Dic, decreasing = TRUE),]
if (FALSE) {
                              dropped       ic      llik df      phi         Dic           Fp
Trust                           Trust 145033.2 -72405.28 28 3.590581 433.3425078 0.000000e+00
LifeStage                   LifeStage 144759.5 -72196.88 46 3.590581 159.6378631 0.000000e+00
s(doy)                         s(doy) 144692.1 -72159.23 47 3.590581  92.3018985 0.000000e+00
fyear                           fyear 144647.4 -72184.55 35 3.590581  47.5246704 0.000000e+00
pass                             pass 144641.9 -72126.19 49 3.590581  42.1084324 3.410605e-13
distance                     distance 144614.2 -72108.32 50 3.590581  14.3320909 2.469640e-06
altitude                     altitude 144607.1 -72108.78 49 3.590581   7.2940432 9.640909e-06
width                           width 144604.6 -72103.53 50 3.590581   4.7448391 3.727199e-04
altitude.lifestage altitude.lifestage 144604.3 -72103.40 50 3.590581   4.4813378 4.288858e-04
s(doy).lifestage     s(doy).lifestage 144600.8 -72105.59 49 3.590581   0.9157497 2.287015e-04
lifestage.pass         lifestage.pass 144600.0 -72101.22 50 3.590581   0.1337689 4.497652e-03
}





