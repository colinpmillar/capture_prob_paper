


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
              fyear = fyear[drop = TRUE]
              })
contrasts(ef $ fyear) <- "contr.sum"
contrasts(ef $ Trust) <- "contr.treatment"
#levels(ef $ Trust) <- c("Tweed", )
ef $ LifeStage <- factor(ef $ LifeStage)



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
  
# reestimate phi
load("intermediate_rData/llsample.rData")

phi <- 2 * (sum(llsample) - logLik(m0)) / 2838
phi
# 3.589

tab <- cbind(dropped = names(forms), summaryMods(mods, m0 = m0, order = FALSE, 
                                            fn = phiBIC, phi = 2.981, n = 2838)[,-1] ) 
tab[order(tab $ Dic, decreasing = TRUE),]
if (FALSE){
                              dropped       ic      llik df   phi        Dic           Fp
Trust                           Trust 174647.5 -87212.43 28 2.981 551.556672 0.000000e+00
LifeStage                   LifeStage 174296.4 -86961.33 47 2.981 200.427334 0.000000e+00
s(doy)                         s(doy) 174213.6 -86915.97 48 2.981 117.647025 0.000000e+00
fyear                           fyear 174179.2 -86946.49 36 2.981  83.284539 0.000000e+00
pass                             pass 174149.9 -86876.17 50 2.981  53.952215 1.110223e-15
distance                     distance 174114.8 -86854.66 51 2.981  18.881634 2.377123e-07
altitude                     altitude 174108.0 -86855.21 50 2.981  12.030573 9.221397e-07
width                           width 174103.3 -86848.89 51 2.981   7.352147 9.376014e-05
altitude.lifestage altitude.lifestage 174102.9 -86848.72 51 2.981   7.017830 1.118159e-04
s(doy).lifestage     s(doy).lifestage 174100.3 -86851.37 50 2.981   4.360578 4.130748e-05
lifestage.pass         lifestage.pass 174097.7 -86846.10 51 2.981   1.776307 1.834205e-03
}


tab2 <- cbind(dropped = names(forms), summaryMods(mods, m0 = m0, order = FALSE, 
                                            fn = phiBIC, phi = 3.5896, n = 2838)[,-1] ) 
tab2[order(tab2 $ Dic, decreasing = TRUE),]
if (FALSE) {
                              dropped       ic      llik df    phi        Dic           Fp
Trust                           Trust 145074.6 -72425.97 28 3.5896 425.690069 0.000000e+00
LifeStage                   LifeStage 144808.6 -72217.44 47 3.5896 159.705658 0.000000e+00
s(doy)                         s(doy) 144741.2 -72179.77 48 3.5896  92.308397 0.000000e+00
fyear                           fyear 144696.5 -72205.12 36 3.5896  47.595544 0.000000e+00
pass                             pass 144691.0 -72146.72 50 3.5896  42.108806 3.408385e-13
distance                     distance 144663.2 -72128.85 51 3.5896  14.332310 2.468911e-06
altitude                     altitude 144656.2 -72129.31 50 3.5896   7.294784 9.634000e-06
width                           width 144653.6 -72124.07 51 3.5896   4.757594 3.701283e-04
altitude.lifestage altitude.lifestage 144653.3 -72123.93 51 3.5896   4.479959 4.291206e-04
s(doy).lifestage     s(doy).lifestage 144649.8 -72126.12 50 3.5896   0.925201 2.275501e-04
lifestage.pass         lifestage.pass 144649.0 -72121.75 51 3.5896   0.127112 4.513279e-03
}





# ------------------------------------------------
# ------------------------------------------------
# 
#  Model summary plots
# 
# ------------------------------------------------
# ------------------------------------------------

source("R/ModelPlottingFunctions.R")

final <- efp(finalf, data = ef, pass = pass)


# refit final model
m0 <- efp(finalf, data = ef, pass = pass)    
  

# predict p for data
# get a gam container
g1 <- gam(G = m0 $ Gsetup)
qr.G <- qr(m0 $ G)
rank.deficient <- qr.G $ pivot[abs(diag(qr.G $ qr)) < 1e-7]
whichkeep <- -rank.deficient
if (!length(whichkeep)) whichkeep <- 1:length(m0 $ coefficients) 
names(g1 $ coefficients[-1 * whichkeep])
g1 $ coefficients[] <- 0
g1 $ coefficients[whichkeep] <- m0 $ coefficients       
g1 $ Vp[] <- 0
diag(g1 $ Vp[]) <- 1e-5
g1 $ Vp[whichkeep, whichkeep] <- m0 $ Vb
g1 $ family <- binomial()


# #################################
#
# setup prediction data
#
# #################################

var.summary <- m0 $ Gsetup $ var.summary
var.type <- sapply(var.summary, function(x) is(x)[1])
var.names <- names(var.summary)

which <- c("LifeStage", "Trust", "fyear", "pass23", "cWater_W", "cElevation_", 
           "cDistance_s", "doy")
fullnames <- data.frame(names = c("Life-stage", 
                          "Organisation", 
                          "Year",
                          "Pass",
                          "Channel width (m)", 
                          "Altitude (m)", 
                          "Distance to sea (km)", 
                          "Day of year"),
                        stringsAsFactors = FALSE)
rownames(fullnames) <- which

var.fullnames <- fullnames[var.names,]

# set base prediction levels
pdata0 <- ifelse(var.type == "numeric", 
                   lapply(var.summary, "[", 2), 
                   lapply(var.summary, function(x) levels(x)[floor(nlevels(x)/2)])
                )
pdata0 $ LifeStage <- "Fry"
pdata0 $ Trust <- "MSS"
pdata0 $ fyear <- "2006"


# set up prediction ranges
pdata1 <- ifelse(var.type == "numeric", 
                   lapply(var.summary, function(x) if (length(x) == 3) seq(x[1], x[3], length=100) else 0), 
                   lapply(var.summary, function(x) levels(x))
                )


# do plots
{
png(file = "figures/pmodel_grid.png", width = 7, height = 7, units = "in", res = 400)

ylim <- c(0.35, 0.82)

par(mar = c(5,3,3,0.5)) # c(bottom, left, top, right)

layout(rbind(c(1,1,2), c(3,4,4), c(5,6,7)))


#   Trust predictions of p
# ----------------------------------------
pdata <- getPlotData("Trust", bootb = bootb)
pdata <- pdata[order(pdata $ p),]
pdata $ Trust <- as.character(pdata $ Trust)
pdata $ Trust[pdata $ Trust == "Other"] <- "Caithness"
factorPlot(pdata, xlab = "Organisation", ylim = ylim, labcex = 0.8)


#   Lifestage predictions of p --- and pass?
# ----------------------------------------
pdata <- getPlotData(c("LifeStage", "pass23"), bootb = bootb)
pdata <- pdata[c(1,3,2,4),]
factorPlot(pdata, xlab = fullnames["LifeStage",], ylim = ylim, labcex = 0.8, yaxislab = FALSE)
text(1:nrow(pdata), ylim[1] - diff(ylim)*.08,  paste(pdata $ LifeStage, ": pass", pdata $ pass23), srt=45, xpd = TRUE, adj = 1, cex = 0.8)


#   DoY predictions of p
# ----------------------------------------
pdata <- getPlotData(c("doy", "LifeStage"), bootb = bootb)
continuousPlot2(pdata, xlab = "Life-stage x day of year", rug = ef $ doy, ylim = ylim, yaxislab = TRUE)


#   Year predictions of p
# ----------------------------------------
pdata <- getPlotData("fyear", bootb = bootb)
factorPlot(pdata, xlab = "Year", ylim = ylim, yaxislab = FALSE)


#   Altitide predictions of p
# ----------------------------------------
pdata <- getPlotData(c("cElevation_", "LifeStage"), bootb = bootb)
pdata $ x <- pdata $ cElevation_ * sd(ef$Elevation_) + mean(ef$Elevation_)
continuousPlot2(pdata, xlab = "Life-stage x Altitude (m)", rug = ef $ Elevation_, ylim = ylim, yaxislab = TRUE)


#   DS predictions of p
# ----------------------------------------
pdata <- getPlotData("cDistance_s", bootb = bootb)
pdata $ x <- pdata $ cDistance_s * sd(ef$Distance_s) + mean(ef$Distance_s)
continuousPlot(pdata, xlab = fullnames["cDistance_s",], rug = ef $ Distance_s, ylim = ylim, yaxislab = TRUE)


#   Width predictions of p
# ----------------------------------------
pdata <- getPlotData("cWater_W", bootb = bootb)
pdata $ x <- pdata $ cWater_W * sd(ef$Water_W) + mean(ef$Water_W)
continuousPlot(pdata, xlab = fullnames["cWater_W",], rug = ef $ Water_W, ylim = ylim, yaxislab = TRUE)


dev.off()

}
