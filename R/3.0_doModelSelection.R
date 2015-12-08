
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
#  Model selection
#
# ------------------------------------------------
# ------------------------------------------------

#efpackage <- devtools::as.package("c:/work/repos/faskally/ef")
#devtools::check(efpackage)
#devtools::install(efpackage)
#devtools::load_all(efpackage)
if (FALSE) {
  httr::set_config(
    httr::use_proxy(url="192.168.41.8", port=80)
  )
  devtools::install_github("faskally/ef@v1.0")
}
library(ef)


library(gmrf)
source("R/ModelSelectionFunctions.R")


# load appropriate data
load("intermediate_rData/phi.rData") # phi

load("intermediate_rData/bigmod.rData") # bigmod

load("intermediate_rData/Qhma.rData") # Qhma
load("intermediate_rData/screenedData.rData") # ef
ef $ LifeStage <- factor(ef $ LifeStage)

# start from bottom

f1s <- c("LifeStage",
         "Trust",
         "fyear",
         "pass23",
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
         "s(HACode, k = 12, bs = 'gmrf', xt = list(penalty = Qhma, rank = 44))")

# up the penalty, nobs should be number of non zero LifeStage:Sample sets times 2.
nobs <- sum(with(ef, tapply(n, list(LifeStage, sampleID), sum)) > 0) * 2
nobs
# 10392

# null model BIC
phiBIC(efp(n ~ 1, data = ef, pass = pass), phi = 2.981, nobs = 10392)
#175361

# step 1
out <- run1Selection(f1s, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392)
out

# step 2
out <- run1Selection(f1s, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = out $ chosen)
out
f1s[out$chosen]

# step 3
f1se <- c(f1s, "LifeStage:Trust")
start <- c(out $ chosen, FALSE)
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
f1se <- f1se[out$chosen]

# step 4
f1se <- c(f1s, "LifeStage:Trust", "LifeStage:fyear", "Trust:fyear")
start <- c(out $ chosen, rep(FALSE, 2))
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
f1s[out$chosen]

# step 5
f1se <- c(f1s, "LifeStage:Trust", "LifeStage:fyear", "Trust:fyear",
               "LifeStage:pass23", "Trust:pass23", "fyear:pass23")
start <- rep(FALSE, length(f1se))
start[which(out$chosen)] <- TRUE
f1s[out$chosen]
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
f1s[out$chosen]


# step 6
f1se <- c(f1s, "LifeStage:Trust",
               "LifeStage:pass23", "Trust:pass23", "fyear:pass23",
               "poly(Water_W, 1):pass23", "poly(Water_W, 1):LifeStage","s(Water_W, k=3)")
start <- rep(FALSE, length(f1se))
start[which(out$chosen)] <- TRUE
f1se[start]
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
f1se[out$chosen]


# step 7
f1se <- c(f1s, "LifeStage:pass23",
               "poly(Water_W, 1):pass23", "poly(Water_W, 1):LifeStage","s(Water_W, k=3)",
               "poly(doy, 1):pass23", "poly(doy, 1):LifeStage","s(doy, k=3)")
start <- rep(FALSE, length(f1se))
start[which(out$chosen)] <- TRUE
f1se[start]
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
last <- f1se[out$chosen]
last

# step 8
f1se <- c(f1s, "LifeStage:pass23",
               "poly(Water_W, 1):pass23", "poly(Water_W, 1):LifeStage","s(Water_W, k=3)",
               "poly(doy, 1):pass23", "poly(doy, 1):LifeStage",
               "s(doy, k=3)", "s(doy, k = 3, by = pass23)", "s(doy, k = 3, by = LifeStage)")
start <- f1se %in% last
f1se[start]
# remove linear doy term as smooth is in now
start[f1se == "poly(doy, 1)"] <- FALSE
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
last <- f1se[out$chosen]
last


# step 9
f1se <- c(f1s, "LifeStage:pass23",
               "poly(Water_W, 1):pass23", "poly(Water_W, 1):LifeStage","s(Water_W, k=3)",
               "poly(doy, 1):pass23", "poly(doy, 1):LifeStage",
               "s(doy, k=3)", "s(doy, k = 3, by = pass23)", "s(doy, k = 3, by = LifeStage)",
               "poly(Distance_s, 1):pass23", "poly(Distance_s, 1):LifeStage","s(Distance_s, k=3)")
start <- f1se %in% last
f1se[start]
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
last <- f1se[out$chosen]
last


# step 10
f1se <- c(f1s, "LifeStage:pass23",
               "poly(Water_W, 1):pass23", "poly(Water_W, 1):LifeStage","s(Water_W, k=3)",
               "poly(doy, 1):pass23", "poly(doy, 1):LifeStage",
               "s(doy, k=3)", "s(doy, k = 3, by = pass23)", "s(doy, k = 3, by = LifeStage)",
               "poly(Distance_s, 1):pass23", "poly(Distance_s, 1):LifeStage","s(Distance_s, k=3)")
start <- f1se %in% last
f1se[start]
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
last <- f1se[out$chosen]
last

# step 11
f1se <- c(f1s, "LifeStage:pass23",
               "poly(Water_W, 1):pass23", "poly(Water_W, 1):LifeStage","s(Water_W, k=3)",
               "poly(doy, 1):pass23", "poly(doy, 1):LifeStage",
               "s(doy, k=3)", "s(doy, k = 3, by = pass23)", "s(doy, k = 3, by = LifeStage)",
               "poly(Distance_s, 1):pass23", "poly(Distance_s, 1):LifeStage","s(Distance_s, k=3)",
               "poly(Elevation_, 1):pass23", "poly(Elevation_, 1):LifeStage","s(Elevation_, k=3)")
start <- f1se %in% last
f1se[start]
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
last <- f1se[out$chosen]
last


# step 12
f1se <- c(f1s, "LifeStage:pass23",
               "poly(Water_W, 1):pass23", "poly(Water_W, 1):LifeStage","s(Water_W, k=3)",
               "poly(doy, 1):pass23", "poly(doy, 1):LifeStage",
               "s(doy, k=3)", "s(doy, k = 3, by = pass23)", "s(doy, k = 3, by = LifeStage)",
               "poly(Distance_s, 1):pass23", "poly(Distance_s, 1):LifeStage","s(Distance_s, k=3)",
               "poly(Elevation_, 1):pass23", "poly(Elevation_, 1):LifeStage","s(Elevation_, k=3)")
start <- f1se %in% last
f1se[start]
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
last <- f1se[out$chosen]
last

# step 13
f1se <- c(f1s, "LifeStage:pass23",
               "poly(Water_W, 1):pass23", "poly(Water_W, 1):LifeStage","s(Water_W, k=3)",
               "poly(doy, 1):pass23",
               "s(doy, k=3)", "s(doy, k = 3, by = pass23)", "s(doy, k = 3, by = LifeStage)",
               "poly(Distance_s, 1):pass23", "poly(Distance_s, 1):LifeStage","s(Distance_s, k=3)",
               "poly(Elevation_, 1):pass23", "poly(Elevation_, 1):LifeStage","s(Elevation_, k=3)")
start <- f1se %in% last
f1se[start]
# remove linear doy term as smooth is in now
start[f1se %in% c("s(doy, k=3)","poly(Elevation_, 1)")] <- FALSE
f1se[start]
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
last <- f1se[out$chosen]
last

# step 14
f1se <- c(f1s, "LifeStage:pass23",
               "poly(Water_W, 1):pass23", "poly(Water_W, 1):LifeStage","s(Water_W, k=3)",
               "poly(doy, 1):pass23",
               "s(doy, k=3)", "s(doy, k = 3, by = pass23)", "s(doy, k = 3, by = LifeStage)",
               "poly(Distance_s, 1):pass23", "poly(Distance_s, 1):LifeStage","s(Distance_s, k=3)",
               "poly(Elevation_, 1):pass23", "poly(Elevation_, 1):LifeStage","s(Elevation_, k=3)")
start <- f1se %in% last
f1se[start]
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
last <- f1se[out$chosen]
last

# step 16
ef $ logGradient <- log(tan(ef $ Slope_deg / 180 * pi) + 0.002)
ef $ clogGradient <- scale(ef $ logGradient)
f1se <- c(f1s, "LifeStage:pass23", "clogGradient",
               "poly(Water_W, 1):pass23", "poly(Water_W, 1):LifeStage","s(Water_W, k=3)",
               "poly(doy, 1):pass23",
               "s(doy, k=3)", "s(doy, k = 3, by = pass23)", "s(doy, k = 3, by = LifeStage)",
               "poly(Distance_s, 1):pass23", "poly(Distance_s, 1):LifeStage","s(Distance_s, k=3)",
               "poly(Elevation_, 1):pass23", "poly(Elevation_, 1):LifeStage","s(Elevation_, k=3)")
start <- f1se %in% last
f1se[start]
out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
last <- f1se[out$chosen]
last



# finished!!

as.formula(paste("n ~", paste(last, collapse = " + ")))
final <-  n ~ LifeStage + Trust + fyear + pass23 + poly(Water_W, 1) +
              poly(Elevation_, 1) + poly(Distance_s, 1) +
              LifeStage:pass23 + s(doy, k = 3, by = LifeStage) +
              poly(Elevation_, 1):LifeStage
mod1 <- efp(final, data = ef, pass = pass)

finalclean <-  n ~ LifeStage + Trust + fyear + pass23 + cWater_W +
                   cElevation_ + cDistance_s +
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage) +
                   cElevation_:LifeStage

ef2 <- within(ef, {
              cDistance_s = scale(Distance_s)
              cWater_W = scale(Water_W)
              cElevation_ = scale(Elevation_)
              })

mod2 <- efp(finalclean, data = ef2, pass = pass)

logLik(mod1)
logLik(mod2)



passtest <-  n ~ LifeStage + Trust + fyear + pass + cWater_W +
                   cElevation_ + cDistance_s +
                   LifeStage:pass + s(doy, k = 3, by = LifeStage) +
                   cElevation_:LifeStage

mod3 <- efp(passtest, data = ef2, pass = pass)

logLik(mod1)
logLik(mod2)
logLik(mod3)







