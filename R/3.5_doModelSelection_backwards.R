
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

# start from the top

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

fes <- c("LifeStage:pass23",
          "poly(Water_W, 1):pass23", "poly(Water_W, 1):LifeStage","s(Water_W, k=3)",
          "poly(doy, 1):pass23",
          "s(doy, k=3)", "s(doy, k = 3, by = pass23)", "s(doy, k = 3, by = LifeStage)",
          "poly(Distance_s, 1):pass23", "poly(Distance_s, 1):LifeStage","s(Distance_s, k=3)",
          "poly(Elevation_, 1):pass23", "poly(Elevation_, 1):LifeStage","s(Elevation_, k=3)",
          "Trust:pass23", "Trust:fyear",
          "s(Upcatch_km, k=3)", "s(sinSlope, k=3)", "s(Urban, k=3)", "s(woodland, k=3)",
          "s(Marsh, k=3)", "s(Other, k=3)")

last <- c("LifeStage", "pass23", "LifeStage:pass23", "Trust", "Trust:pass23", "fyear", "Trust:fyear", "s(HACode, k = 12, bs = 'gmrf', xt = list(penalty = Qhma, rank = 44))",
           "s(Elevation_, k=3)", "s(Upcatch_km, k=3)", "s(Distance_s, k=3)", "s(sinSlope, k=3)",
            "s(Water_W, k=3)", "s(Urban, k=3)", "s(woodland, k=3)", "s(Marsh, k=3)", "s(Other, k=3)", "s(doy, k=3)")

f1se <- c(f1s, fes)

start <- f1se %in% last
f1se[start]
setdiff(last, f1se[start])

out <- run1Selection(f1se, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
out
last <- f1se[out$chosen]
last



