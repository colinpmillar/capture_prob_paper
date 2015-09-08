

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
#  Simulation of density errors
# 
# ------------------------------------------------
# ------------------------------------------------
 
