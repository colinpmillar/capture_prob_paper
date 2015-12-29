

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


# all main effects
fcs <- c("LifeStage",
         "Trust",
         "fyear",
         "pass23",
         "s(HACode, k = 12, bs = 'gmrf', xt = list(penalty = Qhma, rank = 44))")

fcints <- c("LifeStage:pass23")

fl <- function(what) paste0("poly(", what, ", 1)")
fs <- function(what) paste0("s(", what, ", k=3)")
fsint <- function(what1, what2 = c("pass23", "LifeStage")) paste0("s(", what1, ", k=3, by=", what2,")")
flint <- function(what1, what2 = c("pass23", "LifeStage")) paste0("poly(", what1, ", 1):", what2)

cvars <-c("Water_W","Elevation_", "Distance_s", "sinSlope",
         "Upcatch_km", "Urban", "woodland", "Marsh", "Other", "doy")

f1s <- sapply(cvars, fl, USE.NAMES = FALSE)
f2s <- sapply(cvars, fs, USE.NAMES = FALSE)
f3s <- as.vector(sapply(cvars, fsint, USE.NAMES = FALSE))
f4s <- as.vector(sapply(cvars, flint, USE.NAMES = FALSE))





# fit big model
last <- c("LifeStage", "pass23", "LifeStage:pass23", "Trust", "Trust:pass23", "fyear", "Trust:fyear", "s(HACode, k = 12, bs = 'gmrf', xt = list(penalty = Qhma, rank = 44))",
          "s(Elevation_, k=3)", "s(Upcatch_km, k=3)", "s(Distance_s, k=3)", "s(sinSlope, k=3)",
          "s(Water_W, k=3)", "s(Urban, k=3)", "s(woodland, k=3)", "s(Marsh, k=3)", "s(Other, k=3)", "s(doy, k=3)")
m0 <- runModels(rep(TRUE, length(last)), last, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392)
fullf <- apply(!sapply(m0$what, function(x) x == last), 2, function(x) paste(last[x], collapse = " + "))
m0$fullf <- fullf
m0$id="m0"

# try and drop Trust:year
last <- c("LifeStage", "pass23", "LifeStage:pass23", "Trust", "Trust:pass23", "fyear", "s(HACode, k = 12, bs = 'gmrf', xt = list(penalty = Qhma, rank = 44))",
          "s(Elevation_, k=3)", "s(Upcatch_km, k=3)", "s(Distance_s, k=3)", "s(sinSlope, k=3)",
          "s(Water_W, k=3)", "s(Urban, k=3)", "s(woodland, k=3)", "s(Marsh, k=3)", "s(Other, k=3)", "s(doy, k=3)")
m1 <- runModels(rep(TRUE, length(last)), last, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392)
fullf <- apply(!sapply(m1$what, function(x) x == last), 2, function(x) paste(last[x], collapse = " + "))
m1$fullf <- fullf
m1$id="m1"

# Trust:pass23
last <- c("LifeStage", "pass23", "LifeStage:pass23", "Trust", "fyear", "s(HACode, k = 12, bs = 'gmrf', xt = list(penalty = Qhma, rank = 44))",
          "s(Elevation_, k=3)", "s(Upcatch_km, k=3)", "s(Distance_s, k=3)", "s(sinSlope, k=3)",
          "s(Water_W, k=3)", "s(Urban, k=3)", "s(woodland, k=3)", "s(Marsh, k=3)", "s(Other, k=3)", "s(doy, k=3)")
m2 <- runModels(rep(TRUE, length(last)), last, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392)
fullf <- apply(!sapply(m2$what, function(x) x == last), 2, function(x) paste(last[x], collapse = " + "))
m2$fullf <- fullf
m2$id="m2"

# now we can do forwards backwards
form <- c(fcs, fcints, f1s, f2s, f3s)
start <- form %in% last

form[start]
setdiff(last, form[start])
res <- list()
i <- 0

res[[i <- i+1]] <- run1Selection(form, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
print(i); res[[i]]$tab

#[1] 1
#what step       ic      llik df   phi           Dic           Fp
#35 s(HACode, k = 12, bs = 'gmrf', xt = list(penalty = Qhma, rank = 44)) drop 174251.5 -86839.01 62 2.981 -43.270300035 1.893487e-08
#38                                                   s(Elevation_, k=3) drop 174277.0 -86810.18 71 2.981 -17.706086024 6.731862e-01
#40                                                     s(sinSlope, k=3) drop 174281.9 -86812.61 71 2.981 -12.844253687 5.925586e-02
#41                                                   s(Upcatch_km, k=3) drop 174286.1 -86814.70 71 2.981  -8.665924880 7.346814e-03
#44                                                        s(Marsh, k=3) drop 174288.3 -86811.22 72 2.981  -6.380461399 9.036874e-02
#43                                                     s(woodland, k=3) drop 174289.8 -86811.92 72 2.981  -4.961505568 3.842323e-02
#45                                                        s(Other, k=3) drop 174291.2 -86812.63 72 2.981  -3.560121540 1.709273e-02
#30                                            s(doy, k=3, by=LifeStage)  add 174294.3 -86800.34 75 2.981  -0.390711709 7.983792e-05


# loop this?
res[[i <- i+1]] <- run1Selection(form, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = res[[i]]$chosen)
res[[i]]$tab


# until i == 11?
# add in linear interactions
form_ext <- c(fcs, fcints, f1s, f2s, f3s, f4s)
start <- form_ext %in% form[res[[11]]$chosen]

# add in this (linear distance from sea):
start[9] <- TRUE

form_ext[start]
setdiff(form[res[[11]]$chosen], form_ext[start])
res2 <- list()
i <- 0

# run first
res2[[i <- i+1]] <- run1Selection(form_ext, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = start)
res2[[i]]$tab

# loop this?
i
res2[[i <- i+1]] <- run1Selection(form_ext, data = ef, fn = phiBIC, phi = 2.981, nobs = 10392, start = res2[[i]]$chosen)
res2[[i]]$tab

# now table all the fits and add in the model

#save(res, res2, form, form_ext, m0, m1, m2, file = "intermediate_rData/fits.rData")
#load("intermediate_rData/fits.rData")

# columns: Model, DF, BICadj.

# make full model column
for (i in 1:length(res)) {
  fullf <- apply((sapply(res[[i]]$tab$what, function(x) x == form) + res[[i]]$chosen) == 1, 2, function(x) paste(form[x], collapse = " + "))
  res[[i]]$tab$fullf <- fullf
  res[[i]]$tab$id <- paste0("res_",i)
}
for (i in 1:length(res2)) {
  fullf <- apply((sapply(res2[[i]]$tab$what, function(x) x == form_ext) + res2[[i]]$chosen) == 1, 2, function(x) paste(form_ext[x], collapse = " + "))
  res2[[i]]$tab$fullf <- fullf
  res2[[i]]$tab$id <- paste0("res2_",i)
}


# collate info
tab_orig <- do.call(rbind, c(lapply(res, "[[", "tab"), lapply(res2, "[[", "tab"), list(m0, m1, m2)))
str(tab_orig)

tab <- tab_orig[c("fullf", "df", "ic")]
tab <- tab[order(tab$ic, decreasing = FALSE),]
rownames(tab) <- NULL

# simplify text
tab$fullf <- gsub(", k=3", "", tab$fullf)
tab$fullf <- gsub("s[(]HACode, k = 12, bs = [']gmrf['], xt = list[(]penalty = Qhma, rank = 44[)][)]", "hydrometric area", tab$fullf)
tab$fullf <- gsub("poly[(]", "", tab$fullf)
tab$fullf <- gsub(", 1[)]", "", tab$fullf)
tab$fullf <- gsub(", by=LifeStage[)]", "):LifeStage", tab$fullf)


tab$fullf <- gsub("LifeStage", "life-stage", tab$fullf)
tab$fullf <- gsub("Trust", "organisation", tab$fullf)
tab$fullf <- gsub("fyear", "year", tab$fullf)
tab$fullf <- gsub("Water_W", "width", tab$fullf)
tab$fullf <- gsub("Elevation_", "altitude", tab$fullf)
tab$fullf <- gsub("pass23", "pass", tab$fullf)
tab$fullf <- gsub("Distance_s", "distance to sea", tab$fullf)
tab$fullf <- gsub("sinSlope", "gradient", tab$fullf)
tab$fullf <- gsub("Upcatch_km", "upstream catchment area", tab$fullf)
tab$fullf <- gsub("doy", "day of year", tab$fullf)
tab$fullf <- gsub("Urban", "urban", tab$fullf)
tab$fullf <- gsub("Marsh", "marsh", tab$fullf)
tab$fullf <- gsub("Other", "other", tab$fullf)

names(tab) <- c("Model", "d.f.", "BIC")
tab$BIC <- tab$BIC - min(tab$BIC)
tab$BIC <- sprintf("%.2f", tab$BIC)

write.csv(unique(tab), row.names = FALSE, file = "table_of_fits.csv")

head(tab, 20)
