
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



# load appropriate data
load("intermediate_rData/phi.rData") # phi
load("intermediate_rData/screenedData.rData") # ef

source("R/ModelSelectionFunctions.R")





# ------------------------------------------------
# ------------------------------------------------
#
#  Model summary plot panel
#
# ------------------------------------------------
# ------------------------------------------------

source("R/ModelPlottingFunctions.R")

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
contrasts(ef $ Trust) <- "contr.sum"
#levels(ef $ Trust) <- c("Tweed", )


#  best model on all data is:
finalf <-  n ~ LifeStage + Trust + fyear + pass23 + cWater_W +
                   cElevation_ + cDistance_s +
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage) +
                   cElevation_:LifeStage

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


# load bootstraps for CIs
load("intermediate_rData/bootb.rData") # ef

# do plots
{
#png(file = "figures/pmodel_grid.png", width = 7, height = 7, units = "in", res = 400)
png(file = "C:/work/Dropbox/CaptureProbPaper/resubmission/Figure2.png", width = 7, height = 7, units = "in", res = 400)

ylim <- c(0.35, 0.82)

par(mar = c(5,3,3,0.5), oma = c(0, 3, 0, 0)) # c(bottom, left, top, right)

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
factorPlot(pdata, xlab = paste(fullnames["LifeStage",], ": pass"),
           ylim = ylim, labcex = 0.8, yaxislab = FALSE)
text(1:nrow(pdata), ylim[1] - diff(ylim)*.08,  paste(pdata $ LifeStage, ": pass", pdata $ pass23), srt=45, xpd = TRUE, adj = 1, cex = 0.8)

#   DoY predictions of p
# ----------------------------------------
pdata <- getPlotData(c("doy", "LifeStage"), bootb = bootb)
continuousPlot2(pdata, xlab = "Life-stage : day of year", rug = ef $ doy, ylim = ylim, yaxislab = TRUE)


#   Year predictions of p
# ----------------------------------------
pdata <- getPlotData("fyear", bootb = bootb)
factorPlot(pdata, xlab = "Year", ylim = ylim, yaxislab = FALSE)


#   Altitide predictions of p
# ----------------------------------------
pdata <- getPlotData(c("cElevation_", "LifeStage"), bootb = bootb)
pdata $ x <- pdata $ cElevation_ * sd(ef$Elevation_) + mean(ef$Elevation_)
continuousPlot2(pdata, xlab = "Life-stage : altitude (m)", rug = ef $ Elevation_, ylim = ylim, yaxislab = TRUE)


#   DS predictions of p
# ----------------------------------------
pdata <- getPlotData("cDistance_s", bootb = bootb)
pdata $ x <- pdata $ cDistance_s * sd(ef$Distance_s) + mean(ef$Distance_s)
continuousPlot(pdata, xlab = fullnames["cDistance_s",], rug = ef $ Distance_s, ylim = ylim, yaxislab = FALSE)


#   Width predictions of p
# ----------------------------------------
pdata <- getPlotData("cWater_W", bootb = bootb)
pdata $ x <- pdata $ cWater_W * sd(ef$Water_W) + mean(ef$Water_W)
continuousPlot(pdata, xlab = fullnames["cWater_W",], rug = ef $ Water_W, ylim = ylim, yaxislab = FALSE)

# x axis label
mtext("Capture probability", side = 2, outer = TRUE, line = 1.5, font = 1, cex = 1)

# plot lettering
mtext(c("a", "b"), side = 3, outer = TRUE, font = 2, cex = 1,
      line = -2, at = c(0, 2)/3 + 0.05)
mtext(c("c", "d"), side = 3, outer = TRUE, font = 2, cex = 1,
      line = -20, at = c(0, 1)/3 + 0.05)
mtext(c("e", "f", "g"), side = 3, outer = TRUE, font = 2, cex = 1,
      line = -37.75, at = c(0, 1, 2)/3 + 0.05)


dev.off()

}



# ------------------------------------------------
# ------------------------------------------------
#
# summary table for parr effect
#
# ------------------------------------------------
# ------------------------------------------------

pdata <- getPlotData(c("LifeStage", "pass23"), bootb = bootb)
pdata <- pdata[c(1,3,2,4),]
pdata[,c("LifeStage", "pass23", "p", "cil", "ciu")]
if (FALSE) {
  LifeStage pass23         p       cil       ciu
1       Fry      1 0.6749494 0.6684000 0.6814081
3       Fry      2 0.6536586 0.6426738 0.6645524
2      Parr      1 0.7194048 0.7084989 0.7301849
4      Parr      2 0.6757401 0.6573938 0.6932413
}

withinVars <- c("LifeStageParr", "pass232", "LifeStageParr:pass232")
coef(m0)[withinVars]
#        LifeStageParr               pass232 LifeStageParr:pass232
#           0.29888890           -0.09549676           -0.11175150
sqrt(diag(m0 $ Vb[withinVars, withinVars]))
#        LifeStageParr               pass232 LifeStageParr:pass232
#           0.01349444            0.01082034            0.02077458
cov2cor(m0 $ Vb[withinVars, withinVars])
#                      LifeStageParr    pass232 LifeStageParr:pass232
#LifeStageParr             1.0000000 -0.2153085             0.3340658
#pass232                  -0.2153085  1.0000000            -0.5052329
#LifeStageParr:pass232     0.3340658 -0.5052329             1.0000000


# ------------------------------------------------
# ------------------------------------------------
#
# spatial plots of effects
#
# ------------------------------------------------
# ------------------------------------------------


coast <- rgdal::readOGR("mapdata", "britisles")[2,]


# set up prediction data frame
getPlotData <- function(var, func = exp, model = g1) {
  args <- pdata0[!names(pdata0) %in% var]
  args[var] <- pdata1[var]

  pdata <- do.call(expand.grid, args)

  # and predict
  pdata[c("fit", "se")] <- predict(model, newdata = pdata, se.fit = TRUE)
  pdata $ est <- func(pdata $ fit)
  pdata $ cil <- func(pdata $ fit - 2*pdata $ se)
  pdata $ ciu <- func(pdata $ fit + 2*pdata $ se)
  pdata $ var <- paste(var, collapse = ":")

  if (length(var) > 1) {
    pdata[paste(var, collapse = ":")] <- do.call(interaction, c(pdata[var], list(sep = " ")))
  }

  pdata $ x <- as.numeric(pdata[[paste(var, collapse = ":")]])

  pdata
}

# general map plot

mapplot <- function(data, z = "fit", CATCH_ID = NULL, rivers = FALSE, ctm = TRUE,
                    xlim = c(5540, 414105), ylim = c(530223, 1033753), main = "",
                    breaks = NULL, cols = cols, ...) {

  data[c("fit", "se")] <- predict(g1, newdata = data, se = TRUE)
  data $ fit <- data $ fit - mean(data $ fit)
  print(range(data $ fit))

  if (is.null(breaks)) {
    breaks <- quantile(data[[z]], 0:10/10)
    breaks <- unique(breaks)
    breaks[1] <- -Inf
    breaks[length(breaks)] <- Inf
  }
  nbreaks <- length(breaks)
  if (is.null(cols)) {
    cols <- colorRampPalette(c("red", "gold", "green", "blue"))(nbreaks)
  }

  data $ colgrp <- as.numeric(cut(data[[z]], breaks = breaks))
  data <- unique(data[c("NEAR_X", "NEAR_Y", "colgrp")])

  col <- cols[data $ colgrp]

  if (is.null(CATCH_ID)) CATCH_ID <- redctm $ CATCH_ID

  plot(coast, border = grey(0.7), xlim = xlim, ylim = ylim, main = main)
  if (ctm) {
    plot(redctm[redctm $ CATCH_ID %in% CATCH_ID,], border = grey(0.8), add = TRUE)
  }
  if (rivers) {
    plot(redrivs[redrivs $ CATCH_ID %in% CATCH_ID,], add = TRUE, col = "lightblue")
  }
  with(data, {
    #points(NEAR_X, NEAR_Y, bg = paste0(col, "AA"), col = col, pch = 21, ...)
    points(NEAR_X, NEAR_Y, bg = paste0(col, "AA"), col = col, pch = 1, ...)
  })
}



# quick maps
library(sp)
library(gplots)
library(CLdata)

data(redctm)
data(hma)
redctm <- redctm[redctm $ CATCH_ID %in% ef $ CATCH_ID,]
hma <- hma[hma $ HACode %in% ef $ HACode,]

coast @ bbox <- bbox(redctm)

{
#png("figures/covarMaps.png", width = 7, height = 9, res = 600, units = "in")
png("C:/work/Dropbox/CaptureProbPaper/resubmission/Figure3.png", width = 7, height = 9, res = 600, units = "in")




par(mfrow = c(3,3), mar = c(0.5,0,0,0)) # c(bottom, left, top, right)
layout(rbind(c(1,2,3), c(4,5,6), c(7,8,9)), widths = c(1,1,0.5))

cex <- 0.8

# over all prediction

alldat <- ef[c(names(pdata0), "NEAR_X", "NEAR_Y")]
alldat $ year <- as.numeric(pdata0 $ fyear)
alldat $ fyear <- pdata0 $ fyear
alldat $ doy <- pdata0 $ doy
alldat <- unique(alldat)

alldat $ fit <- predict(g1, newdata = alldat, type = "response")

nbreaks <- 15
breaks <- quantile(alldat $ fit, seq(0, 1, length = nbreaks + 1))

breaks[1] <- 0.35
breaks[length(breaks)] <- 0.8

breaks <- seq(0.35, 0.8, by = 0.025)
nbreaks <- length(breaks) - 1

cols <- rev(rich.colors(nbreaks))

alldat $ colgrp <- as.numeric(cut(alldat $ fit, breaks = breaks))
col <- cols[alldat $ colgrp]

plot(coast, border = grey(0.7))
plot(redctm, border = grey(0.8), add = TRUE)
with(alldat, {
  points(NEAR_X, NEAR_Y, col = col, pch = 16, cex = cex)
  points(NEAR_X, NEAR_Y, col = col, pch = 1, cex = cex)
})
text(190000, 570000, "Fitted values", font = 2, adj = 1, cex = 0.9)


## plot for Trust
alldat <- ef[c(names(pdata0), "NEAR_X", "NEAR_Y")]
keep <- c("Trust")

for(i in setdiff(names(pdata0), keep)) {
  alldat[[i]] <- pdata0[[i]]
}
alldat $ fit <- predict(g1, newdata = alldat, type = "response")

alldat $ colgrp <- as.numeric(cut(alldat $ fit, breaks = breaks))
col <- cols[alldat $ colgrp]

plot(coast, border = grey(0.7))
plot(redctm, border = grey(0.8), add = TRUE)
with(alldat, {
  points(NEAR_X, NEAR_Y, col = col, pch = 16, cex = cex)
  points(NEAR_X, NEAR_Y, col = col, pch = 1, cex = cex)
})
text(190000, 570000, "Organisation", font = 2, adj = 1, cex = 0.9)


# now for legend
plot(0,0,xlim = c(0,1), ylim = c(0,1), ann = FALSE, type = "n", axes = FALSE)

x0 <- 0.35; y0 <- 1; dx <- 0.2; dy <- 1/nbreaks
val <- breaks
val <- sprintf("%.2f", rev(val))
val[1:((length(val)-1)/2)*2] <- ""
for (i in 1:nbreaks-1) {
polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - i*dy + c(0, -dy)[c(1,1,2,2)], col = rev(cols)[i+1])
}
text(x0 + 1.1*dx, y0 - 0:nbreaks*dy, val, font = 2, cex = 1, adj = 0)
mtext("Capture probability", side = 2, line = -10)



# new scale for these 4

breaks <- seq(0.59, 0.7, by = 0.01)
nbreaks <- length(breaks) - 1
cols <- rev(rich.colors(nbreaks))


## plot for DS
alldat <- ef[c(names(pdata0), "NEAR_X", "NEAR_Y")]
keep <- c("cDistance_s")

for(i in setdiff(names(pdata0), keep)) {
  alldat[[i]] <- pdata0[[i]]
}
alldat $ fit <- predict(g1, newdata = alldat, type = "response")

alldat $ colgrp <- as.numeric(cut(alldat $ fit, breaks = breaks))
col <- cols[alldat $ colgrp]

plot(coast, border = grey(0.7))
plot(redctm, border = grey(0.8), add = TRUE)
with(alldat, {
  points(NEAR_X, NEAR_Y, col = col, pch = 16, cex = cex)
  points(NEAR_X, NEAR_Y, col = col, pch = 1, cex = cex)
})
text(190000, 570000, "Distance to sea", font = 2, adj = 1, cex = 0.9)



## plot for Altitude
alldat <- ef[c(names(pdata0), "NEAR_X", "NEAR_Y")]
keep <- c("cElevation_")

for(i in setdiff(names(pdata0), keep)) {
  alldat[[i]] <- pdata0[[i]]
}
alldat $ fit <- predict(g1, newdata = alldat, type = "response")

alldat $ colgrp <- as.numeric(cut(alldat $ fit, breaks = breaks))
col <- cols[alldat $ colgrp]

plot(coast, border = grey(0.7))
plot(redctm, border = grey(0.8), add = TRUE)
with(alldat, {
  points(NEAR_X, NEAR_Y, col = col, pch = 16, cex = cex)
  points(NEAR_X, NEAR_Y, col = col, pch = 1, cex = cex)
})
text(190000, 570000, "Altitude", font = 2, adj = 1, cex = 0.9)


# now for legend
plot(0,0,xlim = c(0,1), ylim = c(0,1), ann = FALSE, type = "n", axes = FALSE)

x0 <- 0.35; y0 <- 1; dx <- 0.2; dy <- 1/nbreaks
val <- breaks
val <- sprintf("%.2f", rev(val))
val[1:((length(val)-1)/2)*2] <- ""
val[length(val)] <- ""
for (i in 1:nbreaks-1) {
  polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - i*dy + c(0, -dy)[c(1,1,2,2)], col = rev(cols)[i+1])
}
text(x0 + 1.1*dx, y0 - 0:nbreaks*dy, val, font = 2, cex = 1, adj = 0)
mtext("Capture probability", side = 2, line = -10)



## plot for Large scale spatial

alldat <- ef[c(names(pdata0), "NEAR_X", "NEAR_Y")]
keep <- c("cElevation_", "cDistance_s")

for(i in setdiff(names(pdata0), keep)) {
  alldat[[i]] <- pdata0[[i]]
}
alldat $ fit <- predict(g1, newdata = alldat, type = "response")

alldat $ colgrp <- as.numeric(cut(alldat $ fit, breaks = breaks))
col <- cols[alldat $ colgrp]

plot(coast, border = grey(0.7))
plot(redctm, border = grey(0.8), add = TRUE)
with(alldat, {
  points(NEAR_X, NEAR_Y, col = col, pch = 16, cex = cex)
  points(NEAR_X, NEAR_Y, col = col, pch = 1, cex = cex)
})
text(190000, 570000, "Distance to sea,\n altitude", font = 2, adj = 1, cex = 0.85)



## plot for Width
alldat <- ef[c(names(pdata0), "NEAR_X", "NEAR_Y")]
keep <- c("cWater_W")

for(i in setdiff(names(pdata0), keep)) {
  alldat[[i]] <- pdata0[[i]]
}
alldat $ fit <- predict(g1, newdata = alldat, type = "response")

alldat $ colgrp <- as.numeric(cut(alldat $ fit, breaks = breaks))
col <- cols[alldat $ colgrp]

plot(coast, border = grey(0.7))
plot(redctm, border = grey(0.8), add = TRUE)
with(alldat, {
  points(NEAR_X, NEAR_Y, col = col, pch = 16, cex = cex)
  points(NEAR_X, NEAR_Y, col = col, pch = 1, cex = cex)
})
text(190000, 570000, "Channel width", font = 2, adj = 1, cex = 0.9)



# plot lettering
mtext(c("a", "b"), side = 3, outer = TRUE, font = 2, cex = 1,
      line = -2, at = c(0, 1)/2.5 + 0.05)
mtext(c("c", "d"), side = 3, outer = TRUE, font = 2, cex = 1,
      line = -24.5, at = c(0, 1)/2.5 + 0.05)
mtext(c("e", "f"), side = 3, outer = TRUE, font = 2, cex = 1,
      line = -47.5, at = c(0, 1)/2.5 + 0.05)


dev.off()

}

