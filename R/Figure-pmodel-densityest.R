
if (Sys.info()["user"] == "millaco") {
  setwd("~/work/SMFS-report")    
} else 
if (Sys.info()["user"] == "millarc") {
  setwd("B:/Conservation_Limits/CL_Juvenile_Density/SMFS-report")
} else 
if (Sys.info()["user"] == "Millarc") {
  setwd("C:/work/SMFS-report")
}

# load fits and model data
load("rData/bestdensmodel.rData")
load("rData/densmodeldata.rData")

coast <- rgdal::readOGR("../GIS_shapefiles_for_model/Coastline", "britisles")[2,]


g1 <- best
var.summary <- g1 $ var.summary

library(mgcv)
library(CLmodel)
library(CLdata)

g1 <- best
var.summary <- g1 $ var.summary
var.type <- sapply(var.summary, function(x) is(x)[1])
var.names <- names(var.summary)


which <- c("year", "doy", "Water_W", "Elevation_", "Distance_s", "sinSlope", 
           "Upcatch_km" ,"CTrees", "Urban" ,"NCTrees" ,"Mixed" ,
            "Marsh" ,"Other","Trust", "fyear", "HACode", "CATCH_ID")
fullnames <- data.frame(names = c("Year", 
                          "DoY", 
                          "Width", 
                          "Altitude", 
                          "DS", 
                          "Gradient", 
                          "UCA" ,
                          "Conifer", 
                          "Urban",
                          "Deciduous",
                          "Mixed",
                          "Marsh",
                          "Other",
                          "Organisation",
                          "Year (factor)",
                          "HA", "Catchment"),
                        stringsAsFactors = FALSE)
rownames(fullnames) <- which

var.fullnames <- fullnames[var.names,]


# set base prediction levels
pdata0 <- ifelse(var.type == "numeric", 
                   lapply(var.summary, "[", 2), 
                   lapply(var.summary, function(x) levels(x)[floor(nlevels(x)/2)])
                )
pdata0 $ HACode <- 15
pdata0 $ CATCH_ID <- 1
pdata0 $ fyear <- "2006"
pdata0 $ Trust <- "MSS"

# set up prediction ranges
pdata1 <- ifelse(var.type == "numeric", 
                   lapply(var.summary, function(x) if (length(x) == 3) seq(x[1], x[3], length=100) else 0), 
                   lapply(var.summary, function(x) levels(x))
                )
pdata1 $ fyear <- sort(unique(g1 $ model $ fyear))
pdata1 $ HACode <- unique(g1 $ model $ HACode)
pdata1 $ CATCH_ID <- unique(g1 $ model $ CATCH_ID)

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

data(redctm)
data(hma)
redctm <- redctm[redctm $ CATCH_ID %in% ef $ CATCH_ID,]
hma <- hma[hma $ HACode %in% ef $ HACode,]

hma @ bbox <- bbox(redctm)
hma @ bbox[4] <- 1100000

{
png("figures/pmodel-densityMaps.png", width = 9, height = 9, res = 500, units = "in")



par(mfrow = c(1,1), mar = c(0,0,0,0))


## plot for modelled P
alldat <- ef[c(names(pdata0), "NEAR_X", "NEAR_Y", "T", "p", "Area", "Runs")]
alldat $ fit1 <- with(alldat, T / (Area * (1 - (1-p)^Runs)))


## plot for constant p
alldat $ p <- 0.5277
alldat $ fit2 <- with(alldat, T / (Area * (1 - (1-p)^Runs)))

alldat $ fit <- alldat $ fit2 / alldat $ fit1

nbreaks <- 10
x <- alldat $ fit
nb <- as.integer(nbreaks + 1)
dx <- diff(rx <- range(x, na.rm = TRUE))
breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + dx/1000)

cols <- rev(rich.colors(nbreaks))


alldat $ colgrp <- as.numeric(cut(alldat $ fit, breaks = breaks))
col <- cols[alldat $ colgrp]

plot(hma, border = grey(0.7))
plot(redctm, border = grey(0.8), add = TRUE)
with(alldat, {
  points(NEAR_X, NEAR_Y, col = col, pch = 16)
  points(NEAR_X, NEAR_Y, col = col, pch = 1)
})


# now for legend
x0 <- 421799; y0 <- 947168; dx <- 42000; dy <- 300000/nbreaks
val <- breaks
val <- sprintf("%.2f", rev(val))
for (i in 1:nbreaks-1) {
polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - i*dy + c(0, -dy)[c(1,1,2,2)], col = rev(cols)[i+1])
}
text(x0 + 1.1*dx, y0 - 0:nbreaks*dy, val, font = 2, cex = 1, adj = 0)

dev.off()

}

