


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


# load fits and model data
load("rData/bestpmodel.rData")
load("rData/densmodelData.rData")



# #################################
#
# setup prediction data
#
# #################################

var.summary <- best $ Gsetup $ var.summary
var.type <- sapply(var.summary, function(x) is(x)[1])
var.names <- names(var.summary)

which <- c("LifeStage", "Trust", "fyear", "Elevation_", "Distance_s", "totalN", 
           "Water_W" ,"doy")
fullnames <- data.frame(names = c("Lifestage", 
                          "Organisation", 
                          "Year", 
                          "Altitude", 
                          "DS", 
                          "SalmonPass1", 
                          "Width" ,
                          "DoY"),
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
pdata1 $ fyear <- paste(1997:2013)
pdata1 $ totalN <- seq(0, 1000, length = 100)



# set up prediction data frame
getPlotData <- function(var, func = function(x) exp(x) / (1+exp(x)), model = g1, phi = 2.434) {
  args <- pdata0[!names(pdata0) %in% var]
  args[var] <- pdata1[var]

  pdata <- do.call(expand.grid, args)

  # and predict
  pdata[c("fit", "se")] <- predict(model, newdata = pdata, se.fit = TRUE)
  pdata $ p <- func(pdata $ fit)
  pdata $ cil <- func(pdata $ fit - 2*pdata $ se * sqrt(phi))
  pdata $ ciu <- func(pdata $ fit + 2*pdata $ se * sqrt(phi))
  pdata $ var <- paste(var, collapse = ":")

  if (length(var) > 1) {
    #stop()
    #pdata $ x <- as.numeric(do.call(interaction, c(pdata[var], list(sep = " "))))
    pdata $ x <- as.numeric(pdata[[var[1]]])
    pdata $ cx <- paste(pdata[[var[1]]]) 
  } else {
    pdata $ x <- as.numeric(pdata[[var]])
    pdata $ cx <- paste(pdata[[var]])    
  }
  
  pdata
}



# define the plotting functions
{

factorPlot <- function(pdata, xlab = "", ylim = NULL, yaxislab = TRUE, labcex = 0.9) {

  if (is.null(ylim)) {
    ylim <- range(pdata $ cil, pdata $ ciu)
  }

  plot(pdata $ p, pch = 16, ylim = ylim, axes = FALSE, ann = FALSE)
  segments(x0 = 1:nrow(pdata), y0 = ylim[1], y1 = pdata $ p, col = grey(0.9), lty = 2)
  segments(x0 = 1:nrow(pdata), y0 = pdata $ cil, y1 = pdata $ ciu)
  axis(2, las = 1, labels = yaxislab)
  axis(1, at = 1:nrow(pdata), labels = FALSE, tck = FALSE)
  text(1:nrow(pdata), ylim[1] - diff(ylim)*.08,  pdata[[pdata $ var[1]]], srt=45, xpd = TRUE, adj = 1, cex = labcex)
  title(ylab = "", main = xlab)
  box(bty = "l")
}


continuousPlot <- function(pdata, ylim = NULL, xlab = "", rug = NULL, yaxislab = TRUE) {

  if (is.null(ylim)) {
    ylim <- range(pdata $ cil, pdata $ ciu)
  } 

  plot(pdata $ x, pdata $ p, type = "n", ylim = ylim, axes = FALSE, ann = FALSE)
#  lines(pdata $ x, pdata $ cil, lty = 2)
#  lines(pdata $ x, pdata $ ciu, lty = 2)
  with(pdata, polygon(c(x, rev(x)), c(cil, rev(ciu)), border = NA, col = "lightblue"))
  lines(pdata $ x, pdata $ p)
  axis(2, las = 1, labels = yaxislab)
  axis(1)
  title(ylab = "", main = xlab)
  box(bty = "l")
  # a rug!
  if (!is.null(rug)) {
    dx <- diff(rx <- range(rug, na.rm = TRUE))
    breaks <- seq.int(rx[1], rx[2], length.out = 101)

    dens <- table(cut(rug, breaks))
    cols <- heat.colors(11)

    dz <- diff(rz <- range(dens[dens > 0], na.rm = TRUE))
    zbreaks <- seq.int(rz[1], rz[2], length.out = 11)

    zs <- as.numeric(cut(dens, zbreaks))
    for (i in 1:100) {
      polygon(breaks[i + 0:1][c(1,2,2,1)], ylim[1] + c(0, diff(ylim)*.05)[c(1,1,2,2)], 
              col = cols[zs[i]], border = NA)
    }
  }
}


continuousPlot2 <- function(pdata, ylim = NULL, xlab = "", rug = NULL, yaxislab = TRUE) {

  if (is.null(ylim)) {
    ylim <- range(pdata $ cil, pdata $ ciu)
  } 

  plot(pdata $ x, pdata $ p, type = "n", ylim = ylim, axes = FALSE, ann = FALSE)
#  lines(pdata $ x, pdata $ cil, lty = 2)
#  lines(pdata $ x, pdata $ ciu, lty = 2)

  by <- strsplit(pdata $ var, ":")[[1]][2]
  cols <- c("blue", "red")
  cols <- colorRampPalette(cols)(2)
  colsp <- paste0(cols, "33")
  for (i in 1:2) {
    with(subset(pdata, pdata[[by]] == unique(pdata[[by]])[i]), {
      polygon(c(x, rev(x)), c(cil, rev(ciu)), border = NA, col = colsp[i])
    })
  }
  for (i in 1:2) {
    with(subset(pdata, pdata[[by]] == unique(pdata[[by]])[i]), {
      lines(x, p, col = cols[i])
    })
  }

  axis(2, las = 1, labels = yaxislab)
  axis(1)
  title(ylab = "", main = xlab)
  box(bty = "l")
  # a rug!
  if (!is.null(rug)) {
    dx <- diff(rx <- range(rug, na.rm = TRUE))
    breaks <- seq.int(rx[1], rx[2], length.out = 101)

    dens <- table(cut(rug, breaks))
    cols <- heat.colors(11)

    dz <- diff(rz <- range(dens[dens > 0], na.rm = TRUE))
    zbreaks <- seq.int(rz[1], rz[2], length.out = 11)

    zs <- as.numeric(cut(dens, zbreaks))
    for (i in 1:100) {
      polygon(breaks[i + 0:1][c(1,2,2,1)], ylim[1] + c(0, diff(ylim)*.05)[c(1,1,2,2)], 
              col = cols[zs[i]], border = NA)
    }
  }
}


{
png(file = "figures/pmodel_grid.png", width = 7, height = 7, units = "in", res = 400)

ylim <- c(0.33, .8)

par(mar = c(5,2.5,3,1)) # c(bottom, left, top, right)

layout(rbind(c(1,1,2), c(3,4,4), c(5,6,7)))


#   Trust predictions of p
# ----------------------------------------
pdata <- getPlotData("Trust")
pdata <- pdata[order(pdata $ p),]
pdata $ Trust <- as.character(pdata $ Trust)
pdata $ Trust[pdata $ Trust == "Other"] <- "Caithness"
factorPlot(pdata, xlab = "Organisation", ylim = ylim, labcex = 0.8)


#   Lifestage predictions of p
# ----------------------------------------
pdata <- getPlotData("LifeStage")
factorPlot(pdata, xlab = fullnames["LifeStage",], ylim = ylim, labcex = 0.8, yaxislab = FALSE)


#   DoY predictions of p
# ----------------------------------------
pdata <- getPlotData(c("doy", "LifeStage"))
continuousPlot2(pdata, xlab = "Lifestage x DoY", rug = ef3 $ doy, ylim = ylim, yaxislab = FALSE)


#   Year predictions of p
# ----------------------------------------
pdata <- getPlotData("fyear")
factorPlot(pdata, xlab = "Year", ylim = ylim, yaxislab = FALSE)


#   other predictions of p
# ----------------------------------------
covars <- c("Distance_s", "Elevation_", "Water_W")

for (i in seq_along(covars)) {
  pdata <- getPlotData(covars[i])
  continuousPlot(pdata, xlab = fullnames[covars[i],], rug = ef3[[covars[i]]], ylim = ylim, i == 1)  
} 

dev.off()

}

}
