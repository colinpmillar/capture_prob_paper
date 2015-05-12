


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

which <- c("year", "doy", "Water_W", "Elevation_", "Distance_s", "sinSlope", 
           "Upcatch_km" ,"CTrees", "Urban" ,"NCTrees" ,"Mixed" ,
            "Marsh" ,"Other","Trust", "fyear", "HACode")
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
                          "HA"),
                        stringsAsFactors = FALSE)
rownames(fullnames) <- which

var.fullnames <- fullnames[var.names,]

# set base prediction levels
pdata0 <- ifelse(var.type == "numeric", 
                   lapply(var.summary, "[", 2), 
                   lapply(var.summary, function(x) levels(x)[floor(nlevels(x)/2)])
                )
pdata0 $ Trust <- "MSS"
pdata0 $ fyear <- "2006"

# set up prediction ranges
pdata1 <- ifelse(var.type == "numeric", 
                   lapply(var.summary, function(x) if (length(x) == 3) seq(x[1], x[3], length=100) else 0), 
                   lapply(var.summary, function(x) levels(x))
                )
pdata1 $ fyear <- paste(1997:2013)
pdata1 $ HACode <- sort(unique(ef $ HACode))


# set up prediction data frame
getPlotData <- function(var, func = function(x) exp(x) / (1+exp(x)), model = g1) {
  args <- pdata0[!names(pdata0) %in% var]
  args[var] <- pdata1[var]

  pdata <- do.call(expand.grid, args)

  # and predict
  pdata[c("fit", "se")] <- predict(model, newdata = pdata, se.fit = TRUE)
  pdata $ p <- func(pdata $ fit)
  pdata $ cil <- func(pdata $ fit - 2*pdata $ se)
  pdata $ ciu <- func(pdata $ fit + 2*pdata $ se)
  pdata $ var <- paste(var, collapse = ":")

  if (length(var) > 1) {
    stop()
    #pdata $ x <- as.numeric(do.call(interaction, c(pdata[var], list(sep = " "))))
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

{
png(file = "figures/pmodel_grid.png", width = 7, height = 7, units = "in", res = 500)

ylim <- c(0.1, 0.8)

par(mar = c(5,2.5,3,1)) # c(bottom, left, top, right)

layout(rbind(c(1,1,2), c(3,4,4), c(5,6,7)))


#   Trust predictions of p
# ----------------------------------------
pdata <- getPlotData("Trust")
pdata <- pdata[order(pdata $ p),]
pdata $ Trust <- as.character(pdata $ Trust)
pdata $ Trust[pdata $ Trust == "Other"] <- "Caithness"
factorPlot(pdata, xlab = "Organisation", ylim = ylim, labcex = 0.8)

#   DoY predictions of p
# ----------------------------------------
pdata <- getPlotData("doy")
continuousPlot(pdata, xlab = fullnames["doy",], rug = ef $ doy, ylim = ylim, yaxislab = FALSE)


#   Year predictions of p
# ----------------------------------------
pdata <- getPlotData("fyear")
factorPlot(pdata, xlab = "Year", ylim = ylim)


#   HA predictions of p
# ----------------------------------------
pdata <- getPlotData("HACode")   
# order north to south?
ns <- coordinates(hma)[,2]
names(ns) <- hma $ HACode
ord <- order(ns[paste(pdata $ HACode)])
#ord <- order(pdata $ p)
pdata <- pdata[ord,]
factorPlot(pdata, xlab = fullnames["HACode",], ylim = ylim, labcex = 0.8, yaxislab = FALSE)


#   other predictions of p
# ----------------------------------------
covars <- c("Distance_s", "Water_W", "sinSlope")

for (i in seq_along(covars)) {
  pdata <- getPlotData(covars[i])
  continuousPlot(pdata, xlab = fullnames[covars[i],], rug = ef[[covars[i]]], ylim = ylim, i %in% c(1,4))  
} 

dev.off()

}
}
