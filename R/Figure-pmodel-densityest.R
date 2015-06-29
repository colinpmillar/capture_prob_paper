

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


coast <- rgdal::readOGR("mapdata", "britisles")[2,]


g1 <- best
var.summary <- g1 $ var.summary

library(mgcv)
library(CLmodel)
library(CLdata)

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
getPlotData <- function(var, func = function(x) exp(x) / (1+exp(x)), model = g1) {
  args <- pdata0[!names(pdata0) %in% var]
  args[var] <- pdata1[var]

  pdata <- do.call(expand.grid, args)

  # and predict
  pdata[c("fit", "se")] <- predict(model, newdata = pdata, se.fit = TRUE)
  pdata $ p <- func(pdata $ fit)
  pdata $ cil <- func(pdata $ fit - 2*pdata $ se * sqrt(2.13))
  pdata $ ciu <- func(pdata $ fit + 2*pdata $ se * sqrt(2.13))
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
redctm <- redctm[redctm $ CATCH_ID %in% ef3 $ CATCH_ID,]
hma <- hma[hma $ HACode %in% ef3 $ HACode,]

hma @ bbox <- bbox(redctm)
hma @ bbox[4] <- 1100000



alldat <- ef3[c(names(pdata0), "NEAR_X", "NEAR_Y", "T", "p", "psat", "constantp", "Area", "Runs")]

## values for modelled P
alldat $ fit1 <- with(alldat, T / (Area * (1 - (1-p)^Runs)))


## values for constant p
alldat $ fit2 <- with(alldat, T / (Area * (1 - (1-constantp)^Runs)))

## values for saturated p
alldat $ fit3 <- with(alldat, T / (Area * (1 - (1-psat)^Runs)))

alldat $ fitcp <- alldat $ fit2 / alldat $ fit1
alldat $ fitsp <- alldat $ fit3 / alldat $ fit1



#
##
##   2nd attempt using hex bins
##
#


library(hexbin)
# make a grid of predictions


# set up bins
dat <- hexbin(ef3 $ NEAR_X, ef3 $ NEAR_Y, xbins = 40, IDs = TRUE, shape = 1.5)
xbins <- dat@xbins
shape <- dat@shape
tmp <- hcell2xy(dat, check.erosion = TRUE)
cnt <- dat@count
sx <- xbins/diff(dat@xbnds)
sy <- (xbins * shape)/diff(dat@ybnds)
inner <- 0.5
outer <- (2 * inner)/sqrt(3)
dx <- inner/sx
dy <- outer/(2 * sy)
rad <- sqrt(dx^2 + dy^2)
hexC <- hexcoords(dx, dy, sep = NULL)




{
png("figures/pmodel-densityMaps.png", width = 9, height = 6, res = 500, units = "in")
par(mfrow = c(1,2), mar = c(0,0,0,0))

coast @ bbox <- bbox(redctm)
coast @ bbox[3] <- 440000
coast @ bbox[4] <- 1100000

# set up breaks
breaks <- seq(0.83, 1.17, by = 0.02)
nbreaks <- length(breaks) - 1

# set up colours
cols <- colorRampPalette(c(rich.colors(5)[5], rich.colors(5)[4], "white", rich.colors(5)[2], rich.colors(5)[1]))(nbreaks)
cols[length(cols)/2+0.5] <- rich.colors(15)[8]

# assign colours to cells
means <- tapply(alldat $ fitcp, dat@cID, mean)
range(means)
colgrpcp <- as.numeric(cut(means, breaks = breaks))
colcp <- cols[colgrpcp]
table(factor(colgrpcp, levels = 1:nbreaks))


# do plot 1
plot(coast, border = grey(0.7))
for (i in 1:dat@ncells) {
  polygon(hexC$x + tmp $ x[i], hexC$y + tmp $ y[i], col = colcp[i], border = grey(0.8))
}

if (FALSE) {
  # now for legend
  x0 <- 421799; y0 <- 947168; dx <- 22000; dy <- 300000/nbreaks
  val <- breaks
  val <- sprintf("%.2f", rev(val))
  for (i in 1:nbreaks-1) {
    polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - i*dy + c(0, -dy)[c(1,1,2,2)], col = rev(cols)[i+1]) 
  }
  text(x0 + 1.1*dx, y0 - 0:nbreaks*dy, val, font = 1, cex = 0.8, adj = 0)
} else {
  # now for legend
  x0 <- 440000; y0 <- 947168; dx <- 22000; dy <- 300000/nbreaks
  val <- breaks
  val <- sprintf("%.2f", rev(val))
  for (i in 1:nbreaks-1) {
    polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - i*dy + c(0, -dy)[c(1,1,2,2)], col = rev(cols)[i+1]) 
  }
  text(x0 + 1.1*dx, y0 - 0:nbreaks*dy, val, font = 1, cex = 0.8, adj = 0)  
  tabval <- rev(table(factor(colgrpcp, levels = 1:nbreaks)))
  tabval[tabval == 0] <- ""
  text(x0 - 0.1*dx, y0 - 1:nbreaks*dy + dy/2, tabval, font = 1, cex = 0.8, adj = 1)  
}





# set up breaks
breaks1 <- seq(0.83, 1.17, by = 0.02)
nbreaks1 <- length(breaks1) - 1
# tag on high values
breaks2 <- c(1.17, 1.25, 1.5, 2, 3, 4)
nbreaks2 <- length(breaks2) - 1

breaks <- c(breaks1, breaks2[-1])
nbreaks <- length(breaks) - 1

# set up colours
cols1 <- colorRampPalette(c(rich.colors(5)[5], rich.colors(5)[4], "white", rich.colors(5)[2], rich.colors(5)[1]))(nbreaks1)
cols1[length(cols1)/2+0.5] <- rich.colors(15)[8]

cols2 <- colorRampPalette(c("purple", "hotpink"))(nbreaks2)

cols <- c(cols1, cols2)


# assign colours to cells
alldat $ fitsp_noInf <- alldat $ fitsp
alldat $ fitsp_noInf[alldat $ fitsp_noInf > 1e6] <- NA
means <- tapply(alldat $ fitsp_noInf, dat@cID, mean, na.rm = TRUE)
means[is.nan(means)] <- NA
range(means, na.rm = TRUE)
colgrpsp <- as.numeric(cut(means, breaks = breaks))
colsp <- cols[colgrpsp]
table(colgrpsp)


# do plot 2
plot(coast, border = grey(0.7))
for (i in 1:dat@ncells) {
  polygon(hexC$x + tmp $ x[i], hexC$y + tmp $ y[i], col = colsp[i], border = grey(0.8))
}



if (FALSE) {
  # now for legend
  #x0 <- 421799; y0 <- 947168; dx <- 22000; dy <- 300000/nbreaks
  y0 <- y0 + dy*nbreaks2
  val <- breaks
  val <- sprintf("%.2f", rev(val))
  for (i in 1:nbreaks-1) {
    polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - i*dy + c(0, -dy)[c(1,1,2,2)], col = rev(cols)[i+1]) 
  }
  text(x0 + 1.1*dx, y0 - 0:nbreaks*dy, val, font = 1, cex = 0.8, adj = 0)
} else {
  # now for legend
  #x0 <- 440000; y0 <- 947168; dx <- 22000; dy <- 300000/nbreaks
  y0 <- y0 + dy*nbreaks2
  val <- breaks
  val <- sprintf("%.2f", rev(val))
  for (i in 1:nbreaks-1) {
    polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - i*dy + c(0, -dy)[c(1,1,2,2)], col = rev(cols)[i+1]) 
  }
  text(x0 + 1.1*dx, y0 - 0:nbreaks*dy, val, font = 1, cex = 0.8, adj = 0)  
  tabval <- rev(table(factor(colgrpsp, levels = 1:nbreaks)))
  tabval[tabval == 0] <- ""
  text(x0 - 0.1*dx, y0 - 1:nbreaks*dy + dy/2, tabval, font = 1, cex = 0.8, adj = 1)  
}


mtext(c("a","b"), side = 3, line = -4, outer = TRUE, font = 2, at = c(0, 0.5) + 0.05)
dev.off()
}






## a plot of data coverage
# set up breaks
breaks <- c(0,1,2,5,10,20,40,60,90)
nbreaks <- length(breaks) - 1

# set up colours
cols <- colorRampPalette(c("purple", "cyan"))(nbreaks)
cols <- topo.colors(nbreaks)

# assign colours to cells
means <- table(dat@cID)
range(means)
colgrpcp <- as.numeric(cut(means, breaks = breaks))
colcp <- cols[colgrpcp]
table(colgrpcp)


# do plot 1
coast @ bbox <- bbox(redctm)
coast @ bbox[4] <- 1100000
plot(coast, border = grey(0.7))
for (i in 1:dat@ncells) {
  polygon(hexC$x + tmp $ x[i], hexC$y + tmp $ y[i], col = colcp[i], border = grey(0.8))
}

# now for legend
x0 <- 421799; y0 <- 947168; dx <- 22000; dy <- 300000/nbreaks
val <- breaks
val <- sprintf("%.2f", rev(val))
for (i in 1:nbreaks-1) {
polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - i*dy + c(0, -dy)[c(1,1,2,2)], col = rev(cols)[i+1])
}
text(x0 + 1.1*dx, y0 - 0:nbreaks*dy, val, font = 1, cex = 1, adj = 0)








#
##
##  Old attempt
##
#


if (FALSE) {

png("figures/pmodel-densityMaps.png", width = 9, height = 6, res = 500, units = "in")

nbreaks <- 14
x <- c(alldat $ fitcp) 
nb <- as.integer(nbreaks + 1)
rx <- range(x)
rx <- 1 + c(-1,1)*max(abs(rx - 1))
rx <- c(0.65, 1.35)
dx <- diff(rx)
breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + dx/1000)

cols <- rev(rich.colors(nbreaks))

alldat $ colgrpcp <- as.numeric(cut(alldat $ fitcp, breaks = breaks))
alldat $ colcp <- cols[alldat $ colgrpcp]
table(alldat $ colgrpcp)

par(mfrow = c(1,2), mar = c(0,0,0,0))

plot(hma, border = grey(0.7))
plot(redctm, border = grey(0.8), add = TRUE)
ord <- order(abs(alldat $ fitcp - 1), decreasing = FALSE)
cex <- 0.7
with(alldat[ord,], {
  points(NEAR_X, NEAR_Y, col = paste0(substring(colcp,1,7), "55"), pch = 16, cex = cex)
})
with(alldat[ord,], {
  points(NEAR_X, NEAR_Y, col = paste0(substring(colcp,1,7), "AA"), pch = 1, cex = cex)
})


# now for legend
x0 <- 421799; y0 <- 947168; dx <- 22000; dy <- 300000/nbreaks
val <- breaks
val <- sprintf("%.2f", rev(val))
startb <- 6
for (i in startb:nbreaks-1) {
polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - i*dy + c(0, -dy)[c(1,1,2,2)], col = rev(cols)[i+1])
}
val[1:(startb-1)] <- ""
text(x0 + 1.1*dx, y0 - 0:nbreaks*dy, val, font = 2, cex = 1, adj = 0)



nbreaks <- 20
nb <- as.integer(nbreaks + 1)
dx <- diff(rx <- c(0, 2))
breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + dx/1000)
breaks[1] <- 0
breaks[nbreaks+1] <- Inf

cols <- rev(rich.colors(nbreaks))

alldat $ colgrpsp <- as.numeric(cut(alldat $ fitsp, breaks = breaks))
colsp <- cols[alldat $ colgrpsp]

plot(hma, border = grey(0.7))
plot(redctm, border = grey(0.8), add = TRUE)

ord <- order(abs(alldat $ fitsp - 1), decreasing = FALSE)
cex <- 0.7
with(alldat[ord,], {
  points(NEAR_X, NEAR_Y, col = paste0(substring(colsp,1,7), "55"), pch = 16, cex = cex)
})
with(alldat[ord,], {
  points(NEAR_X, NEAR_Y, col = paste0(substring(colsp,1,7), "AA"), pch = 1, cex = cex)
})

# now for legend
x0 <- 421799; y0 <- 947168; dx <- 22000; dy <- 300000/nbreaks
val <- breaks
val <- sprintf("%.2f", rev(val))
val[1] <- paste0(val[2], "+")
for (i in 1:nbreaks-1) {
polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - i*dy + c(0, -dy)[c(1,1,2,2)], col = rev(cols)[i+1])
}
text(x0 + 1.1*dx, y0 - 0:nbreaks*dy, val, font = 2, cex = 1, adj = 0)


dev.off()

}
