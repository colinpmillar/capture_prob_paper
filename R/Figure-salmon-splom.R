


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



#############################################################
#   
#   splom of covariates
#
#############################################################


library(lattice)
library(gplots)

clsplom <- function(data) {
splom(data,
      colramp = BTC,
      xlab = "",
      diag.panel = function(x, ...) {
        yrng <- current.panel.limits()$ylim
        d <- density(x, na.rm=TRUE)
        d$y <- with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y) )
        panel.lines(d)
        diag.panel.splom(x, ...)
      },
      panel = function(x, y, colramp, zcuts = 25, ...) {
        nb <- as.integer(zcuts + 1)

        dx <- diff(rx <- range(x, na.rm = TRUE))
        xbreaks <- seq.int(rx[1L], rx[2L], length.out = nb)
        xbreaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + dx/1000)
 
        dy <- diff(ry <- range(y, na.rm = TRUE))
        ybreaks <- seq.int(ry[1L], ry[2L], length.out = nb)
        ybreaks[c(1L, nb)] <- c(ry[1L] - dy/1000, ry[2L] + dy/1000)

        if (all(round(x) == x) & diff(range(x)) <= zcuts) {
          xc <- factor(x, levels = min(x):max(x))
          xbreaks <- (min(x)-1):max(x) + 0.5
        } else {
          xc <- cut(x, xbreaks)
        }

        if (all(round(y) == y) & diff(range(y)) <= zcuts) {
          yc <- factor(y, levels = min(y):max(y))
          ybreaks <- (min(y)-1):max(y) + 0.5
        } else {
          yc <- cut(y, ybreaks)
        }

        dens <- as.matrix(table(xc, yc))
        dens[] <- log(dens / sum(dens))
        dens <- dens - max(dens)
        dens[!is.finite(dens)] <- min(dens[is.finite(dens)]) - log(2)

        gxy <- expand.grid(x = xbreaks[-1] - diff(xbreaks), y = ybreaks[-1] - diff(ybreaks))
        panel.levelplot(gxy $ x, gxy $ y, c(dens), col.regions = colramp, ...)
      },
      lower.panel = function(...) NULL,
      pscale=0, varname.cex=0.7, varname.font=2,
    )
}


clsplom2 <- function(data) {
splom(data,
      colramp = BTC,
      xlab = "",
      diag.panel = function(x, ...) {
        yrng <- current.panel.limits()$ylim
        d <- density(x, na.rm=TRUE)
        d$y <- with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y) )
        panel.lines(d)
        diag.panel.splom(x, ...)
      },
      lower.panel = function(x, y, colramp, zcuts = 25, ...) {
        nb <- as.integer(zcuts + 1)

        dx <- diff(rx <- range(x, na.rm = TRUE))
        xbreaks <- seq.int(rx[1L], rx[2L], length.out = nb)
        xbreaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + dx/1000)
 
        dy <- diff(ry <- range(y, na.rm = TRUE))
        ybreaks <- seq.int(ry[1L], ry[2L], length.out = nb)
        ybreaks[c(1L, nb)] <- c(ry[1L] - dy/1000, ry[2L] + dy/1000)

        if (all(round(x) == x) & diff(range(x)) <= zcuts) {
          xc <- factor(x, levels = min(x):max(x))
          xbreaks <- (min(x)-1):max(x) + 0.5
        } else {
          xc <- cut(x, xbreaks)
        }

        if (all(round(y) == y) & diff(range(y)) <= zcuts) {
          yc <- factor(y, levels = min(y):max(y))
          ybreaks <- (min(y)-1):max(y) + 0.5
        } else {
          yc <- cut(y, ybreaks)
        }

        dens <- as.matrix(table(xc, yc))
        dens[] <- log(dens / sum(dens))
        dens <- dens - max(dens)
        dens[!is.finite(dens)] <- min(dens[is.finite(dens)]) - log(2)

        gxy <- expand.grid(x = xbreaks[-1] - diff(xbreaks), y = ybreaks[-1] - diff(ybreaks))
        panel.levelplot(gxy $ x, gxy $ y, c(dens), col.regions = colramp, ...)
      },
      panel = function(...) NULL,
      pscale=0, varname.cex=0.7, varname.font=2,
    )
}


load("rData/densmodelData.rData")

which <- c("year", "doy", "Water_W", "Elevation_", "Distance_s", "sinSlope", 
           "Upcatch_km", "Urban" ,"woodland",
            "Marsh" ,"Other", "NEAR_X", "NEAR_Y", "totalN", "HACode")

fullnames <- data.frame(names = c("Year", 
                          "DoY", 
                          "Width", 
                          "Altitude", 
                          "DS", 
                          "Gradient", 
                          "UCA" ,
                          "Urban",
                          "Woodland",
                          "Marsh",
                          "Other",
                          "Lat", "Lon", "SalmonPass1", "HA"),
                        stringsAsFactors = FALSE)
rownames(fullnames) <- which

pdata <- ef3[which]
names(pdata) <- fullnames[names(pdata),]




p1 <- clsplom(pdata[c("Lat", "Lon", "Year", "DoY", "Width", "Altitude", "DS", "Gradient", "UCA")])
p2 <- clsplom2(pdata[c("Lat", "Lon", "Year", "SalmonPass1", "HA", "Urban","Woodland","Marsh","Other")])

p1
p2


png(file = "figures/data_splom.png", width = 9, height = 9, units = "in", res = 500)
{
dx <- 0.0955
print(p1, position = c(0,   0, 1-dx, 1), more = TRUE)
print(p2, position = c(dx, 0,   1, 1), more = FALSE)
}
dev.off()




 