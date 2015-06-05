



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



#############################################################
#   
#   map of usable data
#
#############################################################

library(spdep)
coast <- rgdal::readOGR("mapdata", "britisles")[2,]
trustpoly <- rgdal::readOGR("mapdata","sa_trust_1110")

# load data
load("rData/modelData.rData")
# keep salmon, 3 pass fishings within the covariate bounds
ef3 <- subset(ef, Runs == 3 & Species == "Salmon" & keep)

{
png(file = "figures/data_map.png", width = 7, height = 9, units = "in", res = 400)

{
cols <- hexbin::BTC(nlevels(ef3 $ Trust))
names(cols) <- levels(ef3 $ Trust)
cols[c("MSS", "SEPA")] <- sapply(c("red", "gold3"), function(x) colorRampPalette(x)(1))
cols[!names(cols) %in% c("MSS", "SEPA")] <- hexbin::BTC(nlevels(ef3 $ Trust)-2, end = 150)

cols1 <- paste0(cols, "77")
cols2 <- paste0(cols, "DD")

par(mar = c(0,0,0,0))
  plot(trustpoly, border = grey(0.5), ylim = c(550000, 970000))
  #plot(coast, border = grey(0.5), add = TRUE)
  
  for (i in 1:nlevels(ef3 $ Trust)) {
    pdat <- subset(ef3, Trust == levels(Trust)[i]) 
    with(unique(pdat[c("NEAR_Y", "NEAR_X")]),
    {
      points(NEAR_X, NEAR_Y, pch = 16, col = cols1[i], cex = .7)    
    })
  }

  for (i in 1:nlevels(ef3 $ Trust)) {
    pdat <- subset(ef3, Trust == levels(Trust)[i]) 
    with(unique(pdat[c("NEAR_Y", "NEAR_X")]),
    {
      points(NEAR_X, NEAR_Y, pch = 1, col = cols2[i], cex = .7)    
    })
  }
 
}

dev.off()
}






#############################################################
#   
#   table of visits
#
#############################################################


load("rData/densmodelData.rData")
library(latticeExtra)
library(grid)
library(hexbin)

{
matrix.plot <- function (x, colorkey = FALSE, ...) {

  di <- dim(x)
  
  xlim = 0.5 + c(0, di[2])
  ylim = 0.5 + c(di[1], 0)
  
  df <- data.frame(y = c(row(x)), x = c(col(x)), z = c(x))

  cols <- function(n) BTC(n, beg = 200, end = 50)
                        
  
  #cols <- function(n) colorRampPalette(c("red","blue"))(n)
  #cols <- function(n) brewer.pal(n, "Spectral")

  zz <- df $ z
  rz <- range(zz, finite = TRUE)
  nn <- 100
  n0 <- min(nn, max(0, round((0 - rz[1])/(rz[2] - rz[1]) * nn)))
  levelplot(
    df$z ~ df$x + df$y,
    sub = "", ylab = "Organisation", xlab = "Year", xlim = xlim, ylim = ylim, 
    aspect = "iso", colorkey = colorkey, 
    col.regions = cols(11), cuts = 10, 
    panel = function(x, y, z, subscripts, at, ..., col.regions) 
      {
        x <- as.numeric(x)
        y <- as.numeric(y)
        numcol <- length(at) - 1
        num.r <- length(col.regions)
        col.regions <- if (num.r <= numcol) {
          rep(col.regions, length = numcol)
        } else {
          col.regions[1 + ((1:numcol - 1) * (num.r - 1))%/%(numcol - 1)]
        }
        zcol <- rep.int(NA_integer_, length(z))
        for (i in seq_along(col.regions)) {
          zcol[!is.na(x) &  !is.na(y) & !is.na(z) & at[i] <= z & z < at[i + 1]] <- i
        }
        wh <- grid::current.viewport()[c("width", "height")]
        wh <- c(grid::convertWidth(wh$width, "inches", valueOnly = TRUE), grid::convertHeight(wh$height, "inches", valueOnly = TRUE)) * par("cra")/par("cin")
        pSize <- wh/di
        pA <- prod(pSize)
        p1 <- min(pSize)
        lwd <- if (p1 < 2 || pA < 6)  0.01 else if (p1 >= 4) 1 else if (p1 > 3)  0.5 else 0.2
        grid.rect(x = x, y = y, width = 1, height = 1, default.units = "native", 
                  gp = gpar(fill = col.regions[zcol],  lwd = lwd, col = grey(0.7)))
        # plot value inside square, scale size depending on maximum dimension
        z <- paste(z)
        z[z=="NA"] <- ""
        grid.text(z, x = x, y = y, default.units = "native", 
                            gp = gpar(col = "black", cex = 0.8 * (10 / max(di))^.4, font = 2))
     },
     scales = list(x = list(at = 1:di[2], label = colnames(x), rot = 45, tck = c(1,0), alternating = 3),
                   y = list(at = 1:di[1], label = rownames(x), tck = c(1,0)), alternating = 3),
     par.settings=list(axis.line = list(col = grey(0.5)),
                       background = list(col = "transparent")),
     axis = function(side, scales, components, ...) {
      if (side %in% c("left", "bottom")) {
        axis.default(side, scales, components = components, ...)
      } else {

        axis.units <- lattice.getOption("axis.units")[["outer"]][[side]]
        axis.settings <- trellis.par.get("axis.components")[[side]]
        lab.unit <- unit(x = axis.settings$pad1 * axis.units$pad1$x, units = axis.units$pad1$units)

        if (side == "top") { 
          zsum <- paste(colSums(x, na.rm = TRUE))
          grid.text(zsum, x = unit(1:di[2], "native"), y = unit(1, "npc") + 2 * lab.unit, 
                    gp = gpar(col = "black", cex = 0.8 * (10 / max(di))^.4, font = 2)) 

        }
        if (side == "right") { 
          zsum <- paste(rowSums(x, na.rm = TRUE))
          grid.text(zsum, y = unit(1:di[1], "native"), x = unit(1, "npc") + 2 * lab.unit, 
                    gp = gpar(col = "black", cex = 0.8 * (10 / max(di))^.4, font = 2)) 

        }


    }
    }
  )
}


tab <- with(ef3, table(year, Trust))
tab[tab == 0] <- NA

tab <- tab[,order(colMeans(tab, na.rm = TRUE), decreasing = TRUE)]

matrix.plot(t(tab))
}

png(file = "figures/data_table_withkey.png", width = 7, height = 9, units = "in", res = 500)
p1 <- matrix.plot(t(tab), colorkey = TRUE)
print(p1)
dev.off()

png(file = "figures/data_table_nokey.png", width = 7, height = 9, units = "in", res = 500)
p1 <- matrix.plot(t(tab), colorkey = FALSE)
print(p1)
dev.off()



