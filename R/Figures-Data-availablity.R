



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

#cols <- hexbin::BTC(5)
cols <- substring(rep(gplots::rich.colors(5), each = 5), 1, 7)
names(cols) <- levels(ef3 $ Trust)

cols1 <- paste0(cols, "77")
cols2 <- paste0(cols, "DD")

pch <- rep(c(15, 16, 17, 18, 25), 5)
cex <- 0.8

par(mar = c(0,0,0,0))
  plot(trustpoly, border = grey(0.5), ylim = c(550000, 970000))
  #plot(coast, border = grey(0.5), add = TRUE)
  
  for (i in 1:nlevels(ef3 $ Trust)) {
    pdat <- subset(ef3, Trust == levels(Trust)[i]) 
    pdat <- unique(pdat[c("NEAR_Y", "NEAR_X")])
    points(pdat $ NEAR_X, pdat $ NEAR_Y, pch = pch[i], col = cols1[i], cex = cex)
  }

#  for (i in 1:nlevels(ef3 $ Trust)) {
#    pdat <- subset(ef3, Trust == levels(Trust)[i]) 
#    pdat <- unique(pdat[c("NEAR_Y", "NEAR_X")])
#    points(pdat $ NEAR_X, pdat $ NEAR_Y, pch = pch[i], bg = "transparent", col = cols2[i], cex = cex)    
#  }
 

dev.off()
}




