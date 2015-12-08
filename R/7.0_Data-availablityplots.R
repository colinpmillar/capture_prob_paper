



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
coast @ bbox <- bbox(trustpoly)


# load data
load("intermediate_rData/screenedData.rData") # ef

load("rData/modelData.rData")
# strip off unused columns
ef <- ef[!names(ef) %in% c(paste0("n_R", 4:6), "pDate")]
# keep salmon, 3 pass fishings within the covariate bounds
ef <- subset(ef, Runs > 2 & Species == "Salmon" & keep & Trust != "Nith")

ef3 <- ef
ef3 $ Trust <- factor(ef3 $ Trust)
levs <- c("MSS", "SEPA", levels(ef3$Trust)[!levels(ef3$Trust) %in% c("MSS", "SEPA")])
ef3 $ Trust <- factor(ef3 $ Trust, levels = levs)


{
#png(file = "figures/data_map.png", width = 13, height = 7, units = "in", res = 400)
png(file = "C:/work/Dropbox/CaptureProbPaper/resubmission/Figure1.png", width = 13, height = 7, units = "in", res = 400)


#cols <- hexbin::BTC(5)
#cols <- rep(gplots::rich.colors(6), each = 4)
cols <- rep(RColorBrewer::brewer.pal(6, "Dark2"), each = 4)
#cols <- rep(c("indianred1", "orange", "chartreuse", "deepskyblue", "darkmagenta", "grey07"), each = 4)
pch <- rep(c(1, 2, 5, 6), 6)

colsym <- data.frame(Trust = levels(ef3 $ Trust),
                     pch = pch,
                     col = cols,
                     cex = 0.7,
                     stringsAsFactors = FALSE)
rownames(colsym) <- colsym $ Trust

colsym[c("MSS","SEPA"),"pch"] <- c(15,16)
colsym[c("MSS","SEPA"),"cex"] <- c(1, 1)
colsym[c("MSS","SEPA"),"col"] <- 1

par(mar = c(0,0,0,0), mfrow = c(1,2))


library(CLdata)
data(hma)
data(redctm)

plot(redctm, border = grey(0.4), ylim = c(550000, 970000))
plot(hma, border = 1, add = TRUE)


  plot(redctm, border = "transparent", ylim = c(550000, 970000))
  plot(trustpoly, border = grey(0.5), add = TRUE)
  #plot(coast, border = grey(0.5), add = TRUE)

  for (i in 1:nrow(colsym)) {
    pdat <- subset(ef3, Trust == colsym$Trust[i])
    pdat <- unique(pdat[c("NEAR_Y", "NEAR_X")])
    points(pdat $ NEAR_X, pdat $ NEAR_Y,
           pch = colsym $ pch[i], col = colsym $ col[i],
           cex = colsym $ cex[i])
  }

# legend
  x0 <- -6000.456; y0 <- 840253.2; dx <- 22000; dy <- -300000/25
  y0 <- y0 + dy*25
  points(rep(x0, 24), y0 - 24:1*dy, cex = colsym$cex, pch = colsym$pch, col = colsym$col, xpd = NA)
  text(x0 + 0.4*dx, y0 - 24:1*dy, colsym $ Trust, font = 1, cex = 1.1, adj = 0, xpd = NA)

 mtext(c("a","b"), side = 3, line = -4, outer = TRUE, font = 2, at = c(0, 0.5) + 0.05, cex = 2)



dev.off()
}




