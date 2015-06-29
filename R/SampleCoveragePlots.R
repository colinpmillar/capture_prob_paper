# ------------------------------------------------
# 
#  Work out crossover between Trust and MSS / SEPA
# 
# ------------------------------------------------


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

# load data
library(CLdata)
library(sp)
load("rData/modelData.rData")

ef <- subset(ef, keep)
ef <- unique(ef[!names(ef) %in% c("Species","Lifestage", "n_R1", "n_R2", "n_R3", "n_R4", "n_R5", "n_R6", "s", "T", "X", "Z", "phi")])

trustpoly <- rgdal::readOGR("mapdata","sa_trust_1110")
trustnames <- read.csv("trust_names.csv")
gis @ proj4string <- trustpoly @ proj4string
# which points are in which trust area:
wk <- sp::over(gis, trustpoly)
ef <- cbind(ef, wk[ef $ Site_OBJECTID,])

# table trusts against the area they fish in
tab <- table(ef $ TrustName, ef $ Trust)

ord <- order(apply(tab, 2, function(x) max(x)/sum(x)))

write.csv(tab[,rev(ord)], file = "TrustFishings_by_TrustRegions.csv")




# first with hydrometric area
tab <- with(ef, table(HAName, Trust))
ord <- order(apply(tab, 1, function(x) max(x)/sum(x)))
write.csv(tab[rev(ord),], file = "TrustFishings_by_HydroArea.csv")


sngl <- rowSums(tab > 1)
range(sngl)
t(t(sngl[sngl < 2]))

#Annan                      1
#Beauly                     1
#Firth of Tay Group         1
#Helmsdale Group            1
#Leven (Durnbartonshire)    1
#Lochy (Invernesshire)      1
#Outer Hebrides             1
#Tweed                      1
#Ythan Group                1

hma @ data $ trusts <- sngl[hma @ data $ HAName]
cols <- c("transparent", "red", "lightblue")
colid <- replace(hma @ data $ trusts, hma @ data $ trusts > 1, 2) + 1
#cols <- c("transparent", "red", colorRampPalette(c("lightblue", "darkblue"))(6))
#colid <- hma @ data $ trusts + 1
sp::plot(hma, col = cols[colid])
xy <- coordinates(hma)
text(xy[,1], xy[,2], label = hma @ data $ trusts, cex = 0.5, font = 2)

png(file = "figures/HMA_single_trust_map.png", width = 9, height = 9, units = "in", res = 400)
plot(hma, col = cols[colid])
dev.off()

png(file = "figures/HMA_single_trust_map_wnos.png", width = 9, height = 9, units = "in", res = 400)
plot(hma, col = cols[colid])
text(xy[,1], xy[,2], label = hma @ data $ trusts, cex = 0.5, font = 2)
dev.off()


# look at sites with wide widths from single sampled HAs

HAs <- names(sngl[sngl < 2])[-(4:5)]
ef_sub <- subset(ef, HAName %in% HAs & Water_W > 20)
by(ef_sub, ef_sub $ HAName[drop = TRUE], function(x) unique(x[c("Site_OBJECTID","Site.Name", "Trust","NEAR_X","NEAR_Y")]))

by(ef_sub, ef_sub $ HAName[drop = TRUE], function(x) unique(x[c("Site_OBJECTID","Trust","NEAR_X","NEAR_Y")]))


# plot these up on a map
tweedwater <- rgdal::readOGR("B:/Env_GIS/Tweed_Catchment/Shapes_Tweed","OSMM_Tweed_InWater")


plot(tweedwater)
points(ef_sub $ NEAR_X, ef_sub $ NEAR_Y, col = "red", cex = 16)

bb <- bbox(cbind(ef_sub $ NEAR_X, ef_sub $ NEAR_Y))




# now with catchments
tab  <- with(ef, table(CATCH_ID, Trust))
sngl <- rowSums(tab > 1)

redctm @ data $ trusts <- unname(sngl[paste(redctm @ data $ CATCH_ID)])
redctm @ data $ trusts[is.na(redctm @ data $ trusts)] <- 0
cols <- c("transparent", "red", "lightblue")
colid <- ifelse(redctm @ data $ trusts == 0, 1, 
           ifelse(redctm @ data $ trusts == 1, 2, 3))
plot(hma, col = grey(0.9), border = grey(0.9))
plot(redctm, col = cols[colid], add = TRUE)

png(file = "figures/catchment_single_trust_map.png", width = 9, height = 9, units = "in", res = 400)
plot(hma, col = grey(0.8), border = grey(0.8))
plot(redctm, col = cols[colid], add = TRUE)
dev.off()


# now for some barcharts by HA
tab <- with(ef, table(HAName, Trust))

df <- data.frame(HAName = rep(rownames(tab), ncol(tab)),
                 Trust = rep(colnames(tab), each = nrow(tab)),
                 val = c(tab))
df <- df[df $ val > 0,]

library(ggplot2)
library(dplyr)


ggplot(df, aes(HAName, y = val, fill = Trust)) +
   geom_bar(stat = "identity") +
   coord_flip() +
   xlab("") +
   ylab("Number of samples") 


library(dplyr)

df3 <- df %>% 
  group_by(HAName) %>% 
  mutate(perc=val/sum(val)) %>%
  mutate(max=max(perc))

levls <- paste(unique(df3 $ HAName[order(df3 $ max)]))
df3 $ HANamesort <- factor(df3 $ HAName, levels = levls)

png(file = "figures/coverage1.png", width = 9, height = 6, units = "in", res = 400)
ggplot(df3, aes(HANamesort, y = perc*100, fill = Trust)) +
   geom_bar(stat = "identity", color = "black") +
   coord_flip() +
   xlab("") +
   ylab("Percentage of samples") +
   scale_colour_grey()
dev.off()

png(file = "figures/coverage2.png", width = 9, height = 6, units = "in", res = 400)
ggplot(df3, aes(HANamesort, y = val, fill = Trust)) +
   geom_bar(stat = "identity", color = "black") +
   coord_flip() +
   xlab("") +
   ylab("Number of samples") 
dev.off()

hma @ data $ lev <- 3

tab <- subset(df3, max == 1)
tab $ HAName

hma @ data $ lev[hma @ data $ HAName %in% tab $ HAName] <- 1

tab <- subset(df3, max > 0.95 & val < 10)
tab $ HAName

hma @ data $ lev[hma @ data $ HAName %in% tab $ HAName] <- 2
hma @ data $ lev[hma @ data $ trusts == 0] <- 0
cols <- c("transparent", "red", "orange", "lightblue")

png(file = "figures/HMA_complex_trust_map.png", width = 9, height = 9, units = "in", res = 400)
plot(hma, col = cols[hma @ data $ lev + 1])
dev.off()


targets <- as.character(hma @ data $ HAName[hma @ data $ lev %in% 1:2])

keep <- c("Site_OBJECTID", "Dataset", "Width", "Trust", "NEAR_X", "NEAR_Y", 
           "Elevation_", "Slope_deg", "Upcatch_km", "Water_A", "Water_W", "Distance_s", 
           "Urban", "CTrees", "NCTrees", "Mixed", "Marsh", "Other",
           "DESCRIPTIO", "barrier", "HACode", "HAName", "doy", "totlanduse") 

wk <- unique(subset(ef, HAName %in% targets & keep)[keep])

write.csv(wk, file = "Sites_1 sampler_per_HydroArea.csv")

