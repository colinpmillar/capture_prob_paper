
if (Sys.info()["user"] == "millaco") {
  setwd("~/work/SMFS-report")    
} else 
if (Sys.info()["user"] == "millarc") {
  setwd("B:/Conservation_Limits/CL_Juvenile_Density/SMFS-report")
} else 
if (Sys.info()["user"] == "Millarc") {
  setwd("C:/work/SMFS-report")
}




#############################################################
#   
#   map of trust capture probabilities
#
#############################################################

# load fits and model data
load("rData/densmodelData.rData")
load("rData/bestpmodel.rData")
library(CLdata)

var.summary <- best $ Gsetup $ var.summary
var.type <- sapply(var.summary, function(x) is(x)[1])

# set base prediction levels
pdata0 <- ifelse(var.type == "numeric", 
                   lapply(var.summary, "[", 2), 
                   lapply(var.summary, function(x) levels(x)[floor(nlevels(x)/2)])
                )
pdata0 $ Trust <- "MSS"
pdata0 $ fyear <- "2006"


alldat <- pdata0
alldat[["Trust"]] <- levels(ef $ Trust)
alldat <- do.call(data.frame, c(alldat, list(stringsAsFactors = FALSE)))
alldat[c("fit", "se")] <- predict(g1, newdata = alldat, se = TRUE)
alldat $ colgrp <- as.numeric(cut(1/(1 + exp(-alldat $ fit)), 11))
alldat $ cols <- gplots::rich.colors(11)[12 - alldat $ colgrp]


trustpoly <- rgdal::readOGR("../GIS_shapefiles_for_model/District maps data","sa_trust_1110")
trustnames <- read.csv("trust_names.csv")
trustpoly @ data <- merge(trustpoly @ data, trustnames, all.x = TRUE)
trustpoly @ data <- merge(trustpoly @ data, alldat[c("Trust", "fit", "cols")], sort = FALSE, all.x = TRUE)
trustpoly @ data <- trustpoly @ data[order(trustpoly @ data $ TrustCode),]
rownames(trustpoly @ data) <- 1:length(trustpoly) - 1


png(file = "figures/pmodel_trust_map.png", width = 9, height = 9, units = "in", res = 500)

plot(trustpoly, col = trustpoly $ cols, ylim = c(550000, 970000), xlim = c(5540, 480000))
plot(redctm[grep("Naver", redctm @ data $ DESCRIPTIO), ], col = subset(alldat, Trust == "Naver") $ cols, add = TRUE)

which <- grep("(Forss)|(Thurso)|(Wester)|(Wick)|(Dunbeath)|(Berriedale)|(Langwell)", redctm @ data $ DESCRIPTIO)
otherpoly <- rgeos::gUnaryUnion(redctm[which, ])
plot(otherpoly, col = subset(alldat, Trust == "Other") $ cols, add = TRUE)

# MSS etc polygons
x0 <- 20978; y0 <- 620000; dx <- 42000; dy <- 30000


j <- 0
for (i in c("MSS", "SEPA")) {
polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - j*dy + c(0, -dy)[c(1,1,2,2)], col = subset(alldat, Trust == i) $ cols)
text(x0 + 1.1 * dx, y0 - j*dy - 0.5*dy, label = i, font = 2, cex = 0.8, adj = 0)
j <- j + 1
}

# now for legend
x0 <- 421799; y0 <- 947168; dx <- 42000; dy <- 30000
vals <- as.numeric(unique(unlist(strsplit(gsub("[(]|[]]", "", levels(cut(1/(1 + exp(-alldat $ fit)), 11))), ","))))
for (i in 0:10) {
polygon(x0 + c(0, dx)[c(1,2,2,1)], y0 - i*dy + c(0, -dy)[c(1,1,2,2)], col = gplots::rich.colors(11)[i+1])
}
text(x0 + 1.1*dx, y0 - 0:11*dy, sprintf("%.2f", rev(vals)), font = 2, cex = 0.8, adj = 0)


dev.off()


#colgrp <- as.numeric(cut(trustpoly @ data $ se, 11))
#cols <- rich.colors(11)[12 - colgrp]

#plot(trustpoly, col = cols)

