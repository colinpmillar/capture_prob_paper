
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


#tmp <- hma
#tmp $ Trust <- "MSS"
#X <- predict(g1, type = "lpmatrix", newdata = tmp)
#b <- coef(g1)
#tmp $ val <- c(X %*% b)

#tmp $ col <- cut(tmp $ val, 6)
#cols <- heat.colors(6)
#plot(tmp, col = cols[tmp $ col], main = i)

# strip off shetland and orkney
hma <- hma[!hma $ HACode %in% c(108, 107),]

alldat <- pdata0
alldat[["HACode"]] <- hma $ HACode
alldat <- do.call(data.frame, c(alldat, list(stringsAsFactors = FALSE)))
alldat[c("fit", "se")] <- predict(g1, newdata = alldat, se = TRUE)
alldat $ colgrp <- as.numeric(cut(1/(1 + exp(-alldat $ fit)), 11))
alldat $ cols <- gplots::rich.colors(11)[12 - alldat $ colgrp]

hma $ cols <- gplots::rich.colors(11)[12 - alldat $ colgrp]

png(file = "figures/pmodel_HA_map.png", width = 9, height = 9, units = "in", res = 500)

plot(hma, col = hma $ cols, ylim = c(550000, 970000), xlim = c(5540, 480000))


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

