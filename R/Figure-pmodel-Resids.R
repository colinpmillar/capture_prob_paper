

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



# residuals
ef4 <- ef3

# get p predictions for ef
ef4 $ pbig <- fitted(best)

# get logLik components
ef4 $ R <- with(ef4, s - 1 - Z)
ef4 $ llsat <- with(ef4, T * log(psat) + T * R * log(1-psat) - T * log(1 - (1-psat)^s) )
ef4 $ llbig <- with(ef4, T * log(pbig) + T * R * log(1-pbig) - T * log(1 - (1-pbig)^s) )

# calculate Deviance components and residuals
ef4 $ devcomp <- with(ef4, 2 * (llsat - llbig))
ef4 $ devcomp[abs(ef4 $ devcomp) < 1e-9] <- 0
ef4 $ devres <- with(ef4, sign(psat - pbig)*sqrt(devcomp))
ef4 $ devres.scaled <- ef4 $ devres / 2.13
plot(ef4 $ devres)
hist(ef4 $ devres)
qqnorm(ef4 $ devres)

# plot
library(lattice)
p1 <- xyplot(devres.scaled ~ I(year + doy/365) | Trust,
             data = ef4, type = c("p", "g"), pch = 16, cex = 0.5)

p1

png(file = "figures/pmodel_resids.png", width = 7, height = 7, units = "in", res = 500)
print(p1)
dev.off()
