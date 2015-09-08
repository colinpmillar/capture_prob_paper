


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




# ------------------------------------------------
# 
#  Model prep
# 
# ------------------------------------------------

best <- efp(X ~ LifeStage + Trust + fyear + poly(Water_W, 1) + poly(Elevation_, 1) +
            poly(Distance_s, 1) + s(doy, k = 3, by = LifeStage), 
            data = ef3, passes = "Runs", hessian = TRUE)    


#summary(best)
best
getScale(ef3, best)
# constant p model
base <- efp(X ~ 1, data = ef3, passes = "Runs", hessian = TRUE)
#base <- efp(X ~ LifeStage, data = ef3, passes = "Runs", hessian = TRUE)


# predict p for data
# get a gam container
g1 <- gam(G = best $ Gsetup)
qr.G <- qr(best $ G)
rank.deficient <- qr.G $ pivot[abs(diag(qr.G $ qr)) < 1e-7]
whichkeep <- -rank.deficient
if (!length(whichkeep)) whichkeep <- 1:length(best $ coefficients) 
names(g1 $ coefficients[-1 * whichkeep])
g1 $ coefficients[] <- 0
g1 $ coefficients[whichkeep] <- best $ coefficients       
g1 $ Vp[] <- 0
diag(g1 $ Vp[]) <- 1e-5
g1 $ Vp[whichkeep, whichkeep] <- best $ Vb
g1 $ family <- binomial()
var.summary <- best $ Gsetup $ var.summary


## check factors
setdiff(levels(ef3 $ Trust), levels(g1 $ model $ Trust))

X <- predict(g1, type = "lpmatrix", newdata = ef3)

ef3 $ p <- 1/(1 + exp(-X %*% coef(g1)))
ef3 $ logitp.se <- sqrt(diag(X %*% g1 $ Vp %*% t(X)))

ef3 $ offset <- with(ef3, log( (1-(1-p)^Runs) * Area) )

psim <- matrix(rnorm(1000 * nrow(ef3)), ncol = 1000, nrow = nrow(ef3))
psim <- psim * ef3 $ logitp.se + c(X %*% coef(g1))
psim <- 1/(1 + exp(-psim))
offsetsim <- with(ef3, (1-(1-psim)^Runs) * Area )
ef3 $ offset.se <- apply(offsetsim, 1, sd)
ef3 $ weights <- 1/ef3 $ offset.se^2

# add in constant p

ef3 $ constantp <- 1/(1 + exp(-coef(base)))


save(best, g1, file = "rData/bestpmodel.rData")
save(ef3, file = "rData/densmodelData.rData")


# ------------------------------------------------
# 
#  Model summaries
# 
# ------------------------------------------------

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
ef4 $ devres.scaled <- ef4 $ devres / 2.11
plot(ef4 $ devres)
hist(ef4 $ devres)
qqnorm(ef4 $ devres)

# plot
library(lattice)
xyplot(devres.scaled ~ I(year + doy/365) | Trust,
       data = ef4, type = c("p", "g"), pch = 16, cex = 0.5)

xyplot(devres.scaled ~ pbig | Trust, group = grepl("GIR",Site.Name),
       data = ef4, type = c("p", "g"), pch = 16, cex = 0.5)


xyplot(devres.scaled ~ I(year + doy/365) | Site.Name,
       data = subset(ef4, Trust == "MSS"), type = c("p", "g"), pch = 16, cex = 0.5)

xyplot(devres.scaled ~ I(year + doy/365) | Site.Name, groups = year,
       data = subset(ef4, grepl("GIR", Site.Name)), type = c("p", "g"), pch = 16, cex = 0.5)


xyplot(pbig ~ psat | cut(llsat, breaks = 50),
       data = ef4, type = c("p", "g"), pch = 16, cex = 0.5)



cols <- c("Site_OBJECTID","Site.Name", "Dataset", "Width", "Date", "Runs", "Area",
          "Trust", "n_R1", "n_R2", "n_R3","psat", "pbig", "devres", "devres.scaled")

head(ef4[order(ef4 $ devres.scaled, decreasing = TRUE),cols], 20)
head(ef4[order(ef4 $ devres.scaled),cols], 20)

subset(ef4,  totalN > 700)

subset(ef4,  Trust == "Nith")




# nosey at width data

tmp <- unique(ef[!names(ef) %in% c("sampleID","Area","Date","fyear","T","X","Z","phi","year",
                         "doy","LifeStage","pass","pass23","n","Runs","Width","Site.Name")])

tail(tmp[order(tmp$ Water_W),], 20)


# data for Karen

db <- "B:/Conservation_Limits/CL_Juvenile_Density/CollateData/db/scottishEF_CLROAME_v01.sqlite3"
con <- DBI::dbConnect(RSQLite::SQLite(), db)
# get data
ef <- DBI::dbReadTable(con, "ef")
DBI::dbDisconnect(con)

ef <- subset(ef, !Dataset %in% c("sfcc","sepa","caithness"))
ef <- ef[c("Site_OBJECTID", "Site.Name", "Dataset", "Date")]

ef <- ef[!is.na(ef $ Site_OBJECTID),]

# get site info
con <- DBI::dbConnect(RSQLite::SQLite(), db)
gis <- DBI::dbReadTable(con, "gis"); DBI::dbDisconnect(con)

# tag on site info data
ef <- cbind(ef, gis[ef $ Site_OBJECTID,])

uef <- unique(ef[c("Site_OBJECTID", "Site.Name", "Dataset", "NEAR_X", "NEAR_Y")])
dim(uef)

head(uef, 50)
write.csv(uef, file = "sitesForKaren.csv")

