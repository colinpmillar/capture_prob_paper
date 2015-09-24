

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

library(ef)

# load appropriate data
load("intermediate_rData/phi.rData") # phi
load("intermediate_rData/phi_new.rData") # phi
load("intermediate_rData/screenedData.rData") # ef

# load final model?
finalf <-  n ~ LifeStage + Trust + fyear + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage
                   
# calculate nicer covariates
ef <- within(ef, {
              cDistance_s = c(scale(Distance_s))
              cWater_W = c(scale(Water_W))
              cElevation_ = c(scale(Elevation_))
              fyear = factor(fyear)
              Trust = factor(Trust)
              LifeStage = factor(LifeStage)
              })
contrasts(ef $ fyear) <- "contr.sum"
contrasts(ef $ Trust) <- "contr.treatment"
#levels(ef $ Trust) <- c("Tweed", )

final <- efp(finalf, data = ef, pass = pass)


finalf_nopass <-  n ~ LifeStage + Trust + fyear + cWater_W + 
                   cElevation_ + cDistance_s + 
                   LifeStage + s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage


final_nopass <- efp(finalf_nopass, data = ef, pass = pass)

mean((fitted(final_nopass) - fitted(final)) / fitted(final))


A1 <- getASummary(final, ef) $ A
A2 <- getASummary(final_nopass, ef) $ A

mean(A2 - A1) / A1)





# ------------------------------------------------
# ------------------------------------------------
# 
#  Simulation of density errors
# 
# ------------------------------------------------
# ------------------------------------------------

# estimate density
library(dplyr)
library(tidyr)

getp <- function(x, data) {
  data $ pijk <- x
  tmp <- data %>% 
         select(sampleID, LifeStage, pass, pijk) %>% 
         as.data.frame(.) %>%
         spread(pass, pijk)

  1 - (1-tmp[["1"]]) * (1-tmp[["2"]]) * (1-tmp[["3"]])
}

getn <- function(data) {
  c(tapply(data $ n, list(data $ LifeStage, data $ sampleID), sum))
}

simp <- function(model, data, nsim = 1000) {
  simb <- MASS::mvrnorm(nsim, coef(model), model $ Vb * 3.58)
  simpijk <- 1/(1 + exp(-model $ Gsetup $ X %*% t(simb)))

  tmp1 <- data %>% select(sampleID, LifeStage, pass) %>% as.data.frame(.)
  apply(simpijk, 2, function(p) {
    tmp1 $ pijk <- p
    tmp2 <- tmp1 %>% spread(pass, pijk)

    1 - (1-tmp2[["1"]]) * (1-tmp2[["2"]]) * (1-tmp2[["3"]])
  })
}

getASummary <- function(model, data, nsim = 1000) {
  ps <- simp(model, data, nsim = nsim)

  tmp <- data %>% 
         select(sampleID, LifeStage, pass, n) %>% 
         as.data.frame(.) %>%
         mutate(pass = paste0("n_", pass)) %>%
         spread(pass, n)

  within(tmp, {
      n = getn(data)
      p = getp(fitted(model), data)
      A = n/p
      varA = n * (1-p)/p^2 + n^2 * apply(1/ps, 1, var)
      cvA = sqrt(varA) / A
    })
}

# get models for summaries

modelledpf <-  n ~ LifeStage + Trust + fyear + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage

constantpf <- n ~ LifeStage

samplewisepf <- n ~ LifeStage * factor(sampleID)

# fit models and summaries variances

modelledp <- efp(modelledpf, data = ef, pass = pass)
modelledp.sum <- getASummary(modelledp, ef)
head(modelledp.sum)

constantp <- efp(constantpf, data = ef, pass = pass)
constantp.sum <- getASummary(constantp, ef)
head(constantp.sum)

plot(constantp.sum $ cvA, modelledp.sum $ cvA)
abline(b=1,a=0)
plot(log(constantp.sum $ cvA), log(modelledp.sum $ cvA))
abline(b=1,a=0)

modelledp.sum[which(log(modelledp.sum$cvA) > 2),]
constantp.sum[which(log(modelledp.sum$cvA) > 2),]
subset(as.data.frame(ef), sampleID == 3147)
table(ef $ Trust)

#samplewisep <- efp(samplewisepf, data = ef, pass = pass)
# get samplewise CVs
if (FALSE) {
  n <- nrow(ef) # data points
  N <- sum(tapply(ef$n, ef$sampleID, sum) > 0)
  sIDs <- sort(unique(ef $ sampleID))
  samplewisep.sum <- constantp.sum
  for (i in 1:N) {
    if (i%%10==0) cat("done", i, "of", N, "     \r"); flush.console()
    samp <- subset(ef, sampleID == sIDs[i])
    if (any(tapply(samp $ n, samp $ LifeStage, sum) == 0)) {
      mod <- efp(n ~ 1, data = samp, pass = pass, hessian = TRUE)
    } else {
      mod <- efp(n ~ LifeStage, data = samp, pass = pass, hessian = TRUE)
    }
    out <- try(getASummary(mod, samp))
    if (class(out) == "try-error") {
      samplewisep.sum[samplewisep.sum$sampleID == sIDs[i],c("cvA", "varA", "A", "p")] <- NA 
    } else {
      samplewisep.sum[samplewisep.sum$sampleID == sIDs[i],] <- out
    }
  }
  cat("Done!\n")
  save(samplewisep.sum, file = "intermediate_rData/samplewisep.sum.rData")
}
load("intermediate_rData/samplewisep.sum.rData")


# ------------------------------------------------
# ------------------------------------------------
# 
#  Violin plots of CVs
# 
# ------------------------------------------------
# ------------------------------------------------

library(lattice)
cv.sum <- data.frame(model = rep(c("Constant p", "Sample-wise p", "Modelled p"), each = nrow(modelledp.sum)),
                     cv = c(constantp.sum $ cvA, samplewisep.sum $ cvA, modelledp.sum $ cvA),
                     LifeStage = rep(constantp.sum $ LifeStage, 3))

x_at <- 10^c(-2:2)
x_labels <- sprintf("%.*f", pmax(-1 * log10(x_at*100), 0), x_at * 100)
#x_labels <- paste("10^", log10(x_at))
p1 <- lattice::bwplot(model ~ cv | LifeStage, 
                data = cv.sum, 
                panel = function(..., box.ratio) {
                  #panel.abline(v = log10(x_at), col = "grey")
                  panel.violin(..., col = "lightblue", box.ratio=box.ratio)
                  panel.bwplot(..., col = 'black', pch='|', box.ratio=.1)
                },
                par.settings = list(box.rectangle = list(col="black", fill = "grey"),
                                    box.umbrella  = list(col="black"),
                                    plot.symbol   = list(col="black", pch=3, cex=0.5)),
                xlab='% CV', strip.left = TRUE, strip = FALSE,
                outer=TRUE, horizontal=TRUE, 
                xlim = c(0.001, 10^3), layout = c(1,2),
                scales=list(x=list(log=10, at=x_at, labels=x_labels, tck=c(1,0))))
print(p1)

png("figures/violinplots.png", width = 7, height = 5, res = 600, units = "in")
print(p1)
dev.off()




# ------------------------------------------------
# ------------------------------------------------
# 
#  Maps of spatial bias
# 
# ------------------------------------------------
# ------------------------------------------------


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
redctm <- redctm[redctm $ CATCH_ID %in% ef $ CATCH_ID,]
hma <- hma[hma $ HACode %in% ef $ HACode,]

hma @ bbox <- bbox(redctm)
hma @ bbox[4] <- 1100000



alldat <- ef[c(names(pdata0), "NEAR_X", "NEAR_Y", "sampleID")]
alldat <- subset(alldat, pass23 == "1")[c("NEAR_X", "NEAR_Y", "sampleID", "LifeStage")]

## values for modelled P
alldat <- merge(alldat, modelledp.sum[c("sampleID", "LifeStage","A")])
alldat <- rename(alldat, fit1 = A)

## values for constant p
alldat <- merge(alldat, constantp.sum[c("sampleID", "LifeStage","A")])
alldat <- rename(alldat, fit2 = A)


## values for saturated p
alldat <- merge(alldat, samplewisep.sum[c("sampleID", "LifeStage","A")])
alldat <- rename(alldat, fit3 = A)


alldat $ fitcp <- alldat $ fit2 / alldat $ fit1
alldat $ fitsp <- alldat $ fit3 / alldat $ fit1


library(hexbin)
# make a grid of predictions

alldat_fry <- subset(alldat, LifeStage == "Fry")
alldat_parr <- subset(alldat, LifeStage == "Parr")

# set up bins
dat <- hexbin(alldat_fry $ NEAR_X, alldat_fry $ NEAR_Y, xbins = 40, IDs = TRUE, shape = 1.5)
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
png("figures/densityBiasMaps.png", width = 9, height = 11, res = 500, units = "in")
par(mfrow = c(2,2), mar = c(0,0,0,0))

coast @ bbox <- bbox(redctm)
coast @ bbox[3] <- 440000
coast @ bbox[4] <- 1100000

# set up breaks
breaks1 <- seq(0.70, 1.10, by = 0.02)
nbreaks1 <- length(breaks1) - 1

breaks <- breaks1
nbreaks <- nbreaks1

# set up colours
cols <- colorRampPalette(c(rich.colors(5)[5], rich.colors(5)[4], "white", rich.colors(5)[2], rich.colors(5)[1]))(nbreaks)
cols[length(cols)/2+0.5] <- rich.colors(15)[8]

# assign colours to cells
means <- tapply(alldat_fry $ fitcp, dat@cID, function(x) mean(x[!is.nan(x)], na.rm=TRUE))
range(means, na.rm = TRUE)
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
# tag on high values
breaks2 <- c(max(breaks1), 1.25, 1.5, 2, 3, 4)
nbreaks2 <- length(breaks2) - 1

breaks <- c(breaks1, breaks2[-1])
nbreaks <- length(breaks) - 1

# set up colours
cols1 <- colorRampPalette(c(rich.colors(5)[5], rich.colors(5)[4], "white", rich.colors(5)[2], rich.colors(5)[1]))(nbreaks1)
cols1[length(cols1)/2+0.5] <- rich.colors(15)[8]

cols2 <- colorRampPalette(c("purple", "hotpink"))(nbreaks2)

cols <- c(cols1, cols2)


# assign colours to cells
alldat_fry $ fitsp_noInf <- alldat_fry $ fitsp
alldat_fry $ fitsp_noInf[alldat_fry $ fitsp_noInf > 1e6 | is.nan(alldat_fry $ fitsp_noInf)] <- NA
means <- tapply(alldat_fry $ fitsp_noInf, dat@cID, function(x) mean(x[!is.nan(x)], na.rm=TRUE))
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


mtext(c("a","b"), side = 3, line = -6, outer = TRUE, font = 2, at = c(0, 0.5) + 0.05)


# set up breaks
breaks <- breaks1
nbreaks <- nbreaks1

# set up colours
cols <- colorRampPalette(c(rich.colors(5)[5], rich.colors(5)[4], "white", rich.colors(5)[2], rich.colors(5)[1]))(nbreaks)
cols[length(cols)/2+0.5] <- rich.colors(15)[8]

# assign colours to cells
means <- tapply(alldat_parr $ fitcp, dat@cID, function(x) mean(x[!is.nan(x)], na.rm=TRUE))
range(means, na.rm = TRUE)
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
# tag on high values
breaks2 <- c(max(breaks1), 1.25, 1.5, 2, 3, 4)
nbreaks2 <- length(breaks2) - 1

breaks <- c(breaks1, breaks2[-1])
nbreaks <- length(breaks) - 1

# set up colours
cols1 <- colorRampPalette(c(rich.colors(5)[5], rich.colors(5)[4], "white", rich.colors(5)[2], rich.colors(5)[1]))(nbreaks1)
cols1[length(cols1)/2+0.5] <- rich.colors(15)[8]

cols2 <- colorRampPalette(c("purple", "hotpink"))(nbreaks2)

cols <- c(cols1, cols2)


# assign colours to cells
alldat_parr $ fitsp_noInf <- alldat_parr $ fitsp
alldat_parr $ fitsp_noInf[alldat_parr $ fitsp_noInf > 1e6 | is.nan(alldat_parr $ fitsp_noInf)] <- NA
means <- tapply(alldat_parr $ fitsp_noInf, dat@cID, function(x) mean(x[!is.nan(x)], na.rm=TRUE))
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


mtext(c("c","d"), side = 3, line = -42, outer = TRUE, font = 2, at = c(0, 0.5) + 0.05)


dev.off()
}


























# ------------------------------------------------
# ------------------------------------------------
# 
#  a plot of data coverage
# 
# ------------------------------------------------
# ------------------------------------------------


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


