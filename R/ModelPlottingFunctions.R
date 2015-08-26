
# set up prediction data frame
getPlotData <- function(var, func = function(x) exp(x) / (1+exp(x)), model = g1, bootb = NULL) {
  args <- pdata0[!names(pdata0) %in% var]
  args[var] <- pdata1[var]

  pdata <- do.call(expand.grid, args)

  # and predict
  pdata[c("fit", "se")] <- predict(model, newdata = pdata, se.fit = TRUE)
  pdata $ p <- func(pdata $ fit)
  pdata $ cil <- func(pdata $ fit - 2*pdata $ se * sqrt(phi))
  pdata $ ciu <- func(pdata $ fit + 2*pdata $ se * sqrt(phi))
  pdata $ var <- paste(var, collapse = ":")

  # use bootstrap / simulated model params if present
  if (!is.null(bootb)) {
    X <- predict(model, type = "lpmatrix", newdata = pdata)
    keep <- grep(paste(var, collapse = "|"), colnames(X))
    bootp <- X[,keep,drop=FALSE] %*% t(bootb[,keep, drop = FALSE])
    bmean <- apply(bootp, 1, mean)
    bcil <- apply(bootp, 1, quantile, 0.025)
    bciu <- apply(bootp, 1, quantile, 0.975)
    pdata $ cil <- func(pdata $ fit + bmean - bciu)
    pdata $ ciu <- func(pdata $ fit + bmean - bcil)
  }


  if (length(var) > 1) {
    #stop()
    #pdata $ x <- as.numeric(do.call(interaction, c(pdata[var], list(sep = " "))))
    pdata $ x <- as.numeric(pdata[[var[1]]])
    pdata $ cx <- paste(pdata[[var[1]]]) 
  } else {
    pdata $ x <- as.numeric(pdata[[var]])
    pdata $ cx <- paste(pdata[[var]])    
  }
  
  pdata
}



# define the plotting functions


factorPlot <- function(pdata, xlab = "", ylim = NULL, yaxislab = TRUE, labcex = 0.9) {

  if (is.null(ylim)) {
    ylim <- range(pdata $ cil, pdata $ ciu)
  }

  plot(pdata $ p, pch = 16, ylim = ylim, axes = FALSE, ann = FALSE)
  segments(x0 = 1:nrow(pdata), y0 = ylim[1], y1 = pdata $ p, col = grey(0.9), lty = 2)
  segments(x0 = 1:nrow(pdata), y0 = pdata $ cil, y1 = pdata $ ciu)
  axis(2, las = 1, labels = yaxislab)
  axis(1, at = 1:nrow(pdata), labels = FALSE, tck = FALSE)
  text(1:nrow(pdata), ylim[1] - diff(ylim)*.08,  pdata[[pdata $ var[1]]], srt=45, xpd = TRUE, adj = 1, cex = labcex)
  title(ylab = "", main = xlab)
  box(bty = "l")
}


continuousPlot <- function(pdata, ylim = NULL, xlab = "", rug = NULL, yaxislab = TRUE) {

  if (is.null(ylim)) {
    ylim <- range(pdata $ cil, pdata $ ciu)
  } 

  plot(pdata $ x, pdata $ p, type = "n", ylim = ylim, axes = FALSE, ann = FALSE)
#  lines(pdata $ x, pdata $ cil, lty = 2)
#  lines(pdata $ x, pdata $ ciu, lty = 2)
  with(pdata, polygon(c(x, rev(x)), c(cil, rev(ciu)), border = NA, col = "lightblue"))
  lines(pdata $ x, pdata $ p)
  axis(2, las = 1, labels = yaxislab)
  axis(1)
  title(ylab = "", main = xlab)
  box(bty = "l")
  # a rug!
  if (!is.null(rug)) {
    dx <- diff(rx <- range(rug, na.rm = TRUE))
    breaks <- seq.int(rx[1], rx[2], length.out = 101)

    dens <- table(cut(rug, breaks))
    cols <- heat.colors(11)

    dz <- diff(rz <- range(dens[dens > 0], na.rm = TRUE))
    zbreaks <- seq.int(rz[1], rz[2], length.out = 11)

    zs <- as.numeric(cut(dens, zbreaks))
    for (i in 1:100) {
      polygon(breaks[i + 0:1][c(1,2,2,1)], ylim[1] + c(0, diff(ylim)*.05)[c(1,1,2,2)], 
              col = cols[zs[i]], border = NA)
    }
  }
}


continuousPlot2 <- function(pdata, ylim = NULL, xlab = "", rug = NULL, yaxislab = TRUE) {

  if (is.null(ylim)) {
    ylim <- range(pdata $ cil, pdata $ ciu)
  } 

  plot(pdata $ x, pdata $ p, type = "n", ylim = ylim, axes = FALSE, ann = FALSE)
#  lines(pdata $ x, pdata $ cil, lty = 2)
#  lines(pdata $ x, pdata $ ciu, lty = 2)

  by <- strsplit(pdata $ var, ":")[[1]][2]
  cols <- c("blue", "red")
  cols <- colorRampPalette(cols)(2)
  colsp <- paste0(cols, "33")
  for (i in 1:2) {
    with(subset(pdata, pdata[[by]] == unique(pdata[[by]])[i]), {
      polygon(c(x, rev(x)), c(cil, rev(ciu)), border = NA, col = colsp[i])
    })
  }
  for (i in 1:2) {
    with(subset(pdata, pdata[[by]] == unique(pdata[[by]])[i]), {
      lines(x, p, col = cols[i])
    })
  }

  axis(2, las = 1, labels = yaxislab)
  axis(1)
  title(ylab = "", main = xlab)
  box(bty = "l")
  # a rug!
  if (!is.null(rug)) {
    dx <- diff(rx <- range(rug, na.rm = TRUE))
    breaks <- seq.int(rx[1], rx[2], length.out = 101)

    dens <- table(cut(rug, breaks))
    cols <- heat.colors(11)

    dz <- diff(rz <- range(dens[dens > 0], na.rm = TRUE))
    zbreaks <- seq.int(rz[1], rz[2], length.out = 11)

    zs <- as.numeric(cut(dens, zbreaks))
    for (i in 1:100) {
      polygon(breaks[i + 0:1][c(1,2,2,1)], ylim[1] + c(0, diff(ylim)*.05)[c(1,1,2,2)], 
              col = cols[zs[i]], border = NA)
    }
  }
}