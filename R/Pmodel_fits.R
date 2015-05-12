

library(CLdata)
library(CLmodel)

library(sp)
library(spdep)
library(rstan)


if (Sys.info()["user"] == "millaco") {
  setwd("~/work/SMFS-report")    
} else 
if (Sys.info()["user"] == "millarc") {
  setwd("B:/Conservation_Limits/CL_Juvenile_Density/SMFS-report")
} else 
if (Sys.info()["user"] == "Millarc") {
  setwd("C:/work/SMFS-report")
}

# load data
load("rData/modelData.rData")


# compile model
m0 <- efp(X ~ 1, data = ef, passes = "Runs")

######################################
##
##   MODELLING Setup
##
######################################


## Q model for spatial regions

{
  hmaadj <- poly2nb(hma, queen = FALSE)
  hmaadj <- nb2mat(hmaadj, style = "B", zero.policy = TRUE)
  # add connections for inner and outer hebs
  hmaadj[hma $ HAName == "Inner Hebrides", hma $ HAName == "Outer Hebrides"] <- 1
  hmaadj[hma $ HAName == "Outer Hebrides", hma $ HAName == "Inner Hebrides"] <- 1
  hmaadj[hma $ HAName == "Inner Hebrides", c(21, 42, 43)] <- 1
  hmaadj[c(21, 42, 43), hma $ HAName == "Inner Hebrides"] <- 1
  hmaadj <- as(hmaadj, "dgTMatrix")

  Q <- hmaadj
  if (any(Q @ x == 0)) stop("something went wrong")
  Q @ x[] <- -1/Q @ x
  diag(Q) <- rowSums(hmaadj)
  Q <- as.matrix(Q)
  colnames(Q) <- rownames(Q) <- hma $ HACode
  Qhma <- Q
}


par(mfrow = c(3,3), mar = c(0,0,3,0)) 
for (i in 1:9 * 2 + 2) {

m0 <- efp(X ~ s(HACode, k = i, bs = "gmrf", xt = list(penalty = Qhma)), 
          data = ef, passes = "Runs", hessian = TRUE, verbose = FALSE)


g1 <- gam(G = m0 $ Gsetup)
qr.G <- qr(m0 $ G)
rank.deficient <- qr.G $ pivot[abs(diag(qr.G $ qr)) < 1e-7]
whichkeep <- -rank.deficient
if (!length(whichkeep)) whichkeep <- 1:length(m0 $ coefficients) 
g1 $ coefficients[] <- 0
g1 $ coefficients[whichkeep] <- m0 $ coefficients       
g1 $ Vp[] <- 0
diag(g1 $ Vp[]) <- 1e-5
g1 $ Vp[whichkeep, whichkeep] <- m0 $ Vb
g1 $ family <- binomial()

tmp <- hma
tmp $ Trust <- "MSS"
X <- predict(g1, type = "lpmatrix", newdata = tmp)
b <- coef(g1)
tmp $ val <- c(X %*% b)

tmp $ col <- cut(tmp $ val, 6)
cols <- heat.colors(6)
plot(tmp, col = cols[tmp $ col], main = i)

}







#############################################################
#   
#   Step 1
#
#############################################################

# now some model selection
f1s <- c("Trust",
         "s(HACode, k = 8, bs = 'gmrf', xt = list(penalty = Qhma))",
         "fyear", 
         "poly(Water_W, 1)",
         "poly(Elevation_, 1)", 
         "poly(Distance_s, 1)",
         "poly(sinSlope, 1)",
         "poly(Upcatch_km, 1)",
         "poly(Urban, 1)",
         "poly(NCTrees, 1)",
         "poly(CTrees, 1)",
         "poly(Mixed, 1)",
         "poly(Marsh, 1)",         
         "poly(Other, 1)",
         "poly(doy, 1)",
         "s(Water_W, k = 3)",
         "s(Elevation_, k = 3)", 
         "s(Distance_s, k = 3)",
         "s(sinSlope, k = 3)",
         "s(Upcatch_km, k = 3)",
         "s(Urban, k = 3)",
         "s(NCTrees, k = 3)",
         "s(CTrees, k = 3)",
         "s(Mixed, k = 3)",
         "s(Marsh, k = 3)",         
         "s(Other, k = 3)",
         "s(doy, k = 3)"
         )

# should maybe return the best model as well summaries...
runModels <- function(chosen, ...) {
  # build a formula list of additions
  formsadd <- lapply(f1s[!chosen], function(x) as.formula(paste0("X ~ ", paste(c(f1s[chosen], x), collapse = " + "))))

  descadd <- f1s[!chosen]
  dropadd <- rep("add", sum(!chosen))

  #build a formula list dropping elements
  if (any(chosen)) {
    if (sum(chosen)==1) {
      formsdrop <- list(X ~ 1)
    } else {
      formsdrop <- lapply(seq(sum(chosen)), function(i) as.formula(paste0("X ~ ", paste(f1s[chosen][-i], collapse = " + ") )))
    }
    descdrop <- f1s[chosen]
    forms <- c(formsadd, formsdrop)
    desc <- c(descadd, descdrop)
    dropadd <- c(dropadd, rep("drop", sum(chosen)))
  } else {
    forms <- formsadd
    desc <- descadd
  }

  mods <- lapply(forms, efp, data = ef, passes = "Runs", ...)

  if (all(!chosen)) {
    m0 <- efp(X ~ 1, data = ef, passes = "Runs", ...)
  } else {
    m0 <- efp(as.formula(paste0("X ~ ", paste(f1s[chosen], collapse = " + "))), data = ef, passes = "Runs", ...)    
  }
  
  tab <- cbind(what = desc, step = dropadd, summaryMods(mods, m0 = m0, order = FALSE)[,-1] )
  tab[order(tab $ Daic),]
}


#############################################################
#   
#   Loop over adding in covariates
#
#############################################################


chosen <- rep(FALSE, length(f1s))
out <- list(Daic = -1)
tol <- 0
outlist <- list() # grow this - bad!
i <- 0
while(out $ Daic[1] < tol & !all(chosen)) {
  out <- runModels(chosen, verbose = FALSE)
  print(out, digits = 3)
  if ( out $ Daic[1] < 0) {
    # then model is an improvement
    chosen[!chosen][as.numeric(rownames(out)[1])] <- TRUE
    print(f1s[chosen])
    outlist[[i <- i + 1]] <- out
  } 
}


f1s[chosen]

# remove  linear terms if there are non linear terms
finalfs <- f1s[chosen]
paste0("X ~ ", paste(finalfs, collapse = " + ")) 
#X ~ Trust + s(HACode, k = 8, bs = 'gmrf', xt = list(penalty = Qhma)) + 
#    fyear + poly(Water_W, 1) + poly(Distance_s, 1) + 
#    s(sinSlope, k = 3) + s(doy, k = 3)


# drop one at a time
forms <- lapply(seq_along(finalfs), function(i) as.formula(paste0("X ~ ", paste(finalfs[-i], collapse = " + ") )))
mods <- lapply(forms, efp, data = ef, passes = "Runs", verbose = FALSE)

m0 <- efp(as.formula(paste0("X ~ ", paste(finalfs, collapse = " + "))), data = ef, passes = "Runs")    
  
tab <- cbind(dropped = finalfs, summaryMods(mods, m0 = m0, order = FALSE)[,-1] ) 
tab[order(tab $ Daic, decreasing = TRUE),]
#                                                   dropped      aic        Daic
#1                                                    Trust 360394.6 889.2142996
#7                                            s(doy, k = 3) 359943.5 438.1212992
#3                                                    fyear 359712.8 207.3613155
#2 s(HACode, k = 8, bs = 'gmrf', xt = list(penalty = Qhma)) 359681.6 176.1508722
#5                                      poly(Distance_s, 1) 359562.8  57.3382802
#4                                         poly(Water_W, 1) 359533.1  27.6788603
#6                                       s(sinSlope, k = 3) 359506.0   0.5996644


best <- efp(X ~ Trust + s(HACode, k = 8, bs = 'gmrf', xt = list(penalty = Qhma)) + 
                fyear + poly(Water_W, 1) + poly(Distance_s, 1) + 
                s(sinSlope, k = 3) + s(doy, k = 3), 
            data = ef, passes = "Runs", verbose = FALSE, hessian = TRUE)    


base <- efp(X ~ 1, 
            data = ef, passes = "Runs", verbose = FALSE, hessian = TRUE)    


summary(best)
best

BIC(best)

BIC(efp(X ~ Trust + s(HACode, k = 10, bs = 'gmrf', xt = list(penalty = Qhma)) + 
                fyear + poly(Water_W, 1) + poly(Distance_s, 1) + 
                s(sinSlope, k = 3) + s(doy, k = 3), 
            data = ef, passes = "Runs", verbose = FALSE))

# k = 10 is slightly better than k = 8 but not enough to worry about here.


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
setdiff(levels(ef $ fyear[drop=TRUE]), levels(g1 $ model[,"fyear"]))
setdiff(levels(ef $ Trust), levels(g1 $ model $ Trust))

X <- predict(g1, type = "lpmatrix", newdata = ef)

ef $ p <- 1/(1 + exp(-X %*% coef(g1)))
ef $ logitp.se <- sqrt(diag(X %*% g1 $ Vp %*% t(X)))

ef $ offset <- with(ef, log( (1-(1-p)^Runs) * Area) )
psim <- matrix(rnorm(1000 * nrow(ef)), ncol = 1000, nrow = nrow(ef))
psim <- psim * ef $ logitp.se + c(X %*% coef(g1))
psim <- 1/(1 + exp(-psim))
offsetsim <- with(ef, (1-(1-psim)^Runs) * Area )
ef $ offset.se <- apply(offsetsim, 1, sd)
ef $ weights <- 1/ef $ offset.se^2

# add in constant p

ef $ constantp <- 1/(1 + exp(-coef(base)))


save(best, g1, file = "rData/bestpmodel.rData")
save(ef, file = "rData/densmodelData.rData")

