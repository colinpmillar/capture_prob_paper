
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
# ------------------------------------------------
# 
#  Variance estimation
# 
# ------------------------------------------------
# ------------------------------------------------

efpackage <- devtools::as.package("c:/work/repos/faskally/ef")
#devtools::check(efpackage)
#devtools::install(efpackage)
devtools::load_all(efpackage)

# load appropriate data
load("intermediate_rData/phi.rData") # phi
load("intermediate_rData/screenedData.rData") # ef

# load final model?
finalf <-  n ~ LifeStage + Trust + fyear + pass23 + cWater_W + 
                   cElevation_ + cDistance_s + 
                   LifeStage:pass23 + s(doy, k = 3, by = LifeStage) + 
                   cElevation_:LifeStage
                   
ef <- within(ef, {
              cDistance_s = scale(Distance_s)
              cWater_W = scale(Water_W)
              cElevation_ = scale(Elevation_)
              fyear = factor(fyear)
              Trust = factor(Trust)
              })
contrasts(ef $ fyear) <- "contr.sum"
contrasts(ef $ Trust) <- "contr.sum"
#levels(ef $ Trust) <- c("Tweed", )
ef $ LifeStage <- factor(ef $ LifeStage)

final <- efp(finalf, data = ef, pass = pass)




# ------------------------------------------------
# ------------------------------------------------
# 
#  Re-estimate phi
# 
# ------------------------------------------------
# ------------------------------------------------
 

load("intermediate_rData/llsample.rData")

Nsamples <- sum(tapply(ef$n, ef$sampleID, sum) > 0)
phi_new <- c(2 * (sum(llsample) - logLik(final)) / Nsamples)
phi_new
save(phi_new, file = "intermediate_rData/phi_new.rData") # phi_new
# 3.590581
load("intermediate_rData/phi_new.rData")

# ------------------------------------------------
# ------------------------------------------------
# 
#  Bootstrap check of confidence intervals
# 
# ------------------------------------------------
# ------------------------------------------------


# bootstrap
if (FALSE) {
  nboot <- 2000
  bootb <- matrix(NA, nboot, length(coef(final)))
  ef <- ef[order(ef$sampleID),]
  sids <- paste(unique(ef $ sampleID))
  ids <- 1:length(sids)
  names(ids) <- sids
  for (i in 1:nboot) {
    cat("                                        \rrunning set ", i, "of", nboot)
    flush.console()
    bids <- sample(sids, length(sids), replace = TRUE)
    bids <- ids[bids]
    dat <- ef[1:6 + rep((bids-1)*6, each = 6),]
    while (nlevels(dat $ Trust[drop=TRUE]) != nlevels(ef $ Trust[drop=TRUE]) |
           nlevels(dat $ fyear[drop=TRUE]) != nlevels(ef $ fyear[drop=TRUE])) {
      bids <- sample(sids, length(sids), replace = TRUE)
      bids <- ids[bids]
      dat <- ef[1:6 + rep((bids-1)*6, each = 6),]    
    }
    mod <- efp(finalf, data = dat, pass = pass)
    bootb[i,] <- coef(mod)  
  }
  save(bootb, file = "intermediate_rData/bootb.rData") # ef
}
load("intermediate_rData/bootb.rData") # ef


H2 <- var(bootb)
dim(H2)
library(Matrix)
image(Matrix(cov2cor(H2)))
H2



# simple
H1 <- final $ Vb

image(Matrix(cov2cor(H1) - cov2cor(H2)))

cbind(sqrt(diag(H1) * phi_new), sqrt(diag(H2)))
t(t(sqrt(diag(H1) * phi_new) / sqrt(diag(H2))))

tmp <- sqrt(diag(H1) * phi_new) / sqrt(diag(H2))
mean(tmp[grep("pass|LifeStage", names(tmp))])
#[1] 1.182284
mean(tmp[!grepl("pass|LifeStage", names(tmp))])
#[1] 0.8546158
 

plot(cbind(sqrt(diag(H1) * phi_new), sqrt(diag(H2))))

simb <- MASS::mvrnorm(1000, coef(final), H1)
X <- final $ Gsetup $ X
simp <- X %*% t(simb)

tmp <- 
  apply(1-simp, 2, function(p) 
       1 - tapply(p, list(ef$LifeStage, ef$sampleID), prod))
dim(tmp)


cbind(coef(final), sqrt(diag(H1)))

Cov <- cov2cor(H1)
diag(Cov) <- 0
image(Matrix(Cov))
A <- Cov[!grepl("year|Trust", rownames(Cov)), !grepl("year|Trust", rownames(Cov))]
ids <- which(apply(A, 2, function(x) any(abs(x) > 0.5)))
A <- A[ids,ids]
tmp <- round(A, 3)
tmp[lower.tri(tmp)] <- 0
Matrix(tmp)
image(Matrix(A))

