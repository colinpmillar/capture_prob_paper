
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

# One way of assessing the adequacy of a model is to compare it with a more general 
# model with the maximum number of parameters that can be estimated. It is referred to
# as the saturated model. In the saturated model there is basically one parameter per
# observation. The deviance assesses the goodness of fit for the model by looking at the
# difference between the log-likelihood functions of the saturated model and the model
# under investigation, i.e. (b y) (b y) sat l , − l , . Here sat b denotes the maximum likelihood
# estimator of the parameter vector of the saturated model, sat β , and b is the maximum
# likelihood estimator of the parameters of the model under investigation, β . The
# maximum likelihood estimator is the estimator that maximises the likelihood function.
# The deviance is defined as
# D = 2{l(bsat , y) − l(b, y)}.


efpackage <- devtools::as.package("c:/work/repos/faskally/ef")
devtools::load_all(efpackage)

# load data
load("rData/modelData.rData")
# strip off unused columns
ef <- ef[!names(ef) %in% c(paste0("n_R", 4:6), "pDate")]
# keep salmon, 3 pass fishings within the covariate bounds
ef <- subset(ef, Runs > 2 & Species == "Salmon" & keep)
ef <- dplyr::as_data_frame(ef)

strReverse <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

ef $ sampleID <- as.numeric(factor(with(ef, 
                  paste0(Trust,
                         strReverse(Date), 
                         Site_OBJECTID))))

# gather
ef <- tidyr::gather(ef, pass, n, n_R1:n_R3)
ef $ pass <- as.numeric(gsub("n_R", "", ef $ pass))
ef $ pass23 <- factor(replace(ef $ pass, ef $ pass > 2, 2))

# set NAs to zero:  This is beacuse in the SFCC database, if there were no fish
# caught, it was practice to enter a zero for the last pass fished and empty cells
# elswhere, hence a zero for pass 3 and then NAs for passes 1 and 2

ef $ n[is.na(ef $ n)] <- 0

# drop the 2 double samples
ef <- subset(ef, !sampleID %in% names(which(table(ef $ sampleID) == 12)))

# reorder columns
ef <- ef[c(52:ncol(ef),1:51)]
ef

# woodland
ef $ woodland <- with(ef, CTrees + NCTrees + Mixed)

# drop out dodgy sites / samples and covariate outliers
ef <- subset(ef, Site_OBJECTID != 2180)
ef <- subset(ef, !(grepl("GIR", Site.Name) & year < 2006))
ef <- subset(ef, !(Date == "13/09/2002" & Site_OBJECTID == 1527))
ef <- subset(ef, !(Date == "13/09/2002" & Site_OBJECTID == 4004))
ef <- subset(ef, !(Date == "27/09/2004" & Site_OBJECTID == 1955))
ef <- subset(ef, !(Date == "03/08/2006" & Site_OBJECTID == 2626))
ef <- subset(ef, keep)


# some summaries
sum(tapply(ef$n, ef$sampleID, sum) >= 0)
# 4333 samples in total
sum(tapply(ef$n, ef$sampleID, sum) > 0)
# of which 2838 had some fish

length(unique(ef $ Site_OBJECTID))
# 2179 unique sites
length(unique(ef $ Site_OBJECTID[ef $ sampleID %in% names(tapply(ef$n, ef$sampleID, sum)[tapply(ef$n, ef$sampleID, sum) > 0])]))
# 1451 unique sites with non-zero samples


# remove samples with zero counts
ef <- subset(ef, sampleID %in% names(tapply(n, sampleID, sum)[tapply(n, sampleID, sum) > 0]))

# how many samples have one lifestage with zeros?
sum(tapply(ef$n, list(ef$sampleID, ef $ LifeStage), sum) == 0)
# 480 samples observe only one lifestage

# some covariate summaries
length(unique(ef $ Trust))
# 25
length(unique(ef $ HACode))
# 44

range(apply(table(ef $ CATCH_ID, ef $ HACode), 2, function(x) sum(x>0)))
# 1 11

length(sort(unique(ef $ year)))
# 17

# some quick model fits
efp(n ~ 1, data = ef, pass = pass)
efp(n ~ LifeStage, data = ef, pass = pass)
efp(n ~ pass23, data = ef, pass = pass)

invlogit(efp(n ~ pass23 - 1, data = ef, pass = pass) $ coeff)

# fit samplewise saturated model
if (FALSE) {
  n <- nrow(ef) # data points
  N <- n / 6    # site visits (samples)
  llsat <- rep(NA, N)
  sIDs <- sort(unique(ef $ sampleID))
  for (i in 1:N) {
    if (i%%10==0) cat("done", i, "of", N, "     \r"); flush.console()
    samp <- subset(ef, sampleID == sIDs[i])  
    mod <- efp(n ~ factor(pass) * LifeStage, data = samp, pass = pass, hessian = FALSE)
    llsat[i] <- logLik(mod)
  }
  cat("Done!\n")
  save(llsat, file = "intermediate_rData/llsat.rData")
}
load("intermediate_rData/llsat.rData")

# fit samplewise screening model
if (FALSE) {
  n <- nrow(ef) # data points
  N <- n / 6    # site visits (samples)
  llscreen <- rep(NA, N)
  sIDs <- sort(unique(ef $ sampleID))
  for (i in 1:N) {
    if (i%%10==0) cat("done", i, "of", N, "     \r"); flush.console()
    samp <- subset(ef, sampleID == sIDs[i])  
    mod <- efp(n ~ pass23 + LifeStage, data = samp, pass = pass, hessian = FALSE)
    llscreen[i] <- logLik(mod)
  }
  cat("Done!\n")
  save(llscreen, file = "intermediate_rData/llscreen.rData")
}
load("intermediate_rData/llscreen.rData")


# calculate deviance components per sample
dscreen <- 2 * (llsat - llscreen)
dscreen[dscreen < 0 & dscreen > -1e-5] <- 0
summary(dscreen)

# calculate deviance per site visit and compare to chisq 1
pchi <- pchisq(dscreen, 1)
summary(pchi)
plot(sort(dscreen))
hist(pchi, nclass = 50) # should be uniform i think?
hist(pchi[dscreen > 0], nclass = 50) # should be uniform i think?
mean(pchi > 0.95)
mean(pchi > 0.99)
# there look to be some outliers...


# these should be uniform...
hist(pchi[pchi > 1e-4], prob = TRUE, nclass = 40)
abline(h = 1, lty = 2)
# these are not very different from uniform so dont drop any
# another test might be to consider the bejamini false discovery rate
mean(1-pchi < 0.05)
mean(1-pchi < 0.01)

mean(p.adjust(1-pchi, "BH") < 0.05)
mean(p.adjust(1-pchi, "BH") < 0.01)

mean(p.adjust(1-pchi, "holm") < 0.05)
mean(p.adjust(1-pchi, "holm") < 0.01)

save(ef, file = "intermediate_rData/screenedData.rData")

