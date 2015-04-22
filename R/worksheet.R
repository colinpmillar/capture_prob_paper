
library(CLmodel)


wk <- subset(ef, Species == "Salmon" & LifeStage == "Fry" & CATCH_ID == 1)

formula <- Z ~ s(year)

g <- gam(formula, data = wk, fit = FALSE)



