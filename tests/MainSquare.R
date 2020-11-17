library(mgcv)
library(fdaPDE)
rm(list = ls())

source("Plots.R")
source("Utils.R")
source("Settings.R")

set <- settings(T)

fams <- c("gamma", "binomial", "poisson", "exponential")
set$FAMILY = fams[1]

set$NSIM = 1

RMSE = NULL

RMSE$GSRPDE = 0
RMSE$GSRtPDE = 0
RMSE$TPS = 0
RMSE$SOAP = 0

print(paste(set$dir_name, set$FAMILY, NSIM))
for(sim in 1:set$NSIM){
  print(paste("Sim #", sim, sep = ""))
  sims <- runsim(set, c(T, T, F, F), sim == 1)
  R <- eval.rmse(sims, sim == 1, set, paste(set$dir_name, "/figures/", set$FAMILY, ".jpeg", sep = ""))
  RMSE$GSRPDE  = RMSE$GSRPDE  + R$GSRPDE
  RMSE$GSRtPDE = RMSE$GSRtPDE + R$GSRtPDE
  RMSE$TPS     = RMSE$TPS     + R$TPS
  RMSE$SOAP    = RMSE$SOAP    + R$SOAP
}

RMSE <- lapply(RMSE, function(x) sqrt(1/NSIM*x))

jpeg(filename = paste(set$dir_name, "/figures/RMSE", set$FAMILY, ".jpeg", sep = ""), width = 1200, height = 1200, res = 300)
boxplot(RMSE, names = c("GSR-PDE", "GSRt-PDE", "TPS", "SOAP"))
dev.off()

save(list = c("RMSE", "sims", "set"), file = paste(set$dir_name, "/", set$FAMILY, ".RData", sep = ""))
