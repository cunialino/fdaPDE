library(mgcv)
library(fdaPDE)
rm(list = ls())
graphics.off()

source("Plots.R")
source("Utils.R")
source("Settings.R")

set <- settings(T)
fams <- c("poisson")#c("gamma", "binomial", "poisson", "exponential")
set$time_locations <- set$time_locations*pi
set$time_mesh <- set$time_mesh*pi
set$space_time_locations[, 1] <- set$space_time_locations[, 1]*pi
set$evalGrid$t <- set$evalGrid$t*pi
# set$lambdaT <- 10^seq(-5, -3, 0.5)
set$NSIM = 5

for(fam in fams){
  set$FAMILY = fam
  
  RMSE = NULL
  
  time = NULL
  time$GSRPDE = 0
  time$GSRtPDE = 0
  time$TPS = 0
  time$SOAP = 0
  
  print(paste(set$dir_name, set$FAMILY, set$NSIM))
  for(sim in 1:set$NSIM){
    print(paste("Sim #", sim, sep = ""))
    sims <- NULL
    sims <- runsim(set, c(T, T, T, T), sim == 1)
    if(sim == 1){
      plot.results(sims, set, paste(set$dir_name, "/figures/", set$FAMILY, ".jpeg", sep = ""))
    }
    R <- eval.rmse(sims, set, covariates = ! is.null((set$betas)))
    
    time$GSRPDE = time$GSRPDE + sims$timeGSRPDE
    time$GSRtPDE = time$GSRtPDE + sims$timeGSRtPDE
    time$TPS = time$TPS + sims$timeTPS
    time$SOAP = time$SOAP + sims$timeSOAP
    
    RMSE$GSRPDE  = c(RMSE$GSRPDE, R$GSRPDE)
    RMSE$GSRtPDE = c(RMSE$GSRtPDE, R$GSRtPDE)
    RMSE$TPS     = c(RMSE$TPS, R$TPS)
    RMSE$SOAP    = c(RMSE$SOAP, R$SOAP)
  }
  

  jpeg(filename = paste(set$dir_name, "/figures/RMSE", set$FAMILY, ".jpeg", sep = ""), width = 1200, height = 1200, res = 300)
  boxplot(RMSE)
  dev.off()

  save(list = c("RMSE", "sims", "set", "time"), file = paste(set$dir_name, "/", set$FAMILY, ".RData", sep = ""))
    
}
