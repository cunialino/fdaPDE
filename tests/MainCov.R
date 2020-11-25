library(mgcv)
library(fdaPDE)
rm(list = ls())
graphics.off()

source("Plots.R")
source("Utils.R")
source("Settings.R")

set <- settings(T)
set$betas <- c(-.2, .3)
fams <- c("gamma") #, "binomial", "poisson", "exponential")
set$f <- function(x, y, t, family = "gamma"){
  
  a = 8
  b = 10
  return(-1/a*(fs.test(x, y, exclude = T)+b )*(t+1))
  
}
set$NSIM = 1

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
    sims <- runsim(set, c(T, F, F, F), sim == 1)
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