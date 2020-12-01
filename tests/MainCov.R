library(mgcv)
library(fdaPDE)
rm(list = ls())
graphics.off()

source("Plots.R")
source("Utils.R")
source("Settings.R")

set <- settings(T)
set$betas <- c(-.2, .3)
set$scale = 0.05 
fams <- c("poisson") #, "poisson")#, "poisson")
set$time_locations <- set$time_locations*pi
set$time_mesh <- set$time_mesh*pi
set$space_time_locations[, 1] <- set$space_time_locations[, 1]*pi
set$evalGrid$t <- set$evalGrid$t*pi
set$mesh <- refine.mesh.2D(set$mesh, maximum_area = 0.025, minimum_angle = 30)
# set$lambdaT <- 10^0 #seq(-5, -3, 0.5)
# set$lambdaS <- 10^seq(-3, 0, .3)
# set$lambdaTs <- set$lambdaT
set$NSIM = 1
# BC = NULL
# BC$BC_indices = set$mesh$nodesmarkers[1]
# BC$BC_values = set$f(set$mesh$nodes[BC$BC_indices, 1], set$mesh$nodes[BC$BC_indices, 2], 0, set$FAMILY)
# set$BC = BC

for(fam in fams){
  set$FAMILY = fam
  
  RMSE = NULL
  betas = NULL
  
  time = NULL
  time$GSRPDE = 0
  time$GSRtPDE = 0
  time$TPS = 0
  time$SOAP = 0
  
  print(paste(set$dir_name, set$FAMILY, set$NSIM))
  for(sim in 1:set$NSIM){
    print(paste("Sim #", sim, sep = ""))
    sims <- NULL
    sims <- runsim(set, c(T, F, F, F), sim)
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
    
    betas$GSRPDE <- rbind(betas$GSRPDE, sims$GSRPDE$beta[, sims$GSRPDE$bestlambda[1], sims$GSRPDE$bestlambda[2]])
    betas$GSRtPDE <- rbind(betas$GSRtPDE, sims$GSRtPDE$beta[, sims$GSRtPDE$bestlambda[1], sims$GSRtPDE$bestlambda[2]])
    betas$TPS <- rbind(betas$TPS, sims$TPS$coeff[2:3])
    betas$SOAP <- rbind(betas$SOAP, sims$SOAP$coeff[2:3])
  }
  

  jpeg(filename = paste(set$dir_name, "/figures/RMSE", set$FAMILY, ".jpeg", sep = ""), width = 1200, height = 1200, res = 300)
  boxplot(RMSE)
  dev.off()

  save(list = c("RMSE", "sims", "set", "time", "betas"), file = paste(set$dir_name, "/", set$FAMILY, ".RData", sep = ""))
    
}
