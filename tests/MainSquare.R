library(mgcv)
library(fdaPDE)
library(RColorBrewer)

rm(list = ls())
graphics.off()

source("Plots.R")
source("Utils.R")
source("Settings.R")

#f = function(x, y, z, FAMILY = "gamma") {
#      a = -1
#      b = 0
#      if (sum(FAMILY == c("gamma", "exponential")) >= 1) {
#        a = .5
#        b = 2
#      }
#      coe = function(x, y) 1 / 2 * sin(5 * pi * x) * exp(-x ^ 2) + 1
#      -1/a*(sin(2 * pi * (coe(y, 1) * x * cos(z - 2) - y * sin(z - 2))) *
#          cos(2 * pi * (coe(y, 1) * x * cos(z - 2 + pi / 2) +
#                        coe(x, 1) * y * sin((z - 2) * pi / 2))) + b)
#    }
f = function(x, y, t, FAMILY = "gamma"){
    a = 1
    b = 0
    c1 = 1
    c2 = 1
    if(sum(FAMILY == c("gamma", "exponential"))>=1){
        a = -1
        b = 2
    }
    #return(1/a*(1/(4*pi*t)*exp(-(x^2+y^2)/(4*pi*t))+b))
    return (1/a*(sin(2*pi*x)*sin(2*pi*y)*exp(-sqrt(8*pi^2)*t)+b))
}
set <- settings(F)
fams <- c("binomial", "poisson", "exponential")
set$f <- f
tau <- 5/sqrt(8*pi^2)
set$time_locations <- set$time_locations*tau
set$time_mesh <- set$time_mesh*tau 
set$space_time_locations[, 1] <- set$space_time_locations[, 1]*tau
set$evalGrid$t <- set$evalGrid$t*tau
# set$lambdaT <- 10^seq(-5, -3, 0.5)
set$NSIM = 50

for (fam in fams) {
  set$FAMILY = fam

  RMSE = NULL

  time = NULL
  time$GSRPDE = 0
  time$GSRtPDE = 0
  time$TPS = 0
  time$SOAP = 0

  avg.coeffs = NULL
  print(paste(set$dir_name, set$FAMILY, set$NSIM))
  for (sim in 1:set$NSIM) {
    print(paste("Sim #", sim, sep = ""))
    sims <- NULL
    sims <- runsim(set, c(T, T, T, T), sim)
    R <- eval.rmse(sims, set, covariates = !is.null((set$betas)))

    time$GSRPDE = time$GSRPDE + sims$timeGSRPDE
    time$GSRtPDE = time$GSRtPDE + sims$timeGSRtPDE
    time$TPS = time$TPS + sims$timeTPS
    time$SOAP = time$SOAP + sims$timeSOAP

    RMSE$GSRPDE = c(RMSE$GSRPDE, R$GSRPDE)
    RMSE$GSRtPDE = c(RMSE$GSRtPDE, R$GSRtPDE)
    RMSE$TPS = c(RMSE$TPS, R$TPS)
    RMSE$SOAP = c(RMSE$SOAP, R$SOAP)
    avg.coeffs$GSRtPDE[[
                   sim]] <- sims$GSRtPDE$fit.FEM.time$coeff[
                                                          ,
                                                          sims$GSRtPDE$bestlambda[
                                                                           1],
                                                          sims$GSRtPDE$bestlambda[
                                                                           2]]
    avg.coeffs$GSRPDE[[
                   sim]] <- sims$GSRPDE$fit.FEM.time$coeff[
                                                         ,
                                                         sims$GSRPDE$bestlambda[
                                                                         1],
                                                         sims$GSRPDE$bestlambda[
                                                                         2]]
    avg.coeffs$TPS[[sim]] <- coef(sims$TPS)
    avg.coeffs$SOAP[[sim]] <- coef(sims$SOAP)
  }
  avg.coeffs <- lapply(avg.coeffs, function(x) {
                                     Reduce("+", x) / length(x)
                                   })
  if (!is.null(avg.coeffs$GSRtPDE)) {
    sims$GSRtPDE$fit.FEM.time <- FEM.time(
        avg.coeffs$GSRtPDE, set$time_mesh, set$FEMbasis, T)
  }
  if (!is.null(avg.coeffs$GSRPDE)) {
    sims$GSRPDE$fit.FEM.time <- FEM.time(
        avg.coeffs$GSRPDE, set$time_mesh, set$FEMbasis, F)
  }
  if (!is.null(avg.coeffs$TPS)) {
    sims$TPS$coefficients <- avg.coeffs$TPS
  }
  if (!is.null(avg.coeffs$SOAP)) {
    sims$SOAP$coefficientsf <- avg.coeffs$SOAP
  }

  plot.results(sims, set,
               paste(set$dir_name, "/figures/", set$FAMILY, ".jpeg", sep = ""))

  jpeg(filename = paste(
           set$dir_name, "/figures/RMSE", set$FAMILY, ".jpeg", sep = ""),
       width = 1200, height = 1200, res = 300)
  boxplot(list(RMSE$GSRtPDE, RMSE$GSRPDE, RMSE$SOAP, RMSE$TPS), names = NULL)
  dev.off()

  save(list = c("RMSE", "sims", "set", "time"),
       file = paste(set$dir_name, "/", set$FAMILY, ".RData", sep = ""))

}

