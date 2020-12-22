library(fdaPDE)
library(RColorBrewer)
library(mgcv)
rm(list = ls())
graphics.off()

source("Settings.R")
source("Utils.R")
source("Plots.R")
set <- settings(F)
# Test function
set$f = function(x, y, z, FAMILY) {
          a1 = 1
          a2 = 4
          (a1 * sin(2 * pi * x) * cos(2 * pi * y) +
           a2 * sin(3 * pi * x)) * cos(pi * z)
        }



fams <- c("poisson")
set$PDE_parameters = 
    list(K = matrix(c(.1, 0, 0, 6), nrow = 2), b = c(0, 0), c = 0)
set$dir_name =
    paste("PDE_", "P", 2 * set$basis_order, "_", "N", set$N, "_", "M",
          set$M, "_", "n", set$n, "_", "m", set$m, sep = "")
if (!dir.exists(set$dir_name)) {
  dir.create(set$dir_name)
  dir.create(paste(set$dir_name, "/figures", sep = ""))
}
set$NSIM <- 20
set$lambdaS <- 10 ^ seq(-5, 0, 0.5)
set$lambdaT <- 10 ^ seq(-7, -2, 0.5)
set$lambdaSs <- 10 ^ seq(-4, -2, 0.5)
set$lambdaTs <- 10 ^ seq(-3, -2, 0.5)

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
  boxplot(list(RMSE$GSRPDE, RMSE$GSRtPDE, RMSE$SOAP, RMSE$TPS), names = NULL)
  dev.off()

  save(list = c("RMSE", "sims", "set", "time"),
       file = paste(set$dir_name, "/", set$FAMILY, ".RData", sep = ""))

}

