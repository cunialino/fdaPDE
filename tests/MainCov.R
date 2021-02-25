library(mgcv)
library(fdaPDE)
library(RColorBrewer)

rm(list = ls())

graphics.off()

setwd("/home/elia/Università/Tesi/fdaPDE/tests")

source("Plots.R")
source("Utils.R")
source("Settings.R")

set <- settings(T)
set$betas <- c(-.2, .3)
set$scale <- 0.1

fams <- c("gamma", "poisson", "exponential")
set$time_locations <- set$time_locations * pi
set$time_mesh <- set$time_mesh * pi
set$space_time_locations[, 1] <- set$space_time_locations[, 1] * pi
set$evalGrid$t <- set$evalGrid$t * pi
set$NSIM <- 25
set$inflfac <- 1
set$maxiters <- 15
for (fam in fams) {
  set$FAMILY <- fam

  RMSE <- NULL
  betas <- NULL
  sigmas <- NULL

  time <- NULL
  time$GSRtPDE <- 0
  time$GSRPDE <- 0
  time$TPS <- 0
  time$SOAP <- 0

  print(paste(set$dir_name, set$FAMILY, set$NSIM))
  avg.coeffs <- NULL
  for (sim in 1:set$NSIM) {
    print(paste("Sim #", sim, sep = ""))
    sims <- NULL
    sims <- runsim(set, c(T, T, T, T), sim)
    R <- eval.rmse(sims, set, covariates = !is.null((set$betas)))

    time$GSRtPDE <- time$GSRtPDE + sims$timeGSRtPDE
    time$GSRPDE <- time$GSRPDE + sims$timeGSRPDE
    time$TPS <- time$TPS + sims$timeTPS
    time$SOAP <- time$SOAP + sims$timeSOAP

    RMSE$GSRtPDE <- c(RMSE$GSRtPDE, R$GSRtPDE)
    RMSE$GSRPDE <- c(RMSE$GSRPDE, R$GSRPDE)
    RMSE$TPS <- c(RMSE$TPS, R$TPS)
    RMSE$SOAP <- c(RMSE$SOAP, R$SOAP)

    if (set$FAMILY == "gamma") {
      sigmas$GSRtPDE <- c(sigmas$GSRtPDE, sims$GSRtPDE$variance.est[
        sims$GSRtPDE$bestlambda[1],
        sims$GSRtPDE$bestlambda[2]
      ])
      sigmas$GSRPDE <- c(sigmas$GSRPDE, sims$GSRPDE$variance.est[
        sims$GSRPDE$bestlambda[1],
        sims$GSRPDE$bestlambda[2]
      ])
      sigmas$TPS <- c(sigmas$TPS, sims$TPS$scale)
      sigmas$SOAP <- c(sigmas$SOAP, sims$SOAP$scale)
    }
    betas$GSRtPDE <- rbind(betas$GSRtPDE, sims$GSRtPDE$beta[
      ,
      sims$GSRtPDE$bestlambda[1],
      sims$GSRtPDE$bestlambda[2]
    ])
    betas$GSRPDE <- rbind(betas$GSRPDE, sims$GSRPDE$beta[
      ,
      sims$GSRPDE$bestlambda[1],
      sims$GSRPDE$bestlambda[2]
    ])
    betas$TPS <- rbind(betas$TPS, sims$TPS$coeff[2:3])
    betas$SOAP <- rbind(betas$SOAP, sims$SOAP$coeff[2:3])

    avg.coeffs$GSRtPDE[[sim]] <- sims$GSRtPDE$fit.FEM.time$coeff[
      ,
      sims$GSRtPDE$bestlambda[1],
      sims$GSRtPDE$bestlambda[2]
    ]
    avg.coeffs$GSRPDE[[sim]] <- sims$GSRPDE$fit.FEM.time$coeff[
      ,
      sims$GSRPDE$bestlambda[1],
      sims$GSRPDE$bestlambda[2]
    ]
    avg.coeffs$TPS[[sim]] <- coef(sims$TPS)
    avg.coeffs$SOAP[[sim]] <- coef(sims$SOAP)
  }
  avg.coeffs <- lapply(avg.coeffs, function(x) {
    Reduce("+", x) / length(x)
  })
  if (!is.null(avg.coeffs$GSRtPDE)) {
    sims$GSRtPDE$fit.FEM.time <- FEM.time(avg.coeffs$GSRtPDE, set$time_mesh,
                                          set$FEMbasis, T)
  }
  if (!is.null(avg.coeffs$GSRPDE)) {
    sims$GSRPDE$fit.FEM.time <- FEM.time(avg.coeffs$GSRPDE, set$time_mesh,
                                         set$FEMbasis, F)
  }
  if (!is.null(avg.coeffs$TPS)) {
    sims$TPS$coefficients <- avg.coeffs$TPS
  }
  if (!is.null(avg.coeffs$SOAP)) {
    sims$SOAP$coefficientsf <- avg.coeffs$SOAP
  }

  plot.results(sims, set, paste(set$dir_name, "/figures/", set$FAMILY,
    ".pdf",
    sep = ""
  ))
  pdf(file = paste(set$dir_name, "/figures/RMSE", set$FAMILY, ".pdf",
    sep =
      ""
  ))
  par(
    oma = c(1, 1, 1, 0), mar = c(1, 1, 1, 1),
    cex.axis = 0.8, font.axis = 2, font.lab = 2
  )
  boxplot(list(RMSE$TPS, RMSE$SOAP, RMSE$GSRtPDE, RMSE$GSRPDE),
    names = c("TPS", "SOAP", "GSRtPDE", "GSRPDE")
  )
  dev.off()

  if (set$FAMILY == "gamma") {
    betas$TPS <- -betas$TPS
    betas$SOAP <- -betas$SOAP
  }
  RMSEb <- function(x) {
    tmp <- x
    tmp[, 1] <- tmp[, 1] + .2
    tmp[, 2] <- tmp[, 2] - .3
    return(sqrt(colMeans(tmp^2)))
  }
  RMSEbeta <- lapply(betas, RMSEb)
  beta1 <- lapply(betas, function(x) {
    x[, 1]
  })
  beta2 <- lapply(betas, function(x) {
    x[, 2]
  })

  pdf(file = paste(set$dir_name, "/figures/BPBeta", set$FAMILY, ".pdf",
    sep
    = ""
  ))
  par(
    mfrow = c(1, 2), oma = c(1, 1, 1, 0), mar = c(1, 1, 1, 1),
    cex.axis = 0.5, font.axis = 2, font.lab = 2
  )
  boxplot(
    list(beta1$TPS, beta1$SOAP, beta1$GSRtPDE, beta1$GSRPDE),
    names = c("TPS", "SOAP", "GSRtPDE", "GSRPDE")
  )
  abline(h = -.2, col = "red")
  boxplot(
    list(beta2$TPS, beta2$SOAP, beta2$GSRtPDE, beta2$GSRPDE),
    names = c("TPS", "SOAP", "GSRtPDE", "GSRPDE")
  )
  abline(h = .3, col = "red")
  dev.off()
  save(
    list = c("RMSE", "sims", "set", "sigmas", "betas"),
    file = paste(set$dir_name, "/", set$FAMILY, ".RData", sep = "")
  )
}
