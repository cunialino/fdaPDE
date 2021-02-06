library(mgcv)
library(fdaPDE)
library(RColorBrewer)

rm(list = ls())
graphics.off()

source("Plots.R")
source("Utils.R")
source("Settings.R")

# For square domain!
# f = function(x, y, t, FAMILY = "gamma"){
#    a = 1
#    b = 0
#    c1 = 1
#    c2 = 1
#    if(sum(FAMILY == c("gamma", "exponential"))>=1){
#        a = -1
#        b = 2
#    }
#    #return(1/a*(1/(4*pi*t)*exp(-(x^2+y^2)/(4*pi*t))+b))
#    return (1/a*(sin(2*pi*x)*sin(2*pi*y)*exp(-sqrt(8*pi^2)*t)+b))
# }
set <- settings(T)
fams <- "gamma" # c("gamma", "binomial", "poisson", "exponential")
# set$f <- f
tau <- pi
set$scale <- .1
set$time_locations <- set$time_locations * tau
set$time_mesh <- set$time_mesh * tau
set$space_time_locations[, 1] <- set$space_time_locations[, 1] * tau
set$evalGrid$t <- set$evalGrid$t * tau
#set$lambdaS  <-  10^seq(-1, 0, 1) #seq(-4, -2, 0.5)
#set$lambdaT <- 10^seq(-1, 0, 1) #seq(-5, -3, 0.5)
# set$lambdaTs <- 10^seq(-5, -3, 0.5)
set$mesh <- refine.mesh.2D(set$mesh, maximum_area = .025, minimum_angle = 30)
set$NSIM <- 5
set$maxiters  <- 100
set$inflfac <- 3
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
    sims <- runsim(set, c(T, F, F, T), sim)
    R <- eval.rmse(sims, set, covariates = !is.null((set$betas)))

    time$GSRtPDE <- time$GSRtPDE + sims$timeGSRtPDE
    time$GSRPDE <- time$GSRPDE + sims$timeGSRPDE
    time$TPS <- time$TPS + sims$timeTPS
    time$SOAP <- time$SOAP + sims$timeSOAP
    print(sims$GSRtPDE$GCV)
    print(sims$GSRtPDE$bestlambda)

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
    avg.coeffs$GSRtPDE[[sim]] <- sims$GSRtPDE$fit.FEM.time$coeff[
      ,
      sims$GSRtPDE$bestlambda[
        1
      ],
      sims$GSRtPDE$bestlambda[
        2
      ]
    ]
    avg.coeffs$GSRPDE[[sim]] <- sims$GSRPDE$fit.FEM.time$coeff[
      ,
      sims$GSRPDE$bestlambda[
        1
      ],
      sims$GSRPDE$bestlambda[
        2
      ]
    ]
    avg.coeffs$TPS[[sim]] <- coef(sims$TPS)
    avg.coeffs$SOAP[[sim]] <- coef(sims$SOAP)
  }
  avg.coeffs <- lapply(avg.coeffs, function(x) {
    Reduce("+", x) / length(x)
  })
  if (!is.null(avg.coeffs$GSRtPDE)) {
    sims$GSRtPDE$fit.FEM.time <- FEM.time(
      avg.coeffs$GSRtPDE,
      set$time_mesh, set$FEMbasis, T
    )
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

  save(
    list = c("RMSE", "sims", "set", "sigmas"),
    file = paste(set$dir_name, "/", set$FAMILY, ".RData", sep = "")
  )
}
