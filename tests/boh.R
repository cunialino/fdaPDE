library(mgcv)
library(fdaPDE)
rm(list = ls())
graphics.off()

source("Plots.R")
source("Utils.R")
source("Settings.R")

set <- settings(T)
set$betas <- c(-.2, .3)
set$lambdaT = 1
set$lambdaS <- 10^seq(-3, 2, .1)
set$scale = 0.1


# field <- set$f(loc[, 1], loc[, 2], t, FAMILY = "gamma")
field <- set$f(set$space_time_locations[, 2], set$space_time_locations[, 3], set$space_time_locations[, 1], "gamma")

desmat <- matrix(0, nrow = nrow(set$space_time_locations), ncol = 2)
desmat[ , 1] <- rbeta(n = nrow(desmat), shape1 = 1.5, shape2 = 2) 
desmat[ , 2] <- rbeta(n = nrow(desmat), shape1 = 3, shape2 = 2)
field = field + desmat %*% set$betas
data = rgamma(n = nrow(set$space_time_locations), shape = -1/field/set$scale,  scale = set$scale)
data <- matrix(data, nrow = nrow(set$loc), ncol = length(set$time_mesh))

sol <- smooth.FEM.time(time_mesh = set$time_mesh, locations = set$loc, observations = data, FEMbasis = set$FEMbasis, covariates = desmat, family = "gamma", 
                  DOF.evaluation = "stochastic", lambda.selection.criterion = "grid", lambda.selection.lossfunction = "GCV", 
                  max.steps.FPIRLS = 15, lambdaT = 1, lambdaS = set$lambdaS, threshold.FPIRLS = 10^-6, FLAG_PARABOLIC = T)

plotgrid <- expand.grid(set$xvec, set$yvec)
evalsol <- eval.FEM.time(sol$fit.FEM, locations = plotgrid, time.instants = 0, lambdaS = sol$bestlambda[1])

mat <- matrix(evalsol[, 1], nrow = length(set$xvec), ncol = length(set$yvec))
image(set$xvec, set$yvec, mat)
contour(set$xvec, set$yvec, mat, add = T)
trumat <- matrix(set$f(plotgrid[, 1], plotgrid[, 2], 0, "gamma"), nrow = length(set$xvec), ncol = length(set$yvec))
image(set$xvec, set$yvec, trumat)
contour(set$xvec, set$yvec, trumat , add = T)

print(sol$beta[, sol$bestlambda])
print(sqrt(mean((trumat - mat)^2, na.rm = T)))








x = NULL
x$ll = 1


test <- function(x){
  x$t <<- 5
  log(x$ll)
}
test(x)
x
