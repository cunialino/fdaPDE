library(mgcv)
library(fdaPDE)
rm(list = ls())
graphics.off()

source("Plots.R")
source("Utils.R")
source("Settings.R")


plotgrid <- expand.grid(set$xvec, set$yvec)

t = 0
trumat <- matrix(set$f(plotgrid[, 1], plotgrid[, 2], t, "gamma"), nrow = length(set$xvec), ncol = length(set$yvec))
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
