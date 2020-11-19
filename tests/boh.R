library(fdaPDE)
library(mgcv)

load("P2_M10_n600_m41/gamma.RData")

evalGrid <- expand.grid(seq(0, 1, 0.05), seq(-1, 4, .025), seq(-1, 1, .025))
evalGrid <- evalGrid[order(evalGrid[ ,1]), ]
names(evalGrid) <- c("t", "x", "y")
time = Sys.time()
boh <- eval.FEM.time(sims$GSRPDE$fit.FEM.time, space.time.locations = evalGrid) #, bary.locations = sims$GSRPDE$bary.locations)
print(Sys.time() - time)

time = Sys.time()
boh3 <- predict.gam(sims$SOAP, newdata = evalGrid, type = "lpmatrix", block.size = 100)%*%coef(sims$SOAP)
print(Sys.time() - time)


space <- evalGrid[which(evalGrid[, 1] == 0), 2:3]
ny <- dim(space[which(space[ ,1] == -1), ])[1]
nx <- dim(space[which(space[ ,2] == -1), ])[1]

image(matrix(bohboh[1:(nx*ny)], nrow = nx, ncol = ny))


