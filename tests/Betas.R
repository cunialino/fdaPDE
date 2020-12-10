load("P2_M11_n405_m11/poisson.RData")
  if(set$FAMILY == "gamma"){
  betas$TPS = -betas$TPS
  betas$SOAP = -betas$SOAP
}
RMSE <- function(x){
  tmp = x
  tmp[, 1] = tmp[, 1] +.2
  tmp[, 2] = tmp[, 2] - .3
  return(sqrt(colMeans(tmp^2)))
}
RMSEbeta <- lapply(betas, RMSE)
beta1 <- lapply(betas, function(x){ x[, 1]})
beta2 <- lapply(betas, function(x){ x[, 2]})

beta1 <- lapply(betas, function(x){ x[, 1]})
beta2 <- lapply(betas, function(x){ x[, 2]})
jpeg(filename = paste(set$dir_name, "/figures/BPBeta", set$FAMILY, ".jpeg", sep = ""), width = 1200, height = 1200, res = 300)
par(mfrow = c(1, 2))
boxplot(beta1, names = NULL)
abline(h = -.2, col = "red")
boxplot(beta2, names = NULL)
abline(h = .3, col = "red")
dev.off()
