runsim <- function(set, whichone = c(T, T, T, T), sim){
  
  plotF <- (sim == 1)
  desmat = NULL
  field = set$f(set$space_time_locations[, 2], set$space_time_locations[, 3], set$space_time_locations[, 1], set$FAMILY)
  if(! is.null(set$betas)){
    desmat <- matrix(0, nrow = nrow(set$space_time_locations), ncol = 2)
    desmat[ , 1] <- rbeta(n = nrow(desmat), shape1 = 1.5, shape2 = 2) 
    desmat[ , 2] <- rbeta(n = nrow(desmat), shape1 = 3, shape2 = 2)+1
    field = field + desmat %*% set$betas
    TPSform <- resp ~ cov1 + cov2 + te(x, y, t, k=c(30,5),d=c(2,1),bs=c("tp","cr"))
    SOAPform <- resp ~ cov1 + cov2 +te(x,y,t,bs=c("sf","cr"),k=c(25,4),d=c(2,1), xt=list(bnd=set$fsb))+
                      te(x,y,t,bs=c("sw","cr"),k=c(25,4),d=c(2,1), xt=list(bnd=set$fsb))
  }
  else{
    TPSform <- resp ~ te(x, y, t, k=c(30,5),d=c(2,1),bs=c("tp","cr"))
    SOAPform <- resp ~ te(x,y,t,bs=c("sf","cr"),k=c(25,4),d=c(2,1), xt=list(bnd=set$fsb))+
                      te(x,y,t,bs=c("sw","cr"),k=c(25,4),d=c(2,1), xt=list(bnd=set$fsb))
  }
  
  if(set$FAMILY == "binomial"){
    data = rbinom(n = nrow(set$space_time_locations),size = 1,  prob = plogis(field))
  }
  if(set$FAMILY == "gamma"){
    data = rgamma(n = nrow(set$space_time_locations), shape = -1/field/set$scale,  scale = set$scale)
  }
  if(set$FAMILY == "gaussian"){
    data <- field + rnorm(nrow(set$space_time_locations), mean = 0, sd = set$scale) 
  }
  if(set$FAMILY == "poisson"){
    data = rpois(n = nrow(set$space_time_locations), lambda = exp(field))
  }
  if(set$FAMILY == "exponential"){
    data = rexp(n = nrow(set$space_time_locations), rate = -field)
  }
  
  mgcvDat = data.frame("resp"=data, "t"=set$space_time_locations[, 1], "x"=set$space_time_locations[, 2], "y"=set$space_time_locations[, 3])
  if(! is.null(set$betas)){
    mgcvDat$cov1 <- desmat[, 1]
    mgcvDat$cov2 <- desmat[, 2]
  }
  
  data = matrix(data, nrow = set$n, ncol = set$m)
  if(sim == 1){
    set$data <<- data
  }
  set$desmat[[sim]] <<- desmat
  mode(data) <- "double"
  if(plotF)
    plot.settings(set$f, set$inv.link, data, set, paste(set$dir_name, "/figures/Settings", set$FAMILY, ".jpeg", sep = ""))
  
  reslist = NULL
  
  if(whichone[1]){
    tt <- Sys.time()
    reslist$GSRPDE <- smooth.FEM.time(time_mesh = set$time_mesh, time_locations = set$time_locations, locations = set$loc,  observations = data, BC = set$BC, 
                                      FEMbasis = set$FEMbasis, covariates = desmat, DOF.evaluation = "stochastic",lambda.selection.lossfunction = "GCV", lambdaS  = set$lambdaS, lambdaT = set$lambdaT,
                                      max.steps.FPIRLS=15,  family=set$FAMILY, mu0=NULL, scale.param=NULL, FLAG_PARABOLIC=T, threshold.FPIRLS = 10^-6, PDE_parameters = set$PDE_parameters)
    reslist$timeGSRPDE = Sys.time() - tt
    print(reslist$timeGSRPDE)
    
    bl <- reslist$GSRPDE$bestlambda
    if(length(set$lambdaS) > 1 && (bl[1] == 1 || bl[1] == length(set$lambdaS)))
       print(paste("Bad lambdaS:", bl[1]))
    if(length(set$lambdaT) > 1 && (bl[2] == 1 || bl[2] == length(set$lambdaT)))
      print(paste("Bad lambdaT:", bl[2]))
       
  }
  else
    reslist$GSRPDE = NULL
  print("Done")
  if(whichone[2]){
    tt <- Sys.time()
    reslist$GSRtPDE <- smooth.FEM.time(time_mesh = set$time_mesh, time_locations = set$time_locations, locations = set$loc,  observations = data,
                                       FEMbasis =set$FEMbasis, covariates = desmat, DOF.evaluation = "stochastic", lambda.selection.lossfunction = "GCV", lambdaS  = set$lambdaSs, lambdaT = set$lambdaTs,
                                       max.steps.FPIRLS=15,  family=set$FAMILY, mu0=NULL, scale.param=NULL, FLAG_PARABOLIC=F, threshold.FPIRLS = 10^-6, PDE_parameters = set$PDE_parameters, BC = set$BC)
    reslist$timeGSRtPDE = Sys.time() - tt
    print(reslist$timeGSRtPDE)
    bl <- reslist$GSRtPDE$bestlambda
    if(length(set$lambdaSs) > 1 && (bl[1] == 1 || bl[1] == length(set$lambdaSs)))
       print(paste("Bad lambdaS:", bl[1]))
    if(length(set$lambdaTs) > 1 && (bl[2] == 1 || bl[2] == length(set$lambdaTs)))
      print(paste("Bad lambdaT:", bl[2]))
  }
  else
    reslist$GSRtPDE = NULL
  print("Done")
  if(whichone[3]){
    tt <- Sys.time()
    if(set$FAMILY == "gamma" || set$FAMILY == "exponential")
      reslist$TPS <- gam(TPSform, method = "GCV.Cp", family = Gamma(link = "inverse"), data = mgcvDat)
    else
      reslist$TPS <- gam(TPSform, method = "GCV.Cp", family = set$FAMILY, data = mgcvDat)
    reslist$timeTPS = Sys.time() - tt
    print(reslist$timeTPS)
  }
  else
    reslist$TPS <- NULL
  
  print("Done")
  if(whichone[4]){
    nmax = 100
    tt <- Sys.time()
    if(set$FAMILY == "gamma" || set$FAMILY == "exponential")
      reslist$SOAP = gam(SOAPform, knots=set$knots, family = Gamma(link = "inverse"), data = mgcvDat)
    else
      reslist$SOAP = gam(SOAPform, knots=set$knots, family = set$FAMILY, data = mgcvDat)
    reslist$timeSOAP = Sys.time() - tt
    print(reslist$timeSOAP)
  }
  else
    reslist$SOAP <- NULL
  print("Done")
  return(reslist)
}

eval.rmse <- function(sols, set, covariates = F){
  
  if(covariates){
    set$evalGrid$cov1 = rep(0, nrow(set$evalGrid))
    set$evalGrid$cov2 = rep(0, nrow(set$evalGrid))
  }
  c = 1
  if(! is.null(sols$TPS) | ! is.null(sols$SOAP))
    if(sum(set$FAMILY == c("gamma", "exponential")) >= 1)
      c = -1
    
  true <- set$f(set$evalGrid$x, set$evalGrid$y, set$evalGrid$t, set$FAMILY)
  
  tt <- Sys.time()
  if(! is.null(sols$GSRPDE))
    gsrpde <- eval.FEM.time(sols$GSRPDE$fit.FEM.time, space.time.locations = set$evalGrid[, 1:3], lambdaS = sols$GSRPDE$bestlambda[1], lambdaT = sols$GSRPDE$bestlambda[2])
  print(Sys.time()-tt)
  if(! is.null(sols$GSRtPDE))
    gsrtpde <- eval.FEM.time(sols$GSRtPDE$fit.FEM.time, space.time.locations = set$evalGrid[,1:3], lambdaS = sols$GSRtPDE$bestlambda[1], lambdaT = sols$GSRtPDE$bestlambda[2])
  
  if(! is.null(sols$TPS))
    tps <- c*predict.gam(sols$TPS, newdata = set$evalGrid) #, block.size = 10, type = "lpmatrix")%*%coef(sols$TPS)
  
  if(! is.null(sols$SOAP))
    soap <- c*predict.gam(sols$SOAP, newdata = set$evalGrid) #, block.size = 10, type = "lpmatrix")%*%coef(sols$SOAP)
  
  RMSE = NULL
  
  if( ! is.null(sols$GSRPDE))
    RMSE$GSRPDE = sqrt(mean((true - gsrpde)^2, na.rm = T))
  if( ! is.null(sols$GSRtPDE))
    RMSE$GSRtPDE = sqrt(mean((true - gsrtpde)^2, na.rm = T))
  if( ! is.null(sols$TPS))
    RMSE$TPS = sqrt(mean((true - tps)^2, na.rm = T))
  if( ! is.null(sols$SOAP))
    RMSE$SOAP = sqrt(mean((true - soap)^2, na.rm = T))
  return(RMSE)
}
