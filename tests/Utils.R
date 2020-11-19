runsim <- function(set, whichone = c(T, T, T, T), plotF = F){
  
  field = set$f(set$space_time_locations[, 2], set$space_time_locations[, 3], set$space_time_locations[, 1], set$FAMILY)
  
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
  
  data = matrix(data, nrow = nrow(set$loc), ncol = length(set$time_locations))
  mode(data) <- "double"
  if(plotF)
    plot.settings(set$f, set$inv.link, data, set, paste(set$dir_name, "/figures/Settings", set$FAMILY, ".jpeg", sep = ""))
  
  reslist = NULL
  
  if(whichone[1]){
    tt <- Sys.time()
    reslist$GSRPDE <- smooth.FEM.time(time_mesh = set$time_mesh, time_locations = set$time_locations, locations = set$loc,  observations = data,
                                      FEMbasis = set$FEMbasis, covariates = NULL, GCV=T, GCVmethod = "Exact", lambdaS  = set$lambdaS, lambdaT = set$lambdaT,
                                      max.steps.FPIRLS=15,  family=set$FAMILY, mu0=NULL, scale.param=NULL, FLAG_PARABOLIC=T, threshold.FPIRLS = 10^-6)
    print(Sys.time() - tt)
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
                                       FEMbasis =set$FEMbasis, covariates = NULL, GCV=T, GCVmethod = "Exact", lambdaS  = set$lambdaSs, lambdaT = set$lambdaTs,
                                       max.steps.FPIRLS=15,  family=set$FAMILY, mu0=NULL, scale.param=NULL, FLAG_PARABOLIC=F, threshold.FPIRLS = 10^-6)
    print(Sys.time() - tt)
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
      reslist$TPS <- gam(resp ~ te(x, y, t, k=c(30,set$M),d=c(2,1),bs=c("tp","cr")), method = "GCV.Cp", family = Gamma(link = "inverse"), data = mgcvDat)
    else
      reslist$TPS <- gam(resp ~ te(x, y, t, k=c(30,set$M),d=c(2,1),bs=c("tp","cr")), method = "GCV.Cp", family = set$FAMILY, data = mgcvDat)
    print(Sys.time() - tt)
  }
  else
    reslist$TPS <- NULL
  
  print("Done")
  if(whichone[4]){
    nmax = 100
    tt <- Sys.time()
    if(set$FAMILY == "gamma" || set$FAMILY == "exponential")
      reslist$SOAP = gam(resp ~ te(x,y,t,bs=c("sf","cr"),k=c(25,4),d=c(2,1), xt=list(bnd=set$fsb))+
                      te(x,y,t,bs=c("sw","cr"),k=c(25,4),d=c(2,1), xt=list(bnd=set$fsb)),
                    knots=set$knots, family = Gamma(link = "inverse"), data = mgcvDat)
    else
      reslist$SOAP = gam(resp ~ te(x,y,t,bs=c("sf","cr"),k=c(25,4),d=c(2,1), xt=list(bnd=set$fsb))+
                      te(x,y,t,bs=c("sw","cr"),k=c(25,4),d=c(2,1), xt=list(bnd=set$fsb)),
                    knots=set$knots, family = set$FAMILY, data = mgcvDat)
    print(Sys.time() - tt)
  }
  else
    reslist$SOAP <- NULL
  print("Done")
  return(reslist)
}

eval.rmse <- function(sols, plotF = F, set, filename){
  TrueMats = NULL
  Mats = NULL
  SepMats = NULL
  TPSMats = NULL
  SOAPMats = NULL
  fine_time = seq(0, 1, .05)
  c = 1
  if(! is.null(sols$TPS) | ! is.null(sols$SOAP))
    if(sum(set$FAMILY == c("gamma", "exponential")) >= 1)
      c = -1
  evalGrid <- expand.grid(fine_time, set$xvec, set$yvec)
  evalGrid <- evalGrid[order(evalGrid[, 1]), ]
  names(evalGrid) <- c("t", "x", "y")
  true <- set$f(evalGrid$x, evalGrid$y, evalGrid$t, set$FAMILY)
  
  tt <- Sys.time()
  if(! is.null(sols$GSRPDE))
    gsrpde <- eval.FEM.time(sols$GSRPDE$fit.FEM.time, space.time.locations = evalGrid, lambdaS = sols$GSRPDE$bestlambda[1], lambdaT = sols$GSRPDE$bestlambda[2])
  print(Sys.time() - tt)
  tt <- Sys.time()
  if(! is.null(sols$GSRtPDE))
    gsrtpde <- eval.FEM.time(sols$GSRtPDE$fit.FEM.time, space.time.locations = evalGrid, lambdaS = sols$GSRtPDE$bestlambda[1], lambdaT = sols$GSRtPDE$bestlambda[2])
  print(Sys.time() - tt)
  tt <- Sys.time()
  if(! is.null(sols$TPS))
    tps <- predict.gam(sols$TPS, newdata = evalGrid, block.size = 100, type = "lpmatrix")%*%coef(sols$TPS)
  print(Sys.time() - tt)
  tt <- Sys.time()
  if(! is.null(sols$SOAP))
    soap <- predict.gam(sols$SOAP, newdata = evalGrid, block.size = 100, type = "lpmatrix")%*%coef(sols$SOAP)
  print(Sys.time() - tt)
  nx <- length(set$xvec)
  ny <- length(set$yvec)
  for(j in 1:length(fine_time)){
    TrueMats[[j]]<- matrix(true[(1+(j-1)*nx*ny):(j*nx*ny)], nrow = nx, ncol = ny)
    if(! is.null(sols$GSRPDE))
       Mats[[j]] <- matrix(gsrpde[(1+(j-1)*nx*ny):(j*nx*ny)], nrow = nx, ncol = ny)
    if(! is.null(sols$GSRtPDE))
       SepMats[[j]] <- matrix(gsrtpde[(1+(j-1)*nx*ny):(j*nx*ny)], nrow = nx, ncol = ny)
    if(! is.null(sols$TPS))
      TPSMats[[j]] <- c*matrix(tps[(1+(j-1)*nx*ny):(j*nx*ny)], nrow = nx, ncol = ny)
    if(! is.null(sols$SOAP))
      SOAPMats[[j]] <- c*matrix(soap[(1+(j-1)*nx*ny):(j*nx*ny)], nrow = nx, ncol = ny)
    # evalMGCV = data.frame("x"=rep(set$xvec, length(set$yvec)), "y"=rep(set$yvec, each = length(set$xvec)), "t"=rep(fine_time[j], length(set$xvec)*length(set$yvec)))
    # TrueMats[[j]] <- matrix(set$f(evalMGCV$x, evalMGCV$y, evalMGCV$t, set$FAMILY), nrow = length(set$xvec), ncol = length(set$yvec))
    # tt <- Sys.time()
    # if(! is.null(sols$GSRPDE)){
    #   Mats[[j]] <- eval.FEM.time(sols$GSRPDE$fit.FEM.time, locations = evalMGCV[, 1:2], time.instants = fine_time[j], lambdaS = sols$GSRPDE$bestlambda[1], lambdaT = sols$GSRPDE$bestlambda[2])
    #   Mats[[j]] <- matrix(Mats[[j]], nrow = length(set$xvec), ncol = length(set$yvec))
    #   
    # }
    # time_fem = time_fem + (Sys.time() - tt)
    # if(! is.null(sols$GSRtPDE)){
    #   SepMats[[j]] <- eval.FEM.time(sols$GSRtPDE$fit.FEM.time, locations = evalMGCV[, 1:2], time.instants = fine_time[j], lambdaS = sols$GSRtPDE$bestlambda[1], lambdaT = sols$GSRtPDE$bestlambda[2])
    #   SepMats[[j]] <- matrix(SepMats[[j]], nrow = length(set$xvec), ncol = length(set$yvec))
    # }
    # tt <- Sys.time()
    # if(! is.null(sols$TPS))
    #   TPSMats[[j]] <- c*matrix(predict.gam(sols$TPS, evalMGCV, type = "lpmatrix")%*%coef(sols$TPS), nrow = length(set$xvec), ncol = length(set$yvec))
    # time_tps = time_tps + (Sys.time() - tt)
    # tt <- Sys.time()
    # if(! is.null(sols$SOAP))
    #   SOAPMats[[j]] <- c*matrix(predict(sols$SOAP, evalMGCV, type = "lpmatrix")%*%coef(sols$SOAP), nrow = length(set$xvec), ncol = length(set$yvec))
    # time_soap = time_soap + (Sys.time() - tt)
  }
  RMSE = NULL
  # RMSE$GSRPDE  = (unlist(TrueMats) - unlist(Mats))^2
  # RMSE$GSRtPDE = (unlist(TrueMats) - unlist(SepMats))^2
  # RMSE$TPS     = (unlist(TrueMats) - unlist(TPSMats))^2
  # RMSE$SOAP    = (unlist(TrueMats) - unlist(SOAPMats))^2
  if( ! is.null(sols$GSRPDE))
    RMSE$GSRPDE = (true - gsrpde)^2
  if( ! is.null(sols$GSRtPDE))
    RMSE$GSRtPDE = (true - gsrtpde)^2
  if( ! is.null(sols$TPS))
    RMSE$TPS = (true - c*tps)^2
  if( ! is.null(sols$SOAP))
    RMSE$SOAP = (true - c*soap)^2
  if(plotF)
    plot.results(TrueMats, Mats, SepMats, TPSMats, SOAPMats, fine_time, set, filename)
  return(RMSE)
}
