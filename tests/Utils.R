inv.link <- function(mu, FAMILY){
  if(FAMILY == "gamma")
    return(-1/mu)
  if(FAMILY == "binomial")
    return(plogis(mu))
}

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
  
  mgcvDat = data.frame("resp"=data, "t"=set$space_time_locations[, 1], "x"=set$space_time_locations[, 2], "y"=set$space_time_locations[, 3])
  
  data = matrix(data, nrow = nrow(set$loc), ncol = length(set$time_locations))
  mode(data) <- "double"
  if(plotF)
    plot.settings(set$f, set$inv.link, data, set, paste(set$dir_name, "/figures/Settings", set$FAMILY, ".jpeg", sep = ""))
  
  reslist = NULL
  
  if(whichone[1]){
    reslist$GSRPDE <- smooth.FEM.time(time_mesh = set$time_mesh, time_locations = set$time_locations, locations = set$loc,  observations = data,
                                      FEMbasis = set$FEMbasis, covariates = NULL, GCV=T, GCVmethod = "Exact", lambdaS  = set$lambdaS, lambdaT = set$lambdaT,
                                      max.steps.FPIRLS=15,  family=set$FAMILY, mu0=NULL, scale.param=NULL, FLAG_PARABOLIC=T, threshold.FPIRLS = 10^-6)
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
    reslist$GSRtPDE <- smooth.FEM.time(time_mesh = set$time_mesh, time_locations = set$time_locations, locations = set$loc,  observations = data,
                                       FEMbasis =set$FEMbasis, covariates = NULL, GCV=T, GCVmethod = "Exact", lambdaS  = set$lambdaSs, lambdaT = set$lambdaTs,
                                       max.steps.FPIRLS=15,  family=set$FAMILY, mu0=NULL, scale.param=NULL, FLAG_PARABOLIC=F, threshold.FPIRLS = 10^-6)
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
    if(set$FAMILY == "gamma" || set$FAMILY == "exopnential")
      reslist$TPS <- gam(resp ~ te(x, y, t, k=c(30,10),d=c(2,1),bs=c("tp","cr")), method = "GCV.Cp", family = Gamma(link = "inverse"), data = mgcvDat)
    else
      reslist$TPS <- gam(resp ~ te(x, y, t, k=c(30,10),d=c(2,1),bs=c("tp","cr")), method = "GCV.Cp", family = set$FAMILY, data = mgcvDat)
  }
  else
    reslist$TPS <- NULL
  
  print("Done")
  if(whichone[4]){
    nmax = 100
    if(set$FAMILY == "gamma" || set$FAMILY == "exopnential")
      reslist$SOAP = gam(resp ~ te(x,y,t,bs=c("sf","cr"),k=c(25,4),d=c(2,1), xt=list(list(bnd=set$fsb,nmax=set$nmax),NULL))+
                      te(x,y,t,bs=c("sw","cr"),k=c(25,4),d=c(2,1), xt=list(list(bnd=set$fsb,nmax=set$nmax),NULL)),
                    knots=set$knots, family = Gamma(link = "inverse"), data = mgcvDat)
    else
      reslist$SOAP = gam(resp ~ te(x,y,t,bs=c("sf","cr"),k=c(25,4),d=c(2,1), xt=list(list(bnd=set$fsb,nmax=set$nmax),NULL))+
                      te(x,y,t,bs=c("sw","cr"),k=c(25,4),d=c(2,1), xt=list(list(bnd=set$fsb,nmax=set$nmax),NULL)),
                    knots=set$knots, family = set$FAMILY, data = mgcvDat)
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
  fine_time = seq(0, 1, length.out = 40)# min(2*set$M, 40))
  c = 1
  if(! is.null(sols$TPS) | ! is.null(sols$SOAP))
    if(set$FAMILY == "gamma")
      c = -1
  for(j in 1:length(fine_time)){
    evalMGCV = data.frame("x"=rep(set$xvec, length(set$yvec)), "y"=rep(set$yvec, each = length(set$xvec)), "t"=rep(fine_time[j], length(set$xvec)*length(set$yvec)))
    TrueMats[[j]] <- matrix(set$f(evalMGCV$x, evalMGCV$y, evalMGCV$t, set$FAMILY), nrow = length(set$xvec), ncol = length(set$yvec))
    if(! is.null(sols$GSRPDE)){
      Mats[[j]] <- eval.FEM.time(sols$GSRPDE$fit.FEM.time, locations = evalMGCV[, 1:2], time.instants = fine_time[j], lambdaS = sols$GSRPDE$bestlambda[1], lambdaT = sols$GSRPDE$bestlambda[2])
      Mats[[j]] <- matrix(Mats[[j]], nrow = length(set$xvec), ncol = length(set$yvec))
      
    }
    if(! is.null(sols$GSRtPDE)){
      SepMats[[j]] <- eval.FEM.time(sols$GSRtPDE$fit.FEM.time, locations = evalMGCV[, 1:2], time.instants = fine_time[j], lambdaS = sols$GSRtPDE$bestlambda[1], lambdaT = sols$GSRtPDE$bestlambda[2])
      SepMats[[j]] <- matrix(SepMats[[j]], nrow = length(set$xvec), ncol = length(set$yvec))
    }
    if(! is.null(sols$TPS))
      TPSMats[[j]] <- c*matrix(predict(sols$TPS, evalMGCV), nrow = length(set$xvec), ncol = length(set$yvec))
    if(! is.null(sols$SOAP))
      SOAPMats[[j]] <- c*matrix(predict(sols$SOAP, evalMGCV), nrow = length(set$xvec), ncol = length(set$yvec))
  }
  RMSE = NULL
  RMSE$GSRPDE  = (unlist(TrueMats) - unlist(Mats))^2
  RMSE$GSRtPDE = (unlist(TrueMats) - unlist(SepMats))^2
  RMSE$TPS     = (unlist(TrueMats) - unlist(TPSMats))^2
  RMSE$SOAP    = (unlist(TrueMats) - unlist(SOAPMats))^2
  if(plotF)
    plot.results(TrueMats, Mats, SepMats, TPSMats, SOAPMats, fine_time, set, filename)
  return(RMSE)
}
