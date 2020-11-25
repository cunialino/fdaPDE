#plot settings
{
  titlecex = .5
  labcex = .3
  pointcex = .3
  textcex = .3
  lwd = .3
}

plot.settings <- function(f, inv.link, data, set, filename){
  TrueMats = NULL
  ParamMats = NULL
  loc_mat <- cbind(rep(set$xvec, each = length(set$yvec)), rep(set$yvec, length(set$xvec)))
  for(j in seq(1, length(set$time_locations), length.out = 5)){
    TrueMats[[j]] <- matrix(set$f(loc_mat[, 1], loc_mat[, 2], set$time_locations[j], set$FAMILY), nrow = length(set$xvec), ncol = length(set$yvec), byrow = T)
    ParamMats[[j]] <- set$inv.link(TrueMats[[j]], set$FAMILY)
  }
  paramPalette = hcl.colors(100, palette = "RdYlGn")
  jpeg(filename = filename, width = 1200, height = 1200, res = 300)
  par(mfrow = c(5, 3), oma = c(.1, 2, 2, .1), mar = c(.1, .1, .1, .1))
  zl = range(unlist(TrueMats), na.rm = T)
  zlP = range(unlist(ParamMats), na.rm = T)
  groups <-  quantile(data, probs = c(0.35, 0.7))
  for(j in seq(1, length(set$time_locations), length.out = 5)){
    levsP = round(seq(min(ParamMats[[j]], na.rm = T), max(ParamMats[[j]], na.rm = T), length.out = 10), digits = 2)
    image(set$xvec,set$yvec,TrueMats[[j]],col=heat.colors(100), xlab = "", ylab = "", axes = F, cex = textcex, zlim = zl)
    contour(set$xvec,set$yvec,TrueMats[[j]],add=TRUE, labcex = labcex, lwd = lwd)
    mtext(paste("t=", round(set$time_locations[j], digits = 2), sep = ""), side = 2, line = 0, las = 2, cex = titlecex)
    if(j == 1){
      mtext("True Field", font = 2, cex = titlecex)
    }
    image(set$xvec,set$yvec,ParamMats[[j]],col=paramPalette, xlab = "", ylab = "", axes = F, cex = textcex, zlim = zlP)
    contour(set$xvec,set$yvec,ParamMats[[j]],add=TRUE, labcex = labcex, levels = levsP, lwd = lwd, zlim = zlP)
    if(j == 1){
      mtext("True Mean", font = 2, cex = titlecex)
    }
    plot(set$loc[which(data[, j] <= groups[1]), ], pch = 19, col = paramPalette[1], cex = pointcex, axes = F, xlim = range(set$xvec), ylim = range(set$yvec))
    points(set$loc[which(data[, j] > groups[1] & data[, j] < groups[2]), ], pch = 19, col = paramPalette[50], cex = pointcex)
    points(set$loc[which(data[, j] >= groups[2]), ], pch = 19, col = paramPalette[100], cex = pointcex)
    box()
    if(j == 1){
      mtext("Data", font = 2, cex = titlecex)
    }
  }
  dev.off()
}

plot.results<- function(sols, set, filename){
  
  TrueMats = NULL
  Mats = NULL
  MatsSep = NULL
  TPSMats = NULL
  SOAPMats = NULL
  c = 1
  if(! is.null(sols$TPS) | ! is.null(sols$SOAP))
    if(sum(set$FAMILY == c("gamma", "exponential")) >= 1)
      c = -1
  evalGrid <- expand.grid(c(0, .25, .5, .75, 1), set$xvec, set$yvec)
  evalGrid <- evalGrid[order(evalGrid[, 1]), ]
  if(! is.null(set$betas)){
    
    evalGrid$cov1 <- rep(0, nrow(evalGrid))
    evalGrid$cov2 <- rep(0, nrow(evalGrid))
  }
  names(evalGrid) <- c("t", "x", "y")
  if(! is.null(sols$GSRPDE))
    gsrpde <- eval.FEM.time(sols$GSRPDE$fit.FEM.time, space.time.locations = evalGrid[, 1:3], lambdaS = sols$GSRPDE$bestlambda[1], lambdaT = sols$GSRPDE$bestlambda[2])
  if(! is.null(sols$GSRtPDE))
    gsrtpde <- eval.FEM.time(sols$GSRtPDE$fit.FEM.time, space.time.locations = evalGrid[,1:3], lambdaS = sols$GSRtPDE$bestlambda[1], lambdaT = sols$GSRtPDE$bestlambda[2])
  if(! is.null(sols$TPS))
    tps <- predict.gam(sols$TPS, newdata = evalGrid) #, block.size = 10, type = "lpmatrix")%*%coef(sols$TPS)
  if(! is.null(sols$SOAP))
    soap <- predict.gam(sols$SOAP, newdata = evalGrid) #, block.size = 10, type = "lpmatrix")%*%coef(sols$SOAP)
  true <- set$f(evalGrid$x, evalGrid$y, evalGrid$t, set$FAMILY)
  nx <- length(set$xvec)
  ny <- length(set$yvec)
  for(j in 1:5){
    TrueMats[[j]]<- matrix(true[(1+(j-1)*nx*ny):(j*nx*ny)], nrow = nx, ncol = ny)
    if(! is.null(sols$GSRPDE))
       Mats[[j]] <- matrix(gsrpde[(1+(j-1)*nx*ny):(j*nx*ny)], nrow = nx, ncol = ny)
    if(! is.null(sols$GSRtPDE))
       MatsSep[[j]] <- matrix(gsrtpde[(1+(j-1)*nx*ny):(j*nx*ny)], nrow = nx, ncol = ny)
    if(! is.null(sols$TPS))
      TPSMats[[j]] <- c*matrix(tps[(1+(j-1)*nx*ny):(j*nx*ny)], nrow = nx, ncol = ny)
    if(! is.null(sols$SOAP))
      SOAPMats[[j]] <- c*matrix(soap[(1+(j-1)*nx*ny):(j*nx*ny)], nrow = nx, ncol = ny)
  }
  if(! is.null(TPSMats)){
    TPSMats = lapply(TPSMats, function(x) x+0*TrueMats[[1]])
  }
  whichones <- c(!is.null(Mats), !is.null(MatsSep), !is.null(TPSMats), !is.null(SOAPMats))
  numplot <- 1 + sum(whichones)
  
  zl =range(TrueMats, Mats, TPSMats, SOAPMats, na.rm = T)
  jpeg(filename = filename, width = 1200, height = 1200, res = 300)
  par(mfrow = c(5, numplot), oma = c(.1, 2, 2, .1), mar = c(.1, .1, .1, .1))
  time_mesh = seq(0, 1, length.out = 5)
  for(j in seq(1, 5)){
    lvR <- range(TrueMats[[j]], Mats[[j]], TPSMats[[j]], SOAPMats[[j]], na.rm = T)
    levs = round(seq(lvR[1], lvR[2], length.out = 10), digits = 1)
    image(set$xvec,set$yvec,TrueMats[[j]],col=heat.colors(100), xlab = "", ylab=paste("Time ", time_mesh[j]), zlim = zl, axes = F, cex = textcex)
    contour(set$xvec,set$yvec,TrueMats[[j]],add=TRUE, labcex = labcex, zlim = zl, levels = levs, lwd = lwd)
    mtext(paste("t=", round(time_mesh[j], digits = 2), sep = ""), side = 2, line = 0, las = 2, cex = titlecex)
    if(j == 1){
      mtext("True Field", font = 2, cex = titlecex)
    }
    if(whichones[1]){
      image(set$xvec,set$yvec,Mats[[j]],col=heat.colors(100), xlab = "", ylab = "", zlim = zl, axes = F, cex = textcex)
      contour(set$xvec,set$yvec,Mats[[j]],add=TRUE, labcex = labcex, zlim = zl, levels = levs, lwd = lwd) 
      if(j == 1){
        mtext("GSR-PDE", font = 2, cex = titlecex)
      }
    }
    if(whichones[2]){
      image(set$xvec,set$yvec,MatsSep[[j]],col=heat.colors(100), xlab = "", ylab = "", zlim = zl, axes = F, cex = textcex)
      contour(set$xvec,set$yvec,MatsSep[[j]],add=TRUE, labcex = labcex, zlim = zl, levels = levs, lwd = lwd) 
      if(j == 1){
        mtext("GSRt-PDE", font = 2, cex = titlecex)
      }
    }
    if(whichones[3]){
      image(set$xvec,set$yvec,TPSMats[[j]],col=heat.colors(100), xlab = "", ylab = "", zlim = zl, axes = F, cex = textcex)
      contour(set$xvec,set$yvec,TPSMats[[j]],add=TRUE, labcex = labcex, zlim = zl, levels = levs, lwd = lwd)
      if(j == 1){
        mtext("TPS", font = 2, cex = titlecex)
      }
    }
    if(whichones[4]){
      image(set$xvec,set$yvec,SOAPMats[[j]],col=heat.colors(100), xlab = "", ylab = "", zlim = zl, axes = F, cex = textcex)
      contour(set$xvec,set$yvec,SOAPMats[[j]],add=TRUE, labcex = labcex, zlim = zl, levels = levs, lwd = lwd)
      if(j == 1){
        mtext("SOAP", font = 2, cex = titlecex)
      }
    }
  }
  dev.off()
}
