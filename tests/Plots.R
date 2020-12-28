#plot settings
{
    #  titlecex = .5
    #  labcex = .3
    #  pointcex = .3
    #  textcex = .3
    #  lwd = .3
    fieldPal = colorRampPalette(brewer.pal(9, "Blues")[2:7])(100)
    paramPal = colorRampPalette(brewer.pal(9, "RdYlGn"))(100)
}

plot.settings <- function(data, set, file) {
    TrueMats = NULL
    ParamMats = NULL
    times <- sort(
                  set$time_locations[seq(1, length(set$time_locations), length.out = 5)])
    data <- data[, seq(1, length(set$time_locations), length.out = 5)]
    evalGrid <- expand.grid(times, set$xvec, set$yvec)
    evalGrid <- evalGrid[order(evalGrid[, 1]), ]
    names(evalGrid) <- c("t", "x", "y")
    true <- set$f(evalGrid$x, evalGrid$y, evalGrid$t, set$FAMILY)
    nx <- length(set$xvec)
    ny <- length(set$yvec)
    for (j in 1:5) {
        TrueMats[[j]] <- matrix(
                                true[(1 + (j - 1) * nx * ny):(j * nx * ny)], nrow = nx, ncol = ny)
        ParamMats[[j]] <- set$inv.link(TrueMats[[j]], set$FAMILY)
    }
    pdf(file = file) #, width = 1200, height = 1200, res = 300)
    par(mfrow = c(5, 3), oma = c(.1, 3, 2, .1), mar = c(.1, .1, .1, .1))
    zl = range(unlist(TrueMats), na.rm = T)
    zlP = range(unlist(ParamMats), na.rm = T)
    if(set$FAMILY != "binomial")
        groups <- quantile(data, probs = c(0.35, 0.7))
    else
        groups <-c(0, 0.9)
    for (j in 1:5) {
        image(set$xvec, set$yvec, TrueMats[[j]], col = fieldPal, xlab = "",
              ylab = "", axes = F, zlim = zl)
        contour(set$xvec, set$yvec, TrueMats[[j]], add = TRUE, zlim = zl)
        mtext(paste("t=", round(times[j], digits = 2), sep = ""), side = 2,
              line = 0, las = 3)
        if (j == 1) {
            mtext("True Field", font = 2)
        }
        image(set$xvec, set$yvec, ParamMats[[j]], col = paramPal, xlab = "",
              ylab = "", axes = F, zlim = zlP)
        contour(set$xvec, set$yvec, ParamMats[[j]], add = TRUE,zlim = zlP)
        if (j == 1) {
            mtext("True Mean", font = 2)
        }
        plot(set$loc[which(data[, j] <= groups[1]), ], pch = 19, col = paramPal[1],
             axes = F, xlim = range(set$xvec),
             ylim = range(set$yvec))
        points(set$loc[which(data[, j] > groups[1] & data[, j] <= groups[2]), ],
               pch = 19, col = paramPal[40])
        points(set$loc[which(data[, j] > groups[2]), ], pch = 19,
               col = paramPal[100])
        lines(set$fsb[[1]], col = "black")
        if (j == 1) {
            mtext("Data", font = 2)
        }
    }
    dev.off()
}

plot.results <- function(sols, set, file) {

    TrueMats = NULL
    Mats = NULL
    MatsSep = NULL
    TPSMats = NULL
    SOAPMats = NULL
    c = 1
    if (!is.null(sols$TPS) | !is.null(sols$SOAP))
        if (sum(set$FAMILY == c("gamma", "exponential")) >= 1) c = -1
    times <- set$time_mesh[seq(1, length(set$time_mesh), length.out = 5)]
    evalGrid <- expand.grid(times, set$xvec, set$yvec)
    evalGrid <- evalGrid[order(evalGrid[, 1]), ]
    names(evalGrid) <- c("t", "x", "y")
    if (!is.null(set$betas)) {

        evalGrid$cov1 <- rep(0, nrow(evalGrid))
        evalGrid$cov2 <- rep(0, nrow(evalGrid))
    }
    if (!is.null(sols$GSRtPDE))
        gsrpde <- eval.FEM.time(
                                sols$GSRtPDE$fit.FEM.time, space.time.locations = evalGrid[, 1:3],
                                lambdaS = sols$GSRtPDE$bestlambda[1],
                                lambdaT = sols$GSRtPDE$bestlambda[2])
    if (!is.null(sols$GSRPDE))
        gsrtpde <- eval.FEM.time(
                                 sols$GSRPDE$fit.FEM.time, space.time.locations = evalGrid[, 1:3],
                                 lambdaS = sols$GSRPDE$bestlambda[1],
                                 lambdaT = sols$GSRPDE$bestlambda[2])
    if (!is.null(sols$TPS))
        tps <- predict.gam(
                           sols$TPS, newdata = evalGrid)  # , block.size = 10, type =
    # "lpmatrix")%*%coef(sols$TPS)
    if (!is.null(sols$SOAP))
        soap <- predict.gam(
                            sols$SOAP, newdata = evalGrid)  # , block.size = 10, type =
    # "lpmatrix")%*%coef(sols$SOAP)
    true <- set$f(evalGrid$x, evalGrid$y, evalGrid$t, set$FAMILY)
nx <- length(set$xvec)
ny <- length(set$yvec)
for (j in 1:5) {
    TrueMats[[j]] <- matrix(
                            true[(1 + (j - 1) * nx * ny):(j * nx * ny)], nrow = nx, ncol = ny)
    if (!is.null(sols$GSRtPDE))
        Mats[[j]] <- matrix(
                            gsrpde[(1 + (j - 1) * nx * ny):(j * nx * ny)], nrow = nx, ncol = ny)
    if (!is.null(sols$GSRPDE))
        MatsSep[[j]] <- matrix(
                               gsrtpde[(1 + (j - 1) * nx * ny):(j * nx * ny)], nrow = nx, ncol = ny)
    if (!is.null(sols$TPS))
        TPSMats[[j]] <- c * matrix(tps[(1 + (j - 1) * nx * ny):(j * nx * ny)],
                                   nrow = nx, ncol = ny)
    if (!is.null(sols$SOAP))
        SOAPMats[[j]] <- c * matrix(soap[(1 + (j - 1) * nx * ny):(j * nx * ny)],
                                    nrow = nx, ncol = ny)
}
if (!is.null(TPSMats)) {
    TPSMats = lapply(TPSMats, function(x) x + 0 * TrueMats[[1]])
}
whichones <- c(
               !is.null(Mats), !is.null(MatsSep), !is.null(TPSMats), !is.null(SOAPMats))
numplot <- 1 + sum(whichones)

zl = range(TrueMats, Mats, TPSMats, SOAPMats, na.rm = T)
pdf(file = file) #, width = 1200, height = 1200, res = 300)
par(mfrow = c(5, numplot), oma = c(.1, 2, 2, .1), mar = c(.1, .1, .1, .1))
for (j in seq(1, 5)) {
    lvR <- range(
                 TrueMats[[j]], Mats[[j]], TPSMats[[j]], SOAPMats[[j]], na.rm = T)
    levs = signif(seq(lvR[1], lvR[2], length.out = 10), digits = 3)
    image(set$xvec, set$yvec, TrueMats[[j]], col = fieldPal, xlab = "",
          ylab = "", zlim = zl, axes = F)
    contour(set$xvec, set$yvec, TrueMats[[j]], add = TRUE,
            zlim = zl, levels = levs)
    mtext(paste("t=", round(times[j], digits = 2), sep = ""), side = 2,
          line = 0, las = 3)
    if (j == 1) {
        mtext("True Field", font = 2)
    }
    if (whichones[1]) {
        image(set$xvec, set$yvec, Mats[[j]], col = fieldPal, xlab = "", ylab = "",
              zlim = zl, axes = F)
        contour(set$xvec, set$yvec, Mats[[j]], add = TRUE,
                zlim = zl, levels = levs)
        if (j == 1) {
            mtext("GSRt-PDE", font = 2)
        }
    }
    if (whichones[2]) {
        image(set$xvec, set$yvec, MatsSep[[j]], col = fieldPal, xlab = "",
              ylab = "", zlim = zl, axes = F)
        contour(set$xvec, set$yvec, MatsSep[[j]], add = TRUE,
                zlim = zl, levels = levs)
        if (j == 1) {
            mtext("GSR-PDE", font = 2)
        }
    }
    if (whichones[3]) {
        image(set$xvec, set$yvec, TPSMats[[j]], col = fieldPal, xlab = "",
              ylab = "", zlim = zl, axes = F)
        contour(set$xvec, set$yvec, TPSMats[[j]], add = TRUE,
                zlim = zl, levels = levs)
        if (j == 1) {
            mtext("TPS", font = 2)
        }
    }
    if (whichones[4]) {
        image(set$xvec, set$yvec, SOAPMats[[j]], col = fieldPal, xlab = "",
              ylab = "", zlim = zl, axes = F)
        contour(set$xvec, set$yvec, SOAPMats[[j]], add = TRUE,
                zlim = zl, levels = levs)
        if (j == 1) {
            mtext("SOAP", font = 2)
        }
    }
}
dev.off()
}

