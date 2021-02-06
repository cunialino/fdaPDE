library(fdaPDE)
library(RColorBrewer)
library(pracma)
library(sp)
library(polyCub)

rm(list = ls())

source("Settings.R")
source("Plots.R")
source("Utils.R")

set <- settings(T)
set$scale <- .1
set$FAMILY <- "gamma"
load("areal.RData")
mesh <- set$mesh
plot(mesh)
title(nrow(mesh$nodes))

nodes <- mesh$nodes
set$time_locations <- set$time_locations * pi
set$time_mesh <- set$time_mesh * pi
set$space_time_locations[, 1] <- set$space_time_locations[, 1] * pi
set$evalGrid$t <- set$evalGrid$t * pi
set$lambdaS <- 10^seq(-6, 2, 0.5)
set$lambdaT <- 1

tri <- set$mesh$triangles
nodi_plot <- NULL
polynodes <- NULL
# If uncommented, every triangle will be a region, otherwise will use the incidence
# matrix in the "areal.RData" file
#Build 10 random reagions of 10 triangles each

incidence_matrix <- diag(nrow(tri))
for (j in 1:ncol(incidence_matrix)) {
  trigs <- which(incidence_matrix[, j] == 1)

  # nodes in the polygon
  nds <- unique(as.vector(tri[trigs, ]))
  # coordinates of the polgon
  cds <- mesh$nodes[nds, ]
  # baricenter of the polygon
  bcs <- colMeans(cds)
  # antclock wise ordering (needed to have the right area in the polygon sp object)
  polynodes[[j]] <- nds[order(atan2(cds[, 2] - bcs[2], cds[, 1] - bcs[1]))]
}

fff <- function(p, t, FAMILY) {
  set$f(p[, 1], p[, 2], rep(t, nrow(p)), FAMILY, exclude = F)
}
polys <- lapply(polynodes, function(z) {
  Polygon(coords = mesh$nodes[c(z, z[1]), ])
})
int <- matrix(NA, nrow = length(polys), ncol = length(set$time_mesh))
areas <- unlist(lapply(polys, function(z) z@area))
for (j in 1:length(set$time_mesh)) {
  int[, j] <- unlist(lapply(polys, function(z) {
    polyCub(z,
      f = fff, method = "SV",
      t = set$time_mesh[j],
      FAMILY = set$FAMILY
    )
  })) / areas
}
desmat <- NULL # rbeta(n=length(int), shape1 = 2, shape2 = 2)
set$loc <- matrix(unlist(lapply(polys, function(z) z@labpt)),
  nrow = length(polys), ncol = 2, byrow = T
)
field <- matrix(int, nrow = length(int), ncol = 1) # - desmat * 3
NSIM <- 5
avg.coeffs <- NULL
RMSE <- NULL
for (sim in 1:NSIM) {
  print(sim)
  data <- rgamma(
    n = length(field), shape = -1 / field / set$scale,
    scale = set$scale
  )
  data <- matrix(data, nrow = length(polys), ncol = set$m)
  storage.mode(data) <- "double"
  sims <- NULL
  sims$GSRtPDE <- smooth.FEM.time(
    time_mesh = set$time_mesh, time_locations = set$time_locations,
    observations = data, FEMbasis = set$FEMbasis, covariates = desmat,
    DOF.evaluation = "stochastic", lambda.selection.lossfunction = "GCV",
    lambdaS = set$lambdaS, lambdaT = set$lambdaT, max.steps.FPIRLS = 15,
    family = set$FAMILY, mu0 = NULL, scale.param = NULL, FLAG_PARABOLIC = T,
    threshold.FPIRLS = 10^-6, incidence_matrix = t(incidence_matrix),
    areal.data.avg = T
  )

  sims$GSRPDE <- smooth.FEM.time(
    time_mesh = set$time_mesh, time_locations = set$time_locations,
    observations = data, FEMbasis = set$FEMbasis, covariates = desmat,
    DOF.evaluation = "stochastic", lambda.selection.lossfunction = "GCV",
    lambdaS = set$lambdaS, lambdaT = set$lambdaT, max.steps.FPIRLS = 15,
    family = set$FAMILY, mu0 = NULL, scale.param = NULL, FLAG_PARABOLIC = F,
    threshold.FPIRLS = 10^-6, incidence_matrix = t(incidence_matrix),
    areal.data.avg = T
  )
  R <- eval.rmse(sims, set, covariates = !is.null((set$betas)))
  RMSE$GSRtPDE <- c(RMSE$GSRtPDE, R$GSRtPDE)
  avg.coeffs$GSRtPDE[[
  sim]] <- sims$GSRtPDE$fit.FEM.time$coeff[
    ,
    sims$GSRtPDE$bestlambda[
      1
    ],
    sims$GSRtPDE$bestlambda[
      2
    ]
  ]
  RMSE$GSRPDE <- c(RMSE$GSRPDE, R$GSRPDE)
  avg.coeffs$GSRPDE[[
  sim]] <- sims$GSRPDE$fit.FEM.time$coeff[
    ,
    sims$GSRPDE$bestlambda[
      1
    ],
    sims$GSRPDE$bestlambda[
      2
    ]
  ]
}
avg.coeffs <- lapply(avg.coeffs, function(x) {
  Reduce("+", x) / length(x)
})
if (!is.null(avg.coeffs$GSRtPDE)) {
  sims$GSRtPDE$fit.FEM.time <- FEM.time(
    avg.coeffs$GSRtPDE, set$time_mesh, set$FEMbasis, T
  )
}
if (!is.null(avg.coeffs$GSRPDE)) {
  sims$GSRPDE$fit.FEM.time <- FEM.time(
    avg.coeffs$GSRPDE, set$time_mesh, set$FEMbasis, F
  )
}

plot.settings(data, set, "arealset.pdf")
plot.results(sims, set, "areal.pdf")
pdf("RMSEareal.pdf") # , width = 1200, height = 1200, res = 300)
boxplot(list(RMSE$GSRtPDE, RMSE$GSRPDE), names = NULL)
dev.off()

pdf("meshAreal.pdf") # , width = 1200, height = 1200, res = 300)
plot(mesh, lwd = 3, cex = 1.9)
for (j in 1:dim(tri)[1]) {
  coords <- mesh$nodes[tri[j, ], ]
  reg <- which(incidence_matrix[j, ] == 1)
  lab <- j
  polygon(coords[, 1], coords[, 2],
    col = brewer.pal(11, "Spectral")[reg[1] %% 11 + 1]
  )
  # text(x=set$loc[j, 1], y=set$loc[j, 2], labels = lab)
}
dev.off()
