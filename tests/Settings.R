# MESH AND DATA SETTINGS

f <- function(x, y, t, FAMILY) {
  if (sum(FAMILY == c("gamma", "exponential")) >= 1) {
    return((-1.5 * sin(2 * pi * x) * cos(2 * pi * y) +
      2 / 5 * sin(3 * pi * x * t) - 2) * (t + 1))
  }

  if (sum(FAMILY == c("binomial", "poisson")) >= 1) {
    return((-2.5 * sin(2 * pi * x) * cos(2 * pi * y) +
      4 / 5 * sin(pi * x * t)) * (t + 1))
  }
}

fC <- function(x, y, t, FAMILY, exclude = T) {
  r <- 0.5
  r0 <- .1
  l <- 3
  R <- x # rep(0, length(x))
  R[y >= 0 & x >= 0] <- y[y >= 0 & x >= 0] - r
  R[x < 0 & y >= 0] <- sqrt((x[x < 0 & y >= 0])^2 + (y[x < 0 & y >= 0])^2) - r
  R[x < 0 & y < 0] <- sqrt((x[x < 0 & y < 0])^2 + (y[x < 0 & y < 0])^2) - r
  R[y < 0 & x >= 0] <- -y[y < 0 & x >= 0] - r
  R[R < -r + r0 | R > r - r0] <- NA
  a <- -1
  b <- 0
  if (sum(FAMILY == c("gamma", "exponential")) >= 1) {
    a <- 8
    b <- 10 # 10
  }
  if(FAMILY == "exponential"){
      a  <- 2
      b  <- 10
  }
  if (FAMILY == "poisson") a <- -2
  K <- (y / 0.1 * as.double((abs(y) <= 0.1 & x > -0.5)) +
    as.double((abs(y) > 0.1 | x <= -0.5)))^2
  res <- numeric(length = length(x))
  for (i in 1:length(x)) {
    if (x[i] >= 0 && y[i] > 0) {
      res[i] <- cos(-t[i]) * (0.25 * pi + x[i]) + (y[i] - 0.5)^2
    }
    if (x[i] >= 0 && y[i] <= 0) {
      res[i] <- cos(-2 * t[i]) * (-0.25 * pi - x[i]) + (-y[i] - 0.5)^2
    }
    if (x[i] < 0 && y[i] > 0) {
      res[i] <- cos(-t[i]) * (-atan(y[i] / x[i]) * 0.5) +
        (sqrt(x[i]^2 + y[i]^2) - 0.5)^2 * K[i]
    }
    if (x[i] < 0 && y[i] <= 0) {
      res[i] <- cos(-2 * t[i]) * (-atan(y[i] / x[i]) * 0.5) +
        (sqrt(x[i]^2 + y[i]^2) - 0.5)^2 * K[i]
    }
    res[i] <- -1 / a * (res[i] + b)
  }
  if (exclude) {
    res[is.na(R) | abs(R) > r - r0 | (x > l & (x - l)^2 + R^2 > (r - r0)^2)] <- NA
  }
  return(res)
}
inv.link <- function(mu, FAMILY) {
  if (FAMILY == "gamma") {
    return(-1 / mu)
  }
  if (FAMILY == "binomial") {
    return(plogis(mu))
  }
  if (FAMILY == "gaussian") {
    return(mu)
  }
  if (FAMILY == "Gaussian") {
    return(mu)
  }
  if (FAMILY == "poisson") {
    return(exp(mu))
  }
  if (FAMILY == "exponential") {
    return(-1 / mu)
  }
}
settings <- function(flagMeshC = T) {
  Settings <- NULL

  Settings$N <- 16
  Settings$M <- 11

  Settings$m <- Settings$M
  Settings$n <- 400

  Settings$PDE_parameters <- NULL
  Settings$scale <- .1
  Settings$inv.link <- inv.link
  if (flagMeshC) Settings$f <- fC else Settings$f <- f
  if (flagMeshC) {
    data("horseshoe2D")
    Settings$xvec <- seq(-1, 3.5, .03)
    Settings$yvec <- seq(-1, 1, .03)
  } else {
    Settings$xvec <- seq(0, 1, by = 0.03)
    Settings$yvec <- seq(0, 1, by = 0.03)
  }

  # Mesh settings

  Settings$basis_order <- 1

  # discrete meshes
  if (flagMeshC) {
    boundary_nodes <- horseshoe2D$boundary_nodes
    boundary_segments <- horseshoe2D$boundary_segments
    Settings$mesh <- create.mesh.2D( boundary_nodes, segments = boundary_segments, order = Settings$basis_order)
    Settings$mesh <- refine.mesh.2D(Settings$mesh, minimum_angle=30, maximum_area=0.025)
    dimeval <- 10000
    xeval <- seq(-1, 3.5, 0.05)
    yeval <- seq(-1, 1, 0.05)
    xyeval <- expand.grid(xeval, yeval)
    Settings$evalGrid <- data.frame(
      t = rep(seq(0, 1, length.out = 20), each = nrow(xyeval)),
      x = rep(xyeval[, 1], 20), y = rep(xyeval[, 2], 20)
    )
  } else {
    tmp <- cbind(
      rep(seq(0, 1, length.out = Settings$N), Settings$N),
      rep(seq(0, 1, length.out = Settings$N), each = Settings$N)
    )
    boundary_nodes <- tmp[which(tmp[, 1] == 0 |
      tmp[, 1] == 1 | tmp[, 2] == 0 | tmp[, 2] == 1), ]
    inner_nodes <- tmp[-which(tmp[, 1] == 0 |
      tmp[, 1] == 1 | tmp[, 2] == 0 | tmp[, 2] == 1), ]
    Settings$mesh <- create.mesh.2D(
      nodes = rbind(boundary_nodes, inner_nodes),
      order = Settings$basis_order
    )
    dimeval <- 1000
    Settings$evalGrid <- data.frame(
      t = runif(n = dimeval, 0, 1), x = runif(n = dimeval, 0, 1),
      y = runif(n = dimeval, 0, 1)
    )
  }

  Settings$time_mesh <- seq(0, 1, length.out = Settings$M)

  # Setting up the datas
  Settings$time_locations <-
    seq(0, 1, length.out = Settings$m) # c(0, sort(runif(Settings$m, 0.0001,
  # 1))) #
  if (flagMeshC) {
    Settings$loc <- cbind(
      runif(2 * Settings$n, min = -1, max = 4),
      runif(2 * Settings$n, min = -1, max = 1)
    )
    ww <- apply(Settings$loc, 1, is.p.in.horseshoe) # ! is.na(fs.test(loc[,1],
    # loc[ ,2], exclude = T))
    Settings$loc <- Settings$loc[ww, ]
    if (nrow(Settings$loc) > Settings$n) {
      Settings$loc <- Settings$loc[1:Settings$n, ]
    }
    Settings$lambdaS <- 10^seq(-4, 2, 1)
    Settings$lambdaT <- 10^seq(-3, 1, 1)

    Settings$lambdaSs <- 10^seq(-4, 2, 1)
    Settings$lambdaTs <- 10^seq(-3, 0, 1)
  } else {
    Settings$loc <- cbind(runif(Settings$n), runif(Settings$n))
    Settings$lambdaS <- 10^seq(-4, 3, .5)
    Settings$lambdaT <- 10^0 # seq(-3, 0, 1)

    Settings$lambdaSs <- 10^seq(-2, 0, 1)
    Settings$lambdaTs <- 10^seq(-3, -1, 1)
  }
  Settings$space_time_locations <-
    cbind(
      rep(Settings$time_locations, each = nrow(Settings$loc)),
      rep(Settings$loc[, 1], length(Settings$time_locations)),
      rep(Settings$loc[, 2], length(Settings$time_locations))
    )
  Settings$FEMbasis <- create.FEM.basis(Settings$mesh)
  Settings$n <- nrow(Settings$loc)
  # Settings$loc = NULL
  if (flagMeshC) {
    Settings$dir_name <-
      paste("P", 2 * Settings$basis_order, "_", "M", Settings$M, "_", "n",
        Settings$n, "_", "m", Settings$m,
        sep = ""
      )
    if (!dir.exists(Settings$dir_name)) {
      dir.create(Settings$dir_name)
      dir.create(paste(Settings$dir_name, "/figures", sep = ""))
    }

    boundaryidx <- which(Settings$mesh$nodesmarkers == 1)
    boundary_nodes <- Settings$mesh$nodes[boundaryidx, ]
    boundary <- rbind(boundary_nodes, boundary_nodes[1, ])
    boundary <- boundary[nrow(boundary):1, ]
    Settings$fsb <- list(list(boundary[, 1], boundary[, 2]))
    names(Settings$fsb[[1]]) <- c("x", "y")
    Settings$nmax <- 100 # FEMbasis$nbasis
    # Settings$knots =
    #    data.frame(
    #               x = Settings$mesh$nodes[boundaryidx, 1],
    #               y = Settings$mesh$nodes[boundaryidx, 2])
    Settings$knots <- data.frame(
      x = rep(seq(-.5, 3, by = .5), 4),
      y = rep(c(-.6, -.3, .3, .6), rep(8, 4))
    )
  } else {
    Settings$dir_name <-
      paste("P", 2 * Settings$basis_order, "_", "N", Settings$N, "_", "M",
        Settings$M, "_", "n", Settings$n, "_", "m", Settings$m,
        sep = ""
      )
    if (!dir.exists(Settings$dir_name)) {
      dir.create(Settings$dir_name)
      dir.create(paste(Settings$dir_name, "/figures", sep = ""))
    }

    boundary <- Settings$mesh$nodes[which(Settings$mesh$nodesmarkers == 1), ]
    idx <- order(atan2(boundary[, 2] - .001, boundary[, 1] - .001))
    boundary <- rbind(boundary[idx, ], c(0, 0))
    Settings$fsb <- list(list(boundary[, 1], boundary[, 2]))
    names(Settings$fsb[[1]]) <- c("x", "y")
    Settings$nmax <- 100
    Settings$knots <-
      data.frame(
        "x" = Settings$mesh$nodes[which(Settings$mesh$nodesmarkers ==
          0), 1],
        "y" = Settings$mesh$nodes[which(Settings$mesh$nodesmarkers ==
          0), 2]
      )
  }
  return(Settings)
}

is.in.horseshoe <- function(x, y) {
  r <- .5
  r0 <- .1

  if ((x - 3)^2 + (y - r)^2 < (r - r0)^2 & x > 3) {
    return(TRUE)
  }
  if (r0^2 <= x^2 + y^2 & x^2 + y^2 <= (2 * r - r0)^2 & x <= 0) {
    return(TRUE)
  }
  if ((x - 3)^2 + (y + r)^2 < (r - r0)^2 & x > 3) {
    return(TRUE)
  }
  if (abs(y) > r0 & abs(y) < 2 * r - r0 & 0 < x & x <= 3) {
    return(TRUE)
  }
  return(FALSE)
}
is.p.in.horseshoe <- function(p) {
  return(is.in.horseshoe(p[1], p[2]))
}
