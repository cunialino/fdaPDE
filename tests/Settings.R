# MESH AND DATA SETTINGS

f <- function(x, y, t, FAMILY){
  if(FAMILY == "gamma")
    return((-1.5*sin(2*pi*x)*cos(2*pi*y) + 2/5*sin(3*pi*x*t)-2)*(t+1))
  if(FAMILY == "binomial")
    return((-2.5*sin(2*pi*x)*cos(2*pi*y) + 4/5*sin(pi*x*t))*(t+1))
}

fC <- function(x, y, t, FAMILY){
  if(FAMILY == "gamma"){
    a = 8
    b = 10  # 10
  }
  if(sum(FAMILY == c("binomial", "gaussian")) >= 1){
    a = -3
    b = 0
  }
  return(-1/a*(fs.test(x, y, exclude = T)+b - sin(pi*t))*(t+1))
}
inv.link <- function(mu, FAMILY){
  if(FAMILY == "gamma"){
    return(-1/mu)
  }
  if(FAMILY == "binomial"){
    return(plogis(mu))
  }
  if(FAMILY == "gaussian")
    return(mu)
}
settings <- function(flagMeshC = T){
  
  Settings = NULL
  Settings$scale = .05
  Settings$inv.link <- inv.link
  if(flagMeshC)
    Settings$f = fC
  else
    Settings$f = f
  if(flagMeshC){
    data("horseshoe2D")
    Settings$xvec = seq(-1, 4, 0.02)
    Settings$yvec = seq(-1, 1, 0.01)
  }
  else{
    Settings$xvec = seq(0, 1, by = 0.01)
    Settings$yvec = seq(0, 1, by = 0.01)
  }
  
  #Mesh settings
  Settings$N = 16
  Settings$M = 10
  
  Settings$basis_order = 1
  
  #discrete meshes
  if(flagMeshC){
    
    boundary_nodes = horseshoe2D$boundary_nodes
    boundary_segments = horseshoe2D$boundary_segments
    locations = horseshoe2D$locations
    Settings$mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments, order = Settings$basis_order)
  }
  else{
    tmp <- cbind(rep(seq(0, 1, length.out = Settings$N), Settings$N), rep(seq(0, 1, length.out = Settings$N), each = Settings$N))
    boundary_nodes = tmp[which(tmp[, 1] == 0 | tmp[, 1] == 1 | tmp[, 2] == 0 | tmp[, 2] == 1), ]
    inner_nodes =  tmp[-which(tmp[, 1] == 0 | tmp[, 1] == 1 | tmp[, 2] == 0 | tmp[, 2] == 1), ]
    Settings$mesh = create.mesh.2D(nodes = rbind(boundary_nodes, inner_nodes), order = Settings$basis_order)
  }
  
  Settings$time_mesh = seq(0, 1, length.out = Settings$M)
  
  #Setting up the datas 
  Settings$m = 41 #Settings$M
  Settings$n = 600
  Settings$time_locations = seq(0, 1, length.out = Settings$m)
  if(flagMeshC){
    loc = cbind(runif(2*Settings$n, min = -1, max = 4), runif(2*Settings$n, min = -1, max = 1))
    ww <- apply(loc, 1, is.p.in.horseshoe) #! is.na(fs.test(loc[,1], loc[ ,2], exclude = T))
    Settings$loc <- loc[ww, ]
    if(nrow(loc) > Settings$n){
      Settings$loc = Settings$loc[1:Settings$n, ]
    }
  }
  else{
    Settings$loc = cbind(runif(Settings$n), runif(Settings$n))
  }
  Settings$space_time_locations = cbind(rep(Settings$time_locations, each=nrow(Settings$loc)), rep(Settings$loc[,1],length(Settings$time_locations)), rep(Settings$loc[,2],length(Settings$time_locations)))
  Settings$FEMbasis = create.FEM.basis(Settings$mesh)
  Settings$lambdaS = 10^0 #seq(0, 1, 0.1)
  Settings$lambdaT = 10^0 #seq(-1, 0, 0.1)
  
  Settings$lambdaSs = 10^0 #seq(1, 1.5, 0.1)
  Settings$lambdaTs = 10^0 #seq(-3, -2.5, 0.1)
  if(flagMeshC){
    Settings$dir_name = paste("P", 2*Settings$basis_order, "_", "M", Settings$M, "_", "n", Settings$n, "_", "m", Settings$m, sep = "")
    if(! dir.exists(Settings$dir_name)){
      dir.create(Settings$dir_name)
      dir.create(paste(Settings$dir_name, "/figures", sep = ""))
    }
    
    boundary = rbind(boundary_nodes, boundary_nodes[1, ])
    boundary <- boundary[nrow(boundary):1, ]
    Settings$fsb=list(list(boundary[,1],boundary[,2]))
    names(Settings$fsb[[1]]) = c("x","y")
    Settings$nmax = 100 #FEMbasis$nbasis
    Settings$knots = data.frame("x"=Settings$mesh$nodes[which(Settings$mesh$nodesmarkers == 0), 1], "y"=Settings$mesh$nodes[which(Settings$mesh$nodesmarkers == 0), 2])
    boundaryidx <- which(Settings$mesh$nodesmarkers == 1)
  }
  else{
    Settings$dir_name = paste("P", 2*Settings$basis_order, "_", "N", Settings$N, "_", "M", Settings$M, "_", "n", Settings$n, "_", "m", Settings$m, sep = "")
    if(! dir.exists(Settings$dir_name)){
      dir.create(Settings$dir_name)
      dir.create(paste(Settings$dir_name, "/figures", sep = ""))
    }
    
    boundary = Settings$mesh$nodes[which(Settings$mesh$nodesmarkers == 1), ]
    idx <- order(atan2(boundary[, 2]-.001, boundary[, 1]-.001))
    boundary = rbind(boundary[idx, ], c(0, 0))
    Settings$fsb=list(list(boundary[,1],boundary[,2]))
    names(Settings$fsb[[1]]) = c("x","y")
    Settings$nmax = 100
    Settings$knots = data.frame("x"=Settings$mesh$nodes[which(Settings$mesh$nodesmarkers == 0), 1], "y"=Settings$mesh$nodes[which(Settings$mesh$nodesmarkers == 0), 2])
  }
  # save(list = c("Settings"), file = paste(Settings$dir_name, "/Settings.RData", sep = ""))
  return(Settings)
}

is.in.horseshoe <- function(x, y){
  r = .5
  r0 = .1
  
  if((x-3)^2 + (y-r)^2 < (r-r0)^2 & x > 3)
    return(TRUE)
  if ( r0^2 <= x^2+y^2 & x^2+y^2 <= (2*r-r0)^2 & x <= 0)
    return(TRUE)
  if((x-3)^2+(y+r)^2 < (r-r0)^2 & x > 3)
    return(TRUE)
  if(abs(y) > r0 & abs(y) < 2*r-r0 & 0 < x & x <= 3)
    return(TRUE)
  return(FALSE)
  
}
is.p.in.horseshoe <- function(p){
  return(is.in.horseshoe(p[1], p[2]))
}
