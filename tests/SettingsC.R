# MESH AND DATA SETTINGS

#plot settings
{
  titlecex = .5
  labcex = .3
  pointcex = .3
  textcex = .3
  lwd = .3
}
data("horseshoe2D")

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


f <- function(x, y, t, FAMILY){
  if(FAMILY == "gamma"){
    a = 8
    b = 10  # 10
  }
  if(FAMILY == "binomial"){
    a = -1
    b = 0
  }
    return(-1/a*(fs.test(x, y, exclude = T)+b)*(t+1))
}
xvec = seq(-1, 4, 0.02)
yvec = seq(-1, 1, 0.01)

#Mesh settings
M = 5

basis_order = 1
boundary_nodes = horseshoe2D$boundary_nodes
boundary_segments = horseshoe2D$boundary_segments
locations = horseshoe2D$locations
mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments, order = basis_order)
time_mesh = seq(0,1,length.out = M)


#Setting up the datas 
m = 21
n = 600
time_locations = seq(0, 1, length.out = m)

loc = cbind(runif(n, min = -1, max = 4), runif(n, min = -1, max = 1))
ww <- ! is.na(fs.test(loc[,1], loc[ ,2], exclude = T))
loc <- loc[ww, ]
# #while(nrow(loc) != n ){
#   new <- cbind(runif(1, min = -1, max = 4), runif(1, min = -1, max = 1))
#   if(is.p.in.horseshoe(new) & ! any(sapply(loc, function(x, new) isTRUE(all.equal(x, new)), new)))
#     loc <- rbind(loc, new)
# }
space_time_locations = cbind(rep(time_locations, each=nrow(loc)), rep(loc[,1],length(time_locations)), rep(loc[,2],length(time_locations)))
FEMbasis = create.FEM.basis(mesh)
lambdaS = 10^seq(-1, 1, 0.1)
lambdaT = 10^seq(-1, 1, 0.1)

