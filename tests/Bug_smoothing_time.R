library(fdaPDE)

#Test copiato da smooth.FEM.time.tests.R, test #2

rm(list=ls())
graphics.off()

data(horseshoe2D)

mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)

FEMbasis = create.FEM.basis(mesh)

f<-function(x,y,t)
{
  K <- (y/0.1*as.double((abs(y)<=0.1 & x>-0.5))+as.double((abs(y)>0.1 | x<=-0.5)))^2
  res=numeric(length =length(x))
  for(i in 1:length(x))
  {
    if(x[i]>=0 && y[i]>0)
      res[i]=cos(t[i])*(0.25*pi+x[i])+(y[i]-0.5)^2
    if(x[i]>=0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-0.25*pi-x[i])+(-y[i]-0.5)^2
    if(x[i]<0 && y[i]>0)
      res[i]=cos(t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
    if(x[i]<0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
  }
  res
}

NumTimeInstants=5
TimePoints=seq(0,pi,length.out =NumTimeInstants)

space_time_locations = cbind(rep(TimePoints,each=nrow(locations)),rep(locations[,1],NumTimeInstants),rep(locations[,2],NumTimeInstants))
sol_exact = f(space_time_locations[,2],space_time_locations[,3],space_time_locations[,1])

ndata = length(sol_exact)

# Create covariates
# Il codice funziona con una sola covariata, il problema è quando ne ho 2 o più
# Il test smooth.FEM.time.tests.R ne ha sola una. 
# 
set.seed(509875)
cov1 = rnorm(ndata, mean = 1, sd = 2)
cov2 = rnorm(ndata, mean = 3, sd = 2)

# Add error to simulate data
set.seed(7893475)
data = sol_exact + 2*cov1 + cov2
data = data + rnorm(length(sol_exact), mean = 0, sd =  0.1)
observations = matrix(data,nrow(locations),NumTimeInstants)

lambdaS = 10^0
lambdaT = 10^0

output_CPP<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                            observations=observations, 
                            covariates = cbind(cov1, cov2),
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT, FLAG_PARABOLIC = T)

# La stima è sbagliata per la seconda covariata !
print(output_CPP$beta)
# Cambiando la linea 177 con : covariates=covariates[(NobsIC+1):nrow(covariates), ]
#