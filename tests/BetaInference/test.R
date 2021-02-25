library(fdaPDE)
rm(list=ls())
graphics.off()

data(horseshoe2D)

mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
locations = refine.mesh.2D(mesh, maximum_area = 0.005)$nodes
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")

FEMbasis = create.FEM.basis(mesh)

ndata = nrow(locations)

# Create covariates
set.seed(509875)
cov1 = rnorm(ndata, mean = 1, sd = 2)
cov2 = sin(locations[,1])

# Exact solution (pointwise at nodes)
DatiEsatti=fs.test(locations[,1], locations[,2]) + 0*cov1 -cov2

for(ind in 1:100)
{
  points(locations[which(round((DatiEsatti-min(DatiEsatti))/(max(DatiEsatti)-min(DatiEsatti))*100)==ind),1],
         locations[which(round((DatiEsatti-min(DatiEsatti))/(max(DatiEsatti)-min(DatiEsatti))*100)==ind),2],
         col=heat.colors(100)[ind], pch=16)
}

# Add error to simulate data
set.seed(543663)
ran = range(DatiEsatti)
data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^seq(-3,3,by=0.25)

output_CPP<-smooth.FEM(locations = locations, observations=data, 
                       covariates = cbind(cov1, cov2),
                       FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid',
                       DOF.evaluation='exact',
                       lambda.selection.lossfunction='GCV')

beta.est <- output_CPP$solution$beta
sd.est <- output_CPP$solution$estimated_sd
z.alpha <- qnorm(1-0.025)
sd.betas <- sd.est*sqrt(diag(output_CPP$solution$CovsS))
ICs <- data.frame(left=beta.est - z.alpha*sd.betas, right=beta.est + z.alpha*sd.betas)
print(ICs)

rm(list=ls())
graphics.off()

data(horseshoe2D)

mesh = create.mesh.2D(nodes=horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
locations = refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)

# plot spatial locations
plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")

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
set.seed(509875)
cov1 = rnorm(ndata, mean = 1, sd = 2)
cov2 = sin(space_time_locations[, 2])

# Add error to simulate data
set.seed(7893475)
data = sol_exact + 2*cov1 - 0*cov2
data = data + rnorm(length(sol_exact), mean = 0, sd =  0.05*diff(range(sol_exact)))
observations = matrix(data,nrow(locations),NumTimeInstants)

# Set smoothing parameter

lambdaS = 10^(-2:0)
lambdaT = 10^(-3:-1)

output_CPP<-smooth.FEM.time(locations = locations, time_mesh = TimePoints, 
                            observations=observations, 
                            covariates = cbind(cov1, cov2),
                            FEMbasis=FEMbasis, lambdaS=lambdaS, lambdaT=lambdaT, 
                            lambda.selection.criterion='grid',
                            DOF.evaluation='stochastic',
                            lambda.selection.lossfunction='GCV', FLAG_PARABOLIC=T)

beta.est <- output_CPP$beta[, output_CPP$bestlambda[1],
                             output_CPP$bestlambda[2]]

sd.est <- output_CPP$stderr[output_CPP$bestlambda[1],
                             output_CPP$bestlambda[2]]
alpha = .05
z.alpha <- qnorm(1-alpha/2)
sd.betas <- sd.est*sqrt(diag(output_CPP$CovsS)/ndata)
ICs <- data.frame(left=beta.est - z.alpha*sd.betas, right=beta.est + z.alpha*sd.betas)
print(ICs)
