#2-flank SPIM model from Augustine et al. 2018.
#This implements the algorithm from this paper that was originally in 
#an R-based sampler in nimble. The ID update has been simplified substantially
#and is now more efficient.

#The observation model is exactly the same as in this paper, with 
#both side caps only possible at double cam stations
#and single side caps depending on whether a station has 1 or 2 cams.

#this sampler works for all 1 cam stations, all 2 cam stations, or a mixture.
#number of cams operable can change across occasions
#you can leave the both side capture prob in the model without any 2 cam stations
#and estimation is still correct (though very slightly slower), p0B is just updated from prior.

#nimble model set up to estimate p0L and p0R separately, but can pool the flanks by switching the 
#commented lines in the model code. Will need to change 

#Code would be faster if rewritten using the 2D M x J data instead of the 3D M x J x K data,
#but using 3D data for maximal flexibility (can add occasion effects, etc.)

source("sim2flank.R")
source("init.2flank.R")
source("LRmatch.R")
library(nimble)
source("NimbleModel 2Flank.R")
source("NimbleFunctions 2Flank.R")

N=50
p0L=0.13
p0R=0.13
p0B=0.2
sigma=0.50
K=10
buff=2 #should be 3sigma or greater
X<- expand.grid(3:8,3:8) #make some traps
J=nrow(X)
#J x K indicator for number of cameras operable 
K2D=matrix(1,nrow=J,ncol=K) #start with all single cam stations
K2D[c(1,5,10,15,20),c(2,4,6,8)]=2 #create some dual cam stations (or could leave all single)

data=sim2flank(N=N,p0L=p0L,p0R=p0R,p0B=p0B,sigma=sigma,K=K,X=X,buff=buff,K2D=K2D)

#True individual numbers for each left and right flank. These must be
#probabilistically associated during MCMC
data$ID.L
data$ID.R

M=100 #data augmentation level. If N posterior hits M, raise M and start over.

#Can initialize to truth for simulated data for testing.
nimbuild=init.2flank(data=data,M=M,initTrue=TRUE)

#We must treat either the left or right flanks as latent
#Code below treats left flank IDs as latent, both and right as known.

#inits for nimble
#use s inits in nimbuild! May not converge if you let nimble do it.
Niminits <- list(z=rep(1,M),s=nimbuild$s,ID.L=nimbuild$ID.L,ID.R=nimbuild$ID.R,
                 y.L.true=nimbuild$y.true[,,,2],#need inits bc latent
                 sigma=sigma*2) #let nimble initialize p0s, don't use true sigma as init, but ballpark start speeds convergence

#constants for Nimble
n.L=nrow(data$y.L.obs)
n.R=nrow(data$y.R.obs)
constants<-list(M=M,J=J,K=K,K2D=K2D,xlim=nimbuild$xlim,ylim=nimbuild$ylim,n.L=n.L,n.R=n.R)

# Supply data to Nimble. Note, y.true is completely latent.
n.B=sum(rowSums(data$y.B.obs)>0)
if(n.B==0){
  z.data=rep(NA,M)
}else{
  z.data=c(rep(1,n.B),rep(NA,M-n.B))
}
#left flank data latent
Nimdata<-list(y.B.true=nimbuild$y.true[,,,1],y.L.true=array(NA,dim=c(M,J,K)),y.R.true=nimbuild$y.true[,,,3],
              ID.L=rep(NA,n.L),z=z.data,X=as.matrix(X))

# set parameters to monitor
parameters<-c('psi','p0B','p0L','p0R','sigma','N')
parameters2<-c("ID.L")

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters,monitors2=parameters2,thin=1,thin2=5,useConjugacy = TRUE)

###*required* sampler replacement for left flank IDs
conf$removeSampler("y.L.true")
conf$addSampler(target = paste0("y.L.true[1:",M,",1:",J,",1:",K,"]"),
                type = 'IDLSampler',control = list(K2D=K2D,n.fixed=n.B,n.L=n.L,prop.scale=2),silent = TRUE)

###*optional* sampler replacements:
#Update activity X and Y locations at the same time. Longer adapt interval than default may be helpful since
#data assigned to each activity center is constantly updating. Hard to tune.
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""), #do not adapt covariance bc s's not deterministically linked to individuals
                  type = 'RW_block',control=list(adaptive=TRUE,adaptScaleOnly=TRUE,adaptInterval=500),silent = TRUE)
}

#AF slice update can improve mixing substantially when strong covariance between p0x and sigma
#but it is slower per iter. Typically more useful with sparse data sets? With detection covariates, 
#you can put AF slice on the p0x and sigma intercepts.
# conf$removeSampler(c("p0B","p0L","p0R","sigma"))
# conf$addSampler(target = c("p0B","p0L","p0R","sigma"),
#                 type = 'AF_slice',
#                 control = list(adaptive=TRUE),
#                 silent = TRUE)
#Block RW faster than AF, but mixing worse
# conf$removeSampler(c("p0B","p0L","p0R","sigma"))
# conf$addSampler(target = c("p0B","p0L","p0R","sigma"),
#                 type = 'RW_block',
#                 control = list(adaptive=TRUE),
#                 silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model
start.time2<-Sys.time()
Cmcmc$run(2000,reset=FALSE) #can continue MCMC by rerunning this line with reset=FALSE
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples=as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))


#posterior sample match probs
mvSamples2=as.matrix(Cmcmc$mvSamples2)
postL=mvSamples2[,1:n.L]

#Calculate posterior prob that left flank l matches right flank r
postprobs=matrix(NA,nrow=n.L,ncol=n.R)
for(l in 1:n.L){
  for(r in 1:n.R){
    tmp=postL[,l]==nimbuild$ID.R[r]
    postprobs[l,r]=mean(tmp[-1])
  }
}

#Posterior match probability that left flanks are matched to correct right flanks (for simulated data)
for(i in 1:n.L){
  idx=which(data$ID.R==data$ID.L[i])
  if(length(idx)>0){
    print(postprobs[i,idx])
  }else{
    print(NA) #was no matching right flank
  }
}
#Posterior match probability that right flanks are matched to correct left flanks (for simulated data)
for(i in 1:n.R){
  idx=which(data$ID.L==data$ID.R[i])
  if(length(idx)>0){
    print(postprobs[idx,i])
  }else{
    print(NA) #was no matching left flank
  }
}
