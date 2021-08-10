#2-flank SPIM model from Augustine et al. 2018.
#This implements the algorithm from this paper that was originally in 
#an R-based sampler in nimble. The ID update has been simplified substantially
#and is now more efficient.

#The observation model is exactly the same as in the paper, with 
#both-side caps only possible at double cam stations
#and single side cap prob depending on whether a station has 1 or 2 cams.

#This sampler works for all 1 cam stations, all 2 cam stations, or a mixture.
#The number of cams operable can change across occasions
#You can leave the both-side capture prob in the model without any 2 cam stations
#and estimation is still correct (though very slightly slower), p0B is just updated from the prior.

#nimble model set up to estimate p0L and p0R as a pooled parameter, but they can be estimated 
#separately by switching the commented lines in the model code. Will need to change 

#Code would be faster if rewritten using the 2D M x J data instead of the 3D M x J x K data,
#but using 3D data for maximal flexibility (can add occasion effects, etc.)

#To add occasion-level effects, "pd" needs to be switched to 3D (adding occasion dimension) and
#the custom ID updates need to be modified to reflect that.

#Pay very close attention to the "n.fixed" part of the code. This allows for known flank matches of individuals
#never captured by both flanks simultaneously. Otherwise, we assume we only know the flank matches of n.B 
#individuals, the number captured in the both-flank detection process.
#Possibilities are
#1) At least 1 both-side captures and no extra known flank matches
#2) At least 1 both-side captures and some additional known flank matches
#3) No both-side captures and no known flank matches
#4) No both-side captures and some known flank matches
#Must get this correct to ensure data are simulated and initialized correctly.

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
buff=2 #should be 3sigma or greater (at least for real data, unless using a habitat mask)
X<- expand.grid(3:8,3:8) #make some traps
J=nrow(X)
#J x K indicator for number of cameras operable 
K2D=matrix(1,nrow=J,ncol=K) #start with all single cam stations
# K2D[c(1,5,10,15,20,25,30,35),1:10]=2 #create some dual cam stations (or could leave all single)
#If you have no trap-occasions with 2 cams, data simulator assumes no flanks can be deterministically linked.

#You may know some flank linkages without both-side captures. If NA, we assume we do not.
# n.fixed=5 #minimum number of matched flanks for simulation
n.fixed=NA #supply NA for only n.B flank matches to be known (if there are any).
#must supply this for simulation so we know to retain left and right known flanks with no captures
data=sim2flank(N=N,p0L=p0L,p0R=p0R,p0B=p0B,sigma=sigma,K=K,X=X,buff=buff,K2D=K2D,n.fixed=n.fixed)
n.fixed=data$n.fixed #update in case n.B > n.fixed (n.B is number of both-side captured individuals)
#*Do not change* n.fixed from here down because all 0 left and/or right capture histories are retained
# for n.fixed individuals

#True individual numbers for each left and right flank. If they match,
#the left and right flanks belong to the same individual. These must be
#probabilistically associated during MCMC.
data$ID.L
data$ID.R

#Observed both flank capture data is complete
str(data$y.B.obs)

#Left and right observed data. In correct order for 1:n.B individuals but not
#for individuals never captured on both flanks. Can enforce that left only and right
#only flanks are deterministically linked if known by other means than a both-side capture
str(data$y.L.obs)
str(data$y.R.obs)

M=100 #data augmentation level. If N posterior hits M, raise M and start over.

#initialize latent states, including flank matches.
#Can initialize to truth for simulated data for testing purposes.
#Otherwise, left and right flanks are initialized by 
#crudely minimizing the distance between the mean left and right
#side capture locations. 
#If n.fixed is NA, will fix flanks for n.B individuals if there are any, will not fix any otherwise
#plug in n.fixed here if different from n.B
nimbuild=init.2flank(data=data,M=M,n.fixed=n.fixed,initTrue=FALSE)

#crappy plot
par(mfrow=c(1,1),ask=FALSE)
plot(X,pch=4,xlim=nimbuild$xlim,ylim=nimbuild$ylim)
points(nimbuild$s[nimbuild$z==1,],pch=16)
y2D=apply(nimbuild$y.true,c(1,2),sum)
for(i in 1:M){
  trapcaps=which(y2D[i,]>0)
  for(j in 1:length(trapcaps)){
    lines(x=c(nimbuild$s[i,1],X[trapcaps[j],1]),y=c(nimbuild$s[i,2],X[trapcaps[j],2]))
  }
}


#inits for nimble
#use s inits in nimbuild! May not converge if you let nimble do it.
Niminits <- list(z=nimbuild$z,s=nimbuild$s,ID.L=nimbuild$ID.L,ID.R=nimbuild$ID.R,
                 y.L.true=nimbuild$y.true[,,,2],
                 y.R.true=nimbuild$y.true[,,,3],
                 sigma=sigma,p0S=0.13,psi=0.5,p0B=0.5) #let nimble initialize p0s, don't use true sigma as init, but ballpark start speeds convergence

#constants for Nimble
n.L=nrow(data$y.L.obs)
n.R=nrow(data$y.R.obs)
constants<-list(M=M,J=J,K=K,K2D=K2D,xlim=nimbuild$xlim,ylim=nimbuild$ylim,n.L=n.L,n.R=n.R)

# Supply data to Nimble. Only individuals with deterministially-matched flanks have known z's.
if(n.fixed==0){
  z.data=rep(NA,M)
}else{
  z.data=c(rep(1,n.fixed),rep(NA,M-n.fixed))
}
#left and right flank data treated as fully latent
#Known flank matches will not be updated during MCMC.
Nimdata<-list(y.B.true=nimbuild$y.true[,,,1],y.L.true=array(NA,dim=c(M,J,K)),y.R.true=array(NA,dim=c(M,J,K)),
              ID.L=rep(NA,n.L),ID.R=rep(NA,n.R),z=z.data,X=as.matrix(X))

# set parameters to monitor
# parameters<-c('psi','p0B','p0L','p0R','sigma','N') #left and right side cap prob not shared
parameters<-c('psi','p0B','p0S','sigma','N')
parameters2<-c("ID.L","ID.R")

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
#thinning parameters 2 more to reduce posterior size
conf <- configureMCMC(Rmodel,monitors=parameters,monitors2=parameters2,thin=1,thin2=5,useConjugacy = TRUE)

###*required* sampler replacements
conf$removeSampler("y.L.true")
conf$removeSampler("y.R.true")
conf$addSampler(target = paste0("y.L.true[1:",M,",1:",J,",1:",K,"]"),
                type = 'IDLSampler',control = list(K2D=K2D,n.fixed=n.fixed,n.L=n.L,prop.scale=1),silent = TRUE)
conf$addSampler(target = paste0("y.R.true[1:",M,",1:",J,",1:",K,"]"),
                type = 'IDRSampler',control = list(K2D=K2D,n.fixed=n.fixed,n.R=n.R,prop.scale=1),silent = TRUE)
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
#you can put AF slice on the p0x and sigma intercepts. Don't replace p0B sampler if no double cam traps
# conf$removeSampler(c("p0S","sigma"))
# conf$addSampler(target = c(p0S","sigma"),
#                 type = 'AF_slice',
#                 control = list(adaptive=TRUE),
#                 silent = TRUE)
#Block RW faster than AF, but mixing worse. May be most efficient?
conf$removeSampler(c("p0S","sigma"))
conf$addSampler(target = c("p0S","sigma"),
                type = 'RW_block',
                control = list(adaptive=TRUE),
                silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for easier debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model
start.time2<-Sys.time()
#short run for demonstration
Cmcmc$run(5000,reset=FALSE) #can continue MCMC by rerunning this line with reset=FALSE
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

#extract and plot posterior samples
library(coda)
mvSamples=as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))


#posterior sample match probs. not removing burnin here
mvSamples2=as.matrix(Cmcmc$mvSamples2)
postL=mvSamples2[,1:n.L] #which 1:M individual is left flank matched to on each iteration>
postR=mvSamples2[,(n.L+1):(n.R+n.L)] #same for right flanks

#Calculate posterior prob that left flank l matches right flank r
postprobs=matrix(NA,nrow=n.L,ncol=n.R)
for(l in 1:n.L){
  for(r in 1:n.R){
    tmp=postL[,l]==postR[,r]
    postprobs[l,r]=mean(tmp[-1])
  }
}

#Posterior match probability that left flanks are matched to correct right flanks (for simulated data)
#Higher match probs indicate more certainty in matches and less precision lost relative to true data
#**match probs should be 1 for 1:n.fixed individuals. If not, something went wrong.**
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

#posterior for number of matching flank pairs
n.match=rep(NA,nrow(postL)-1)
for(i in 2:nrow(postL)){
  n.match[i]=sum(postL[i,]%in%postR[i,])
}
n.match
sum(data$ID.L%in%data$ID.R)
