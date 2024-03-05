#this version allows a discrete habitat mask
#AND only uses the 2D individual x trap data omitting the trap dimension
#this can be much faster if you don't consider time effects on detection
#and camera operation that changes camera number at the same site over time (e.g., failure of 1 of 2 cams)
source("sim2flank.y2D.HabMask.R")
source("init.2flank.y2D.HabMask.R")
source("LRmatch.R")
library(nimble)
source("NimbleModel 2Flank y2D HabMask.R")
source("NimbleFunctions 2Flank y2D HabMask.R")
source("sSamplerGrid.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

N <- 50
p0L <- 0.13
p0R <- 0.13
p0B <- 0.2
sigma <- 0.50
K <- 10

#making a grid that doesn't have a lower left origin at (0,0)
buff <- 2 #should be 3sigma or greater (at least for real data, unless using a habitat mask)
X <- expand.grid(13:18,13:18) #make some traps

##Make discrete habitat mask
res <- 1 #grid resolution
#state space. Must start at (0,0)
xlim.orig <- range(X[,1]) + c(-buff,buff)
ylim.orig <- range(X[,2]) + c(-buff,buff)
sub.x <- xlim.orig[1] #this is how far we have to shift the grid to start lower left corner at (0,0)
sub.y <- ylim.orig[1]
#new x and y limits with (0,0) origin
xlim <- xlim.orig - sub.x
ylim <- ylim.orig - sub.y
#xlim must fit n.cells.x grid cells exactly, same for ylim and n.cells.x
#if so, these will be integers
diff(xlim)/res
diff(ylim)/res

#unique x and y grid cell values.
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res)
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#store grid centroids in dSS, one row per grid cell
dSS <- as.matrix(expand.grid(x.vals,y.vals))
n.cells <- nrow(dSS)
#enumerate the grid cells
cells <- matrix(1:n.cells,nrow=n.cells.x,ncol=n.cells.y)


#OK, InSS=1 if grid cell is habitat, 0 otherwise
InSS <- rep(1,n.cells)
#Let's remove some arbitrary cells
#I'm going to remove the corner cells
InSS[dSS[,1]==x.vals[1]&dSS[,2]==y.vals[1]] <- 0
InSS[dSS[,1]==x.vals[n.cells.x]&dSS[,2]==y.vals[n.cells.y]] <- 0
InSS[dSS[,1]==x.vals[n.cells.x]&dSS[,2]==y.vals[1]] <- 0
InSS[dSS[,1]==x.vals[1]&dSS[,2]==y.vals[n.cells.y]] <- 0

image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y))

#finally, adjust the traps
X <- cbind(c(X[,1]-sub.x),c(X[,2]-sub.y))
points(X,pch=4)

#store grid objects here
grid.objects <- list(dSS=dSS,cells=cells,xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,
                     InSS=InSS,res=res)

J <- nrow(X)

#vector indicating which traps had 1 vs. 2 cameras. 
J.cams <- rep(1,J) #start with all 1
J.cams[seq(1,J,2)] <- 2 #set some to 2 (or not)
J.cams

#vector of trap operation for each trap, 1=operable, 2=inoberable. 
#Must use model version that uses 3D data to account for stations that change from single to double cams at the same site
#e.g., failure of 1 of the cameras.
K1D <- rep(K,J) #setting to perfect operation here

#You may know some flank linkages without both-side captures. If NA, we assume we do not.
# n.fixed <- 5 #minimum number of matched flanks for simulation
n.fixed <- NA #supply NA for only n.B flank matches to be known (if there are any).
#must supply this for simulation so we know to retain left and right known flank individuals with no captures
data <- sim2flank.y2D.HabMask(N=N,p0L=p0L,p0R=p0R,p0B=p0B,sigma=sigma,K=K,X=X,buff=buff,
                              J.cams=J.cams,K1D=K1D,n.fixed=n.fixed,
                          grid.objects=grid.objects)
n.fixed <- data$n.fixed #update in case n.B > n.fixed (n.B is number of both-side captured individuals)
#*Do not change* n.fixed from here down because all 0 left and/or right capture histories are retained
# for n.fixed individuals

#make sure simulated activity centers are inside the habitat mask
image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y))
points(data$s,col="lightblue")

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

M <- 100 #data augmentation level. If N posterior hits M, raise M and start over.

#initialize latent states, including flank matches.
#Can initialize to truth for simulated data for testing purposes.
#Otherwise, left and right flanks are initialized by 
#crudely minimizing the distance between the mean left and right
#side capture locations. 
#If n.fixed is NA, will fix flanks for n.B individuals if there are any, will not fix any otherwise
#plug in n.fixed here if different from n.B
#expecting grid.objects to be an element of data
nimbuild <- init.2flank.HabMask(data=data,M=M,n.fixed=n.fixed,initTrue=FALSE)

#crappy plot - ACs should not be in nonhabitat
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y))
points(X,pch=4,xlim=nimbuild$xlim,ylim=nimbuild$ylim)
points(nimbuild$s[nimbuild$z==1,],pch=16,col="lightblue")
y2D <- apply(nimbuild$y.true,c(1,2),sum)
for(i in 1:M){
  trapcaps <- which(y2D[i,]>0)
  for(j in 1:length(trapcaps)){
    lines(x=c(nimbuild$s[i,1],X[trapcaps[j],1]),y=c(nimbuild$s[i,2],X[trapcaps[j],2]))
  }
}
points(nimbuild$s[nimbuild$z==0,],pch=1,col="lightblue")

#inits for nimble
#use s inits in nimbuild! May not converge if you let nimble do it.
Niminits <- list(z=nimbuild$z,s=nimbuild$s,ID.L=nimbuild$ID.L,ID.R=nimbuild$ID.R,
                 y.L.true=nimbuild$y.true[,,2],
                 y.R.true=nimbuild$y.true[,,3],
                 sigma=sigma,p0S=0.13,psi=0.5,p0B=0.5) #let nimble initialize p0s, don't use true sigma as init, but ballpark start speeds convergence

#constants for Nimble
n.L <- nrow(data$y.L.obs)
n.R <- nrow(data$y.R.obs)
constants <- list(M=M,J=J,K1D=K1D,J.cams=J.cams,
                  xlim=nimbuild$xlim,ylim=nimbuild$ylim,n.L=n.L,n.R=n.R,res=res)

# Supply data to Nimble. Only individuals with deterministially-matched flanks have known z's.
if(n.fixed==0){
  z.data <- rep(NA,M)
}else{
  z.data <- c(rep(1,n.fixed),rep(NA,M-n.fixed))
}
#left and right flank data treated as fully latent
#Known flank matches will not be updated during MCMC.
Nimdata <- list(y.B.true=nimbuild$y.true[,,1],
                ID.L=rep(NA,n.L),ID.R=rep(NA,n.R),z=z.data,X=as.matrix(X),
                dummy.data=rep(1,M),cells=cells,InSS=InSS) #grid stuff

# set parameters to monitor
# parameters<-c('psi','p0B','p0L','p0R','sigma','N') #left and right side cap prob not shared
parameters <- c('psi','p0B','p0S','sigma','N')
parameters2 <- c("ID.L","ID.R") #monitor these with different thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)

#It's faster to not let nimble configure samplers for y.true. Only requesting nimble configure samplers
#that are not assigned below. If you add/remove parameters, make changes here.
config.nodes <- c("p0B","p0S","sigma","psi","z")
#thinning parameters 2 more to reduce posterior size
conf <- configureMCMC(Rmodel,monitors=parameters,monitors2=parameters2,
                      nodes=config.nodes,
                      thin=1,thin2=5,useConjugacy = TRUE)

###*required* sampler replacements
###Don't need to remove if not assigned by nimble
# conf$removeSampler("y.L.true")
# conf$removeSampler("y.R.true")
conf$addSampler(target = paste0("y.L.true[1:",M,",1:",J,"]"),
                type = 'IDLSampler',control = list(K1D=K1D,J.cams=J.cams,n.fixed=n.fixed,n.L=n.L,prop.scale=1),silent = TRUE)
conf$addSampler(target = paste0("y.R.true[1:",M,",1:",J,"]"),
                type = 'IDRSampler',control = list(K1D=K1D,J.cams=J.cams,n.fixed=n.fixed,n.R=n.R,prop.scale=1),silent = TRUE)

#required s sampler
# conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
pi.cell <- InSS/sum(InSS) #used to propose cells for z=0 inds. Uniformity across valid grid cells
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)
for(i in 1:M){
  calcNodes <- Rmodel$getDependencies(paste("s[",i,",1:2]"))
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSamplerGrid',control=list(i=i,pi.cell=pi.cell,xlim=xlim,ylim=ylim,
                                                     n.cells.x=n.cells.x,n.cells.y=n.cells.y,res=res),silent = TRUE)
}

#block these is posterior correlation is high. Often will be.
conf$removeSampler(c("p0S","sigma"))
conf$addSampler(target = c("p0S","sigma"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for easier debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model
start.time2 <- Sys.time()
#short run for demonstration
Cmcmc$run(5000,reset=FALSE) #can continue MCMC by rerunning this line with reset=FALSE
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

#extract and plot posterior samples
library(coda)
mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[500:nrow(mvSamples),]))


# can see final activity centers are still in valid habitat
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(InSS,length(x.vals),length(y.vals)))
points(Cmodel$s,pch=16,col="lightblue")

#posterior sample match probs. not removing burnin here
mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
postL <- mvSamples2[,1:n.L] #which 1:M individual is left flank matched to on each iteration>
postR <- mvSamples2[,(n.L+1):(n.R+n.L)] #same for right flanks

#Calculate posterior prob that left flank l matches right flank r
postprobs <- matrix(NA,nrow=n.L,ncol=n.R)
for(l in 1:n.L){
  for(r in 1:n.R){
    tmp <- postL[,l]==postR[,r]
    postprobs[l,r] <- mean(tmp[-1])
  }
}

#Posterior match probability that left flanks are matched to correct right flanks (for simulated data)
#Higher match probs indicate more certainty in matches and less precision lost relative to true data
#**match probs should be 1 for 1:n.fixed individuals. If not, something went wrong.**
for(i in 1:n.L){
  idx <- which(data$ID.R==data$ID.L[i])
  if(length(idx)>0){
    print(postprobs[i,idx])
  }else{
    print(NA) #was no matching right flank
  }
}
#Posterior match probability that right flanks are matched to correct left flanks (for simulated data)
for(i in 1:n.R){
  idx <- which(data$ID.L==data$ID.R[i])
  if(length(idx)>0){
    print(postprobs[idx,i])
  }else{
    print(NA) #was no matching left flank
  }
}

#posterior for number of matching flank pairs
n.match <- rep(NA,nrow(postL)-1)
for(i in 2:nrow(postL)){
  n.match[i] <- sum(postL[i,]%in%postR[i,])
}
plot(mcmc(n.match[-1]))
sum(data$ID.L%in%data$ID.R) #truth

