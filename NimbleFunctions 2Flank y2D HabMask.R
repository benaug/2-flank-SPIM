dCell <- nimbleFunction(
  run = function(x = double(0), InSS = integer(0),log = integer(0)) {
    returnType(double(0))
    if(InSS==1){
      logProb <- 0
    }else{
      logProb <- -Inf
    }
    return(logProb)
  }
)
#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),InSS = integer(0)) {
    returnType(double(0))
    return(0)
  }
)

IDdummyfun <- nimbleFunction(
  run = function(ID.L=double(1),ID.R=double(1)){ 
    returnType(double(0))
    return(0)
  }
)
GetD2 <- nimbleFunction(
  run = function(s=double(1),X=double(2), z=double(0)){ 
    returnType(double(1))
    J <- nimDim(X)[1]
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      return(d2)
    }
  }
)
#make pd 3D here and in PMF functions if adding occasion effects to model. 
#Must modify custom ID updates to reflect this change.
GetDetectionProb <- nimbleFunction(
  run = function(d2=double(1), p0=double(0),sigma=double(0), z=double(0)){ 
    returnType(double(1))
    J <- nimDim(d2)[1]
    pd <- rep(0,J)
    if(z==0) return(pd)
    if(z==1){
      pd <- p0*exp(-d2/(2*sigma^2))
      return(pd)
    }
  }
)
dBernoulliVectorBoth <- nimbleFunction(
  run = function(x = double(1), pd=double(1), K1D=double(1), J.cams=double(1), z = double(0),log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      J <- nimDim(K1D)[1]
      logProb <- 0
      for(j in 1:J){
        if(J.cams[j]==2){
          logProb <- logProb + dbinom(x[j],prob=pd[j],size=K1D[j],log=TRUE)
        }
      }
      return(logProb)
    }
  }
)
dBernoulliVectorSingle <- nimbleFunction(
  run = function(x = double(1), pd=double(1), K1D=double(1), J.cams=double(1), z = double(0),log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      J <- nimDim(K1D)[1]
      logProb <- 0
      for(j in 1:J){
        if(J.cams[j]==1){
          logProb <- logProb + dbinom(x[j],prob=pd[j],size=K1D[j],log=TRUE)
        }else{
          logProb <- logProb + dbinom(x[j],prob=2*pd[j]-pd[j]^2,size=K1D[j],log=TRUE)
        }
      }
      return(logProb)
    }
  }
)
#dummy function to make this work in parallel
rBernoulliVectorBoth <- nimbleFunction(
  run = function(n=integer(0),pd=double(1), K1D=double(1), J.cams=double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(K1D)[1]
    y.true <- rep(0,J)
    return(y.true)
  }
)
#dummy function to make this work in parallel
rBernoulliVectorSingle <- nimbleFunction(
  run = function(n=integer(0),pd=double(1), K1D=double(1), J.cams=double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(K1D)[1]
    y.true <- rep(0,J)
    return(y.true)
  }
)

#Left flank update
IDLSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K1D <- control$K1D
    J.cams <- control$J.cams
    n.L <- control$n.L
    n.fixed <- control$n.fixed
    prop.scale <- control$prop.scale #can scale the distance proposal, but a value of 1 should be a good choice.
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    if(n.L>0){ #skip if no lefts
      z <- model$z
      s <- model$s
      y.L.true <- model$y.L.true
      ID.L <- model$ID.L
      pd.L <- model$pd.L
      sigma <- model$sigma
      M <- nimDim(y.L.true)[1]
      J <- nimDim(y.L.true)[2]
      
      #precalculate log likelihoods. Can pull from nimble, but sticking with this 
      ll.y.L <- array(0,dim=c(M,J))
      for(i in 1:M){
        if(z[i]==1){
          for(j in 1:J){
            if(J.cams[j]==1){
              ll.y.L[i,j] <- dbinom(y.L.true[i,j],prob=pd.L[i,j],size=K1D[j],log=TRUE)
            }else{ #2 cams
              ll.y.L[i,j] <- dbinom(y.L.true[i,j],prob=2*pd.L[i,j]-pd.L[i,j]^2,size=K1D[j],log=TRUE)
            }
          }
        }
      }
      ll.y.L.cand <- ll.y.L
      for(l in (n.fixed+1):n.L){
        y.L.true.cand <- y.L.true
        this.i <- ID.L[l]
        #select an individual to swap it to
        d2 <- (s[this.i,1]-s[,1])^2+(s[this.i,2]-s[,2])^2
        propprobs <- exp(-d2/(2*(prop.scale*sigma[1])^2)) #distance-based proposal. must reference index for sigma
        propprobs[this.i] <- 0 #don't choose focal
        propprobs[z==0] <- 0 #must select an individual with z==1
        if(n.fixed>0){
          propprobs[1:n.fixed] <- 0 #cannot select an individual with fixed flanks
        }
        propprobs <- propprobs/sum(propprobs)
        cand.i <- rcat(1,prob=propprobs)
        
        cand.ID.idx <- which(ID.L==cand.i)#which ID.L index is this individual. Will be empty if not in ID.L
        y.L.true.cand[this.i,] <- y.L.true[cand.i,]
        y.L.true.cand[cand.i,] <- y.L.true[this.i,]
        #update ll.y.L
        for(j in 1:J){
          if(J.cams[j]==1){
            ll.y.L.cand[this.i,j] <- dbinom(y.L.true.cand[this.i,j],prob=pd.L[this.i,j],size=K1D[j],log=TRUE)
            ll.y.L.cand[cand.i,j] <- dbinom(y.L.true.cand[cand.i,j],prob=pd.L[cand.i,j],size=K1D[j],log=TRUE)
          }else{ #2 cams
            ll.y.L.cand[this.i,j] <- dbinom(y.L.true.cand[this.i,j],prob=2*pd.L[this.i,j]-pd.L[this.i,j]^2,size=K1D[j],log=TRUE)
            ll.y.L.cand[cand.i,j] <- dbinom(y.L.true.cand[cand.i,j],prob=2*pd.L[cand.i,j]-pd.L[cand.i,j]^2,size=K1D[j],log=TRUE)
          }
        }
        #calculate backwards proposal probs
        d2 <- (s[cand.i,1]-s[,1])^2+(s[cand.i,2]-s[,2])^2
        backprobs <- exp(-d2/(2*(prop.scale*sigma[1])^2)) #distance-based proposal
        backprobs[cand.i] <- 0 #don't choose focal
        backprobs[z==0] <- 0 #must select an individual with z==1
        if(n.fixed>0){
          backprobs[1:n.fixed] <- 0 #cannot select an individual with fixed flanks
        }
        backprobs <- backprobs/sum(backprobs)
        lp_initial <- sum(ll.y.L[this.i,])+sum(ll.y.L[cand.i,])
        lp_proposed <- sum(ll.y.L.cand[this.i,])+sum(ll.y.L.cand[cand.i,])
        log_MH_ratio <- (lp_proposed+log(backprobs[this.i])) - (lp_initial+log(propprobs[cand.i]))
        accept <- decide(log_MH_ratio)
        if(accept){
          y.L.true[this.i,] <- y.L.true.cand[this.i,]
          y.L.true[cand.i,] <- y.L.true.cand[cand.i,]
          ll.y.L[this.i,] <- ll.y.L.cand[this.i,]
          ll.y.L[cand.i,] <- ll.y.L.cand[cand.i,]
          #swap flank indices
          ID.L[l] <- cand.i
          #if cand.i was in ID.L, update this ID index
          tmp <- nimDim(cand.ID.idx)[1]
          if(tmp>0){
            ID.L[cand.ID.idx] <- this.i
          }
        }
      }

      #put everything back into model$stuff
      model$y.L.true <<- y.L.true
      model$ID.L <<- ID.L
      model$calculate(calcNodes) #update logprob
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)

#Right flank update
IDRSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K1D <- control$K1D
    J.cams <- control$J.cams
    n.R <- control$n.R
    n.fixed <- control$n.fixed
    prop.scale <- control$prop.scale
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    if(n.R>0){ #skip if no lefts
      z <- model$z
      s <- model$s
      y.R.true <- model$y.R.true
      ID.R<- model$ID.R
      pd.R <- model$pd.R
      M <- nimDim(y.R.true)[1]
      J <- nimDim(y.R.true)[2]
      sigma <- model$sigma
      #precalculate log likelihoods. Can pull from nimble, but sticking with this.
      ll.y.R <- array(0,dim=c(M,J))
      for(i in 1:M){
        if(z[i]==1){
          for(j in 1:J){
            if(J.cams[j]==1){
              ll.y.R[i,j] <- dbinom(y.R.true[i,j],prob=pd.R[i,j],size=K1D[j],log=TRUE)
            }else{
              ll.y.R[i,j] <- dbinom(y.R.true[i,j],prob=2*pd.R[i,j]-pd.R[i,j]^2,size=K1D[j],log=TRUE)
            }
          }
        }
      }
      ll.y.R.cand <- ll.y.R
      for(l in (n.fixed+1):n.R){
        y.R.true.cand <- y.R.true
        this.i <- ID.R[l]
        #select an individual to swap it to
        d2 <- (s[this.i,1]-s[,1])^2+(s[this.i,2]-s[,2])^2
        propprobs <- exp(-d2/(2*(prop.scale*sigma[1])^2)) #distance-based proposal
        propprobs[this.i] <- 0 #don't choose focal
        propprobs[z==0] <- 0 #must select an individual with z==1
        if(n.fixed>0){
          propprobs[1:n.fixed] <- 0 #cannot select an individual with fixed flanks
        }
        propprobs <- propprobs/sum(propprobs)
        cand.i <- rcat(1,propprobs)
        cand.ID.idx <- which(ID.R==cand.i)#which ID.R index is this individual. Will be empty if not in ID.R
        y.R.true.cand[this.i,] <- y.R.true[cand.i,]
        y.R.true.cand[cand.i,] <- y.R.true[this.i,]
        #update ll.y.R
        for(j in 1:J){
          if(J.cams[j]==1){
            ll.y.R.cand[this.i,j] <- dbinom(y.R.true.cand[this.i,j],prob=pd.R[this.i,j],size=K1D[j],log=TRUE)
            ll.y.R.cand[cand.i,j] <- dbinom(y.R.true.cand[cand.i,j],prob=pd.R[cand.i,j],size=K1D[j],log=TRUE)
          }else{
            ll.y.R.cand[this.i,j] <- dbinom(y.R.true.cand[this.i,j],prob=2*pd.R[this.i,j]-pd.R[this.i,j]^2,size=K1D[j],log=TRUE)
            ll.y.R.cand[cand.i,j] <- dbinom(y.R.true.cand[cand.i,j],prob=2*pd.R[cand.i,j]-pd.R[cand.i,j]^2,size=K1D[j],log=TRUE)
          }
        }
        #calculate backwards proposal probs
        d2 <- (s[cand.i,1]-s[,1])^2+(s[cand.i,2]-s[,2])^2
        backprobs <- exp(-d2/(2*(prop.scale*sigma[1])^2)) #distance-based proposal
        backprobs[cand.i] <- 0 #don't choose focal
        backprobs[z==0] <- 0 #must select an individual with z==1
        if(n.fixed>0){
          backprobs[1:n.fixed] <- 0 #cannot select an individual with fixed flanks
        }
        backprobs <- backprobs/sum(backprobs)
        lp_initial <- sum(ll.y.R[this.i,])+sum(ll.y.R[cand.i,])
        lp_proposed <- sum(ll.y.R.cand[this.i,])+sum(ll.y.R.cand[cand.i,])
        log_MH_ratio <- (lp_proposed+log(backprobs[this.i])) - (lp_initial+log(propprobs[cand.i]))
        accept <- decide(log_MH_ratio)
        
        if(accept){
          y.R.true[this.i,] <- y.R.true.cand[this.i,]
          y.R.true[cand.i,] <- y.R.true.cand[cand.i,]
          ll.y.R[this.i,] <- ll.y.R.cand[this.i,]
          ll.y.R[cand.i,] <- ll.y.R.cand[cand.i,]
          #swap flank indices
          ID.R[l] <- cand.i
          #if cand.i was in ID.R, update this ID index
          tmp <- nimDim(cand.ID.idx)[1]
          if(tmp>0){
            ID.R[cand.ID.idx] <- this.i
          }
        }
      }
      #put everything back into model$stuff
      model$y.R.true <<- y.R.true
      model$ID.R <<- ID.R
      model$calculate(calcNodes) #update logprob
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)