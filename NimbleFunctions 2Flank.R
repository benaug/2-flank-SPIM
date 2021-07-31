IDdummyfun <- nimbleFunction(
  run = function(ID.L=double(1)){ 
    returnType(double(0))
    return(0)
  }
)
GetD2 <- nimbleFunction(
  run = function(s=double(1),X=double(2), z=double(0)){ 
    returnType(double(1))
    J=nimDim(X)[1]
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      return(d2)
    }
  }
)
GetDetectionProb <- nimbleFunction(
  run = function(d2=double(1), p0=double(0),sigma=double(0), z=double(0)){ 
    returnType(double(1))
    J=nimDim(d2)[1]
    pd=rep(0,J)
    if(z==0) return(pd)
    if(z==1){
      pd <-p0*exp(-d2/(2*sigma^2))
      return(pd)
    }
  }
)
dBernoulliVectorBoth <- nimbleFunction(
  run = function(x = double(2), pd=double(1), K2D=double(2), z = double(0),log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      J=nimDim(K2D)[1]
      K=nimDim(K2D)[2]
      logProb=0
      for(j in 1:J){
        for(k in 1:K){
          if(K2D[j,k]==2){
            logProb=logProb+dbinom(x[j,k],prob=pd[j],size=1,log=TRUE)
          }
        }
      }
      return(logProb)
    }
  }
)
dBernoulliVectorSingle <- nimbleFunction( #don't need to use ID, but nimble must think it is necessary somewhere in model structure
  run = function(x = double(2), pd=double(1), K2D=double(2), z = double(0),log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      J=nimDim(K2D)[1]
      K=nimDim(K2D)[2]
      logProb=0
      for(j in 1:J){
        for(k in 1:K){
          if(K2D[j,k]==1){
            logProb=logProb+dbinom(x[j,k],prob=pd[j],size=1,log=TRUE)
          }else if(K2D[j,k]==2){
            logProb=logProb+dbinom(x[j,k],prob=2*pd[j]-pd[j]^2,size=1,log=TRUE)
          }
        }
      }
      return(logProb)
    }
  }
)
#dummy function
rBernoulliVectorBoth <- nimbleFunction(
  run = function(n=integer(0),pd=double(1), K2D=double(2), z = double(0)) {
    returnType(double(2))
    J=nimDim(K2D)[1]
    K=nimDim(K2D)[2]
    y.true=matrix(0,J,K)
    return(y.true)
  }
)
#dummy function
rBernoulliVectorSingle <- nimbleFunction(
  run = function(n=integer(0),pd=double(1), K2D=double(2), z = double(0)) {
    returnType(double(2))
    J=nimDim(K2D)[1]
    K=nimDim(K2D)[2]
    y.true=matrix(0,J,K)
    return(y.true)
  }
)

IDLSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K2D <- control$K2D
    n.L <- control$n.L
    n.fixed <- control$n.fixed
    prop.scale <- control$prop.scale
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    if(n.L>0){ #skip if no lefts
      z <- model$z
      s <- model$s
      y.L.true <- model$y.L.true
      ID.L<- model$ID.L
      pd.L <- model$pd.L
      M <- nimDim(y.L.true)[1]
      J <- nimDim(y.L.true)[2]
      K <- nimDim(y.L.true)[3]

      #precalculate log likelihoods. Can pull from nimble, but structure depends on nimble model structure (vectorized or not)
      ll.y.L <- array(0,dim=c(M,J,K))
      for(i in 1:M){
        if(z[i]==1){
          for(j in 1:J){
            for(k in 1:K){
              if(K2D[j,k]==1){
                ll.y.L[i,j,k]=dbinom(y.L.true[i,j,k],prob=pd.L[i,j],size=1,log=TRUE)
              }else if(K2D[j,k]==2){
                ll.y.L[i,j,k]=dbinom(y.L.true[i,j,k],prob=2*pd.L[i,j]-pd.L[i,j]^2,size=1,log=TRUE)
              }
            }
          }
        }
      }
      ll.y.L.cand <- ll.y.L
      for(l in (n.fixed+1):n.L){
        y.L.true.cand <- y.L.true
        this.i=ID.L[l]
        #select an individual to swap it to
        d2=(s[this.i,1]-s[,1])^2+(s[this.i,2]-s[,2])^2
        propprobs=exp(-d2/(2*(prop.scale*model$sigma)^2)) #distance-based proposal
        propprobs[this.i]=0 #don't choose focal
        propprobs[z==0]=0 #must select an individual with z==1
        propprobs[1:n.fixed]=0 #cannot select an individual with fixed flanks
        propprobs=propprobs/sum(propprobs)
        cand.i=rcat(1,propprobs)
        if(this.i!=cand.i){ #trying this, can choose yourself
          cand.ID.idx=which(ID.L==cand.i)#which ID.L index is this individual. Will be empty if not in ID.L
          y.L.true.cand[this.i,,] <- y.L.true[cand.i,,]
          y.L.true.cand[cand.i,,] <- y.L.true[this.i,,]
          #update ll.y.L
          for(j in 1:J){
            for(k in 1:K){
              if(K2D[j,k]==1){
                ll.y.L.cand[this.i,j,k]=dbinom(y.L.true.cand[this.i,j,k],prob=pd.L[this.i,j],size=1,log=TRUE)
                ll.y.L.cand[cand.i,j,k]=dbinom(y.L.true.cand[cand.i,j,k],prob=pd.L[cand.i,j],size=1,log=TRUE)
              }else if(K2D[j,k]==2){
                ll.y.L.cand[this.i,j,k]=dbinom(y.L.true.cand[this.i,j,k],prob=2*pd.L[this.i,j]-pd.L[this.i,j]^2,size=1,log=TRUE)
                ll.y.L.cand[cand.i,j,k]=dbinom(y.L.true.cand[cand.i,j,k],prob=2*pd.L[cand.i,j]-pd.L[cand.i,j]^2,size=1,log=TRUE)
              }
            }
          }
          #calculate backwards proposal probs
          d2=(s[cand.i,1]-s[,1])^2+(s[cand.i,2]-s[,2])^2
          backprobs=exp(-d2/(2*(prop.scale*model$sigma)^2)) #distance-based proposal
          backprobs[cand.i]=0 #don't choose focal
          backprobs[z==0]=0 #must select an individual with z==1
          backprobs[1:n.fixed]=0 #cannot select an individual with fixed flanks
          backprobs=backprobs/sum(backprobs)
          lp_initial <- sum(ll.y.L[this.i,,])+sum(ll.y.L[cand.i,,])
          lp_proposed <- sum(ll.y.L.cand[this.i,,])+sum(ll.y.L.cand[cand.i,,])
          log_MH_ratio <- (lp_proposed+log(backprobs[this.i])) - (lp_initial+log(propprobs[cand.i]))
          accept <- decide(log_MH_ratio)
          if(accept){
            y.L.true[this.i,,]=y.L.true.cand[this.i,,]
            y.L.true[cand.i,,]=y.L.true.cand[cand.i,,]
            ll.y.L[this.i,,]=ll.y.L.cand[this.i,,]
            ll.y.L[cand.i,,]=ll.y.L.cand[cand.i,,]
            #swap flank indices
            ID.L[l]=cand.i
            #if cand.i was in ID.L, update this ID index
            tmp=nimDim(cand.ID.idx)[1]
            if(tmp>0){
              ID.L[cand.ID.idx]=this.i
            }
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