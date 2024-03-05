e2dist<-function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim2flank.y2D.HabMask<-
  function(N=NA,p0L=NA,p0R=NA,p0B=NA,sigma=NA,K=NA,X=X,buff=NA,J.cams=NA,K1D=NA,n.fixed=NA,
           grid.objects=grid.objects){
    J <- nrow(X)
    # simulate a population of activity centers
    s <- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    
    # simulate a population of activity centers
    #chose cells first at random
    s.cell <- sample(which(grid.objects$InSS==1),N,replace=TRUE)
    #assign spatial uniform location inside selected cell
    s <- matrix(NA,N,2)
    s.xlim <- s.ylim <- matrix(NA,N,2)
    for(i in 1:N){
      s.xlim[i,1] <- grid.objects$dSS[s.cell[i],1]-res/2
      s.xlim[i,2] <- grid.objects$dSS[s.cell[i],1]+res/2
      s.ylim[i,1] <- grid.objects$dSS[s.cell[i],2]-res/2
      s.ylim[i,2] <- grid.objects$dSS[s.cell[i],2]+res/2
      s[i,] <- cbind(runif(1,s.xlim[i,1],s.xlim[i,2]),
                     runif(1,s.ylim[i,1],s.ylim[i,2]))
    }
    
    #capture individuals
    D <- e2dist(s,X)
    expstuff <- exp(-D*D/(2*sigma*sigma))
    pdL <- p0L*expstuff
    pdR <- p0R*expstuff
    pdB <- p0B*expstuff
    y.L <- y.R <- y.B <- matrix(0,N,J)
    for(i in 1:N){
      for(j in 1:J){
        if(J.cams[j]==1){
          y.L[i,j] <- rbinom(1,K1D[j],pdL[i,j])
          y.R[i,j] <- rbinom(1,K1D[j],pdR[i,j])
        }else{ #2 cams
          #P(A or B)=P(A)+P(B)-P(A and B)
          y.L[i,j] <- rbinom(1,K1D[j],2*pdL[i,j]-pdL[i,j]^2) #single side p. two chances for capture with 2 cameras
          y.R[i,j] <- rbinom(1,K1D[j],2*pdR[i,j]-pdR[i,j]^2)
          y.B[i,j] <- rbinom(1,K1D[j],pdB[i,j])
        }
      }
    }
    
    ######Process data#############
    ID.B <- which(apply(y.B,1,sum)>0)
    n.B <- length(ID.B)
    ID.L <- which(rowSums(y.L)>0)
    ID.R <- which(rowSums(y.R)>0)
    ID.L <- setdiff(ID.L,ID.B)
    ID.R <- setdiff(ID.R,ID.B)
    ID.L <- c(ID.B,ID.L)
    ID.R <- c(ID.B,ID.R)
    if(!is.na(n.fixed)){
      if(n.fixed<n.B){
        print("More both-side captured individuals than n.fixed. Changing n.fixed to n.B")
        n.fixed <- n.B
      }
      #add uncaptured but known left and right flanks
      ID.L <- sort(unique(c(ID.L,1:n.fixed)))
      ID.R <- sort(unique(c(ID.R,1:n.fixed)))
    }
    
    
    Nknown <- length(ID.B)
    n <- sum(apply(y.B+y.L+y.R,1,sum)>0)
    
    #remove uncaptured individuals, put ID.B at top of y.L and y.R
    y.B.obs <- y.B[ID.B,]
    y.L.obs <- y.L[ID.L,]
    y.R.obs <- y.R[ID.R,]
    n.B <- length(ID.B)
    n.L <- length(ID.L)
    n.R <- length(ID.R)
    if(n.B==1){
      y.B.obs <- matrix(y.B.obs,1,J)
    }
    if(n.L==1){
      y.L.obs <- matrix(y.L.obs,1,J)
    }
    if(n.R==1){
      y.R.obs <- matrix(y.R.obs,1,J)
    }
    if(is.na(n.fixed)){
      n.fixed <- n.B
    }
    out <-list(y.B.obs=y.B.obs,y.L.obs=y.L.obs,y.R.obs=y.R.obs,X=X,K=K,buff=buff,
              y.B=y.B,y.L=y.L,y.R=y.R,s=s,s.cell=s.cell,n=n,ID.B=ID.B,ID.L=ID.L,ID.R=ID.R,n.fixed=n.fixed,
              grid.objects=grid.objects)

    return(out)
  }