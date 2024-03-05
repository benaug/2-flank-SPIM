e2dist<-function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

getCellR = function(u,res,cells,xlim,ylim){
  inout=1*(u[1]>xlim[1]&u[1]<xlim[2]&u[2]>ylim[1]&u[2]<ylim[2])
  if(inout==1){
    this.cell=cells[trunc(u[1]/res)+1,trunc(u[2]/res)+1]
  }else{
    this.cell=0
  }
  return(this.cell)
}

init.2flank.HabMask <- function(data=data,M=M,n.fixed=NA,initTrue=FALSE){
  J <- dim(data$y.L.obs)[2]
  K <- dim(data$y.L.obs)[3]
  X <- data$X
  n.B <- nrow(data$y.B.obs)
  if(is.na(n.fixed)|(n.B==n.fixed)){
    if(n.B>0){
      print("Assuming the only known flank matches are for the n.B both side captured individuals.")
      n.fixed <- n.B
    }else{
      print("Assuming no known flank linkages because no both side captures and n.fixed not supplied or is 0.")
      n.fixed <- 0
    }
  }else{
    if(n.fixed<n.B)stop("n.fixed should be greater than n.B. We should know flank matches for n.B individuals.")
    n.diff <- n.fixed-n.B
    print("Assuming 1:n.fixed individual flank matches are known instead of 1:n.B because n.fixed supplied and is larger than n.B.")
  }
  
  if(initTrue){#initialize flanks to true matches?
    if(!all(c("y.B","y.L","y.R")%in%names(data)))stop("To init to truth, true data y.B, y.L, and y.R must be in data object.")
    y.B <- data$y.B.obs
    y.L <- data$y.L.obs
    y.R <- data$y.R.obs
    ID.B.true <- data$ID.B
    ID.L.true <- data$ID.L
    ID.R.true <- data$ID.R
    n.B <- nrow(y.B)
    n.L <- nrow(y.L)
    n.R <- nrow(y.R)
    
    y.true <- array(0,dim=c(M,J,K,3))
    ID.L <- rep(NA,n.L)
    ID.R <- rep(NA,n.R)
    if(n.fixed>0){
      if(n.B>0){
        y.true[1:n.B,,,1] <- y.B[1:n.B,,]
      }
      y.true[1:n.fixed,,,2] <- y.L[1:n.fixed,,]
      y.true[1:n.fixed,,,3] <- y.R[1:n.fixed,,]
      ID.L[1:n.fixed] <- ID.R[1:n.fixed] <- 1:n.fixed
    }
    if(n.fixed>0){
      ID.remain <- unique(c(ID.L.true[-c(1:n.fixed)],ID.R.true[-c(1:n.fixed)]))
    }else{
      ID.remain <- unique(c(ID.L.true,ID.R.true))
    }
    for(i in 1:length(ID.remain)){
      idx <- which(ID.L.true==ID.remain[i])
      if(length(idx)>0){
        y.true[i+n.fixed,,,2] <- y.L[idx,,]
        ID.L[idx] <- i+n.fixed
        
      }
      idx <- which(ID.R.true==ID.remain[i])
      if(length(idx)>0){
        y.true[i+n.fixed,,,3] <- y.R[idx,,]
        ID.R[idx] <- i+n.fixed
      }
    }
    #used to test below
    n.B <- sum(rowSums(y.B)>0)
    n.L <- sum(rowSums(y.L)>0)
    n.R <- sum(rowSums(y.R)>0)
  }else{ #initialize without knowing truth
    y.true <- array(0,dim=c(M,J,K,3))
    #observed data, not using known data to initialize
    y.B <- data$y.B.obs
    y.L <- data$y.L.obs
    y.R <- data$y.R.obs
    #caps per type
    n.B <- nrow(y.B)
    n.L <- nrow(y.L)
    n.R <- nrow(y.R)
    #initialize lefts and rights to "both side" order
    if(n.B>0){
      y.true[1:n.B,,,1] <- y.B
    }

    #initialized unknown flank matches using spatial locations of captures
    flank.init <- LRmatch(M=M, left=apply(y.L,c(1,2),sum), nleft=n.L-n.fixed,
                         right=apply(y.R,c(1,2),sum),nright=n.R-n.fixed, X=X, Nfixed=n.fixed)
    ID.L <- flank.init$ID.L
    ID.R <- flank.init$ID.R
    
    for(i in 1:length(ID.L)){
      y.true[ID.L[i],,,2] <- y.L[i,,]
    }
    for(i in 1:length(ID.R)){
      y.true[ID.R[i],,,3] <- y.R[i,,]
    }
  }
  #for simulated data, test to make sure initialization algorithm is correct
  if(all(c("y.B","y.L","y.R")%in%names(data))){
    if(n.B>0){
      for(i in 1:n.B){
        if(!all(y.true[i,,,1]==y.B[i,,]))stop("Error initializing y.B")
      }
    }
    for(i in 1:n.L){
      if(!all(y.true[ID.L[i],,,2]==y.L[i,,]))stop("Error initializing y.L")
    }
    for(i in 1:n.R){
      if(!all(y.true[ID.R[i],,,3]==y.R[i,,]))stop("Error initializing y.R")
    }
  }
  
  #initialize s inside habitat mask
  list2env(data$grid.objects, envir = .GlobalEnv)
  X <- as.matrix(data$X)
  s <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
  y.true2D <- apply(y.true,c(1,2),sum)
  for(i in 1:M){
    idx <- which(y.true2D[i,]>0)
    if(length(idx)>0){
      if(length(idx)>1){
        s[i,] <- colMeans(X[idx,])
      }else{
        s[i,] <- X[idx,]
      }
    }
  }
  #move any initialized outside state space
  for(i in 1:M){
    s.cell <- getCellR(s[i,],res,cells,xlim,ylim)
    if(InSS[s.cell]==0){#not in SS, move to nearest cell
      dists <- sqrt((dSS[s.cell,1]-dSS[,1])^2+(dSS[s.cell,2]-dSS[,2])^2)
      dists[InSS==0] <- Inf
      pick <- which(dists==min(dists))[1] #if more than 1, just use first
      s[i,] <- dSS[pick,]
    }
  }
  #initialize z to only captured individuals. Can let nimble initialize z if you want.
  z <- 1*(rowSums(y.true)>0)
  return(list(y.true=y.true,ID.L=ID.L,ID.R=ID.R,xlim=xlim,ylim=ylim,s=s,z=z))
}