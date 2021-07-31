e2dist<-function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.2flank=function(data=data,M=M,initTrue=FALSE){
  J=dim(data$y.L.obs)[2]
  K=dim(data$y.L.obs)[3]
  if(initTrue){
    if(!all(c("y.B","y.L","y.R")%in%names(data)))stop("To init to truth, true data y.B, y.L, and y.R must be in data object.")
        #trying to get these in the right order...
    y.B=data$y.B.obs
    y.L=data$y.L.obs
    y.R=data$y.R.obs
    ID.B.true=data$ID.B
    ID.L.true=data$ID.L
    ID.R.true=data$ID.R
    n.B=nrow(y.B)
    n.L=nrow(y.L)
    n.R=nrow(y.R)
    
    y.true=array(0,dim=c(M,J,K,3))
    ID.L=rep(NA,n.L)
    ID.R=rep(NA,n.R)
    if(n.B>0){
      for(i in 1:n.B){
        y.true[i,,,1]=y.B[i,,]
        y.true[i,,,2]=y.L[i,,]
        y.true[i,,,3]=y.R[i,,]
        ID.L[i]=ID.R[i]=i
      }
    }
    if(n.B>0){
      ID.remain=unique(c(ID.L.true[-c(1:n.B)],ID.R.true[-c(1:n.B)]))
    }else{
      ID.remain=unique(c(ID.L.true,ID.R.true))
    }
    for(i in 1:length(ID.remain)){
      idx=which(ID.L.true==ID.remain[i])
      if(length(idx)>0){
        y.true[i+n.B,,,2]=y.L[idx,,]
        ID.L[idx]=i+n.B
        
      }
      idx=which(ID.R.true==ID.remain[i])
      if(length(idx)>0){
        y.true[i+n.B,,,3]=y.R[idx,,]
        ID.R[idx]=i+n.B
      }
    }
    #used to test below
    n.B=sum(rowSums(y.B)>0)
    n.L=sum(rowSums(y.L)>0)
    n.R=sum(rowSums(y.R)>0)
  }else{
    y.true=array(0,dim=c(M,J,K,3))
    #observed data
    y.B=data$y.B.obs
    y.L=data$y.L.obs
    y.R=data$y.R.obs
    #caps per type
    n.B=nrow(y.B)
    n.L=nrow(y.L)
    n.R=nrow(y.R)
    #initialize lefts and rights to "both side" order
    if(n.B>0){
      y.true[1:n.B,,,1]=y.B
    }
   
    flank.init=LRmatch(M=M, left=apply(y.L,c(1,2),sum), nleft=n.L-n.B, right=apply(y.R,c(1,2),sum),nright=n.R-n.B, X=X, Nfixed=n.B)
    
    # M=M
    # left=apply(y.L,c(1,2),sum)
    # nleft=n.L-n.B
    # right=apply(y.R,c(1,2),sum)
    # nright=n.R-n.B
    # X=X
    # Nfixed=n.B
    
    ID.L=flank.init$ID.L
    ID.R=flank.init$ID.R
    for(i in 1:length(ID.L)){
      y.true[ID.L[i],,,2]=y.L[i,,]
    }
    for(i in 1:length(ID.R)){
      y.true[ID.R[i],,,3]=y.R[i,,]
    }
  }
  #test
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
  X=as.matrix(data$X)
  xlim=range(X[,1])+c(-data$buff,data$buff)
  ylim=range(X[,2])+c(-data$buff,data$buff)
  #initialize s
  s=cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
  y.true2D=apply(y.true,c(1,2),sum)
  for(i in 1:M){
    idx=which(y.true2D[i,]>0)
    if(length(idx)>0){
      if(length(idx)>1){
        s[i,]=colMeans(X[idx,])
      }else{
        s[i,]=X[idx,]
      }
    }
  }
  return(list(y.true=y.true,ID.L=ID.L,ID.R=ID.R,xlim=xlim,ylim=ylim,s=s))
}