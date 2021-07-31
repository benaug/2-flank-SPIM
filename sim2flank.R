e2dist<-function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim2flank<-
  function(N=NA,p0L=NA,p0R=NA,p0B=NA,sigma=NA,K=NA,X=X,buff=NA,K2D=NA){
    library(abind)
    J<- nrow(X)
    # simulate a population of activity centers
    s=cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    #capture individuals
    D=e2dist(s,X)
    expstuff=exp(-D*D/(2*sigma*sigma))
    pdL=p0L*expstuff
    pdR=p0R*expstuff
    pdB=p0B*expstuff
    y.L <-y.R <- y.B <-array(0,dim=c(N,J,K))
    for(i in 1:N){
      for(j in 1:J){
        for(k in 1:K){
          #P(A or B)=P(A)+P(B)-P(A and B)
          y.L[i,j,k]=rbinom(1,1,(K2D[j,k]==1)*pdL[i,j]+(K2D[j,k]==2)*(2*pdL[i,j]-pdL[i,j]^2)) #single side p. two chances for capture with 2 cameras
          y.R[i,j,k]=rbinom(1,1,(K2D[j,k]==1)*pdR[i,j]+(K2D[j,k]==2)*(2*pdR[i,j]-pdR[i,j]^2)) #single side p
          y.B[i,j,k]=rbinom(1,1,pdB[i,j])*(K2D[j,k]==2)*1  #y.B side lambda multiplied by indicator for 2 traps at site-occasion
        }
      }
    }
    
    ######Process data#############
    ID.B=which(apply(y.B,1,sum)>0)
    ID.L=which(rowSums(y.L)>0)
    ID.R=which(rowSums(y.R)>0)
    ID.L=setdiff(ID.L,ID.B)
    ID.R=setdiff(ID.R,ID.B)
    ID.L=c(ID.B,ID.L)
    ID.R=c(ID.B,ID.R)
    
    Nknown=length(ID.B)
    n=sum(apply(y.B+y.L+y.R,1,sum)>0)
    
    #remove uncaptured individuals, put ID.B at top of y.L and y.R
    y.B.obs=y.B[ID.B,,]
    y.L.obs=y.L[c(ID.L),,]
    y.R.obs=y.R[c(ID.R),,]
    n.B=length(ID.B)
    n.L=length(ID.L)
    n.R=length(ID.R)
    if(n.B==1){
      y.B.obs=array(y.B.obs,dim=c(1,J,K))
    }
    if(n.L==1){
      y.L.obs=array(y.L.obs,dim=c(1,J,K))
    }
    if(n.R==1){
      y.R.obs=array(y.R.obs,dim=c(1,J,K))
    }
    
    out<-list(y.B.obs=y.B.obs,y.L.obs=y.L.obs,y.R.obs=y.R.obs,X=X,K=K,buff=buff,
              y.B=y.B,y.L=y.L,y.R=y.R,s=s,n=n,ID.B=ID.B,ID.L=ID.L,ID.R=ID.R)

    return(out)
  }