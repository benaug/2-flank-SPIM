NimModel <- nimbleCode({
  #detection function priors
  p0B ~ dunif(0,1) #not using here
  #if you do not want to share p0L and p0R
  # p0L~dunif(0,1)
  # p0R~dunif(0,1)
  # #if you want to share p0L and p0R
  p0S ~ dunif(0,1)
  p0L <- p0S
  p0R <- p0S
  #Make sure upper bound on sigma uniform prior is appropriate
  sigma ~ dunif(0,20)
  #data augmentation prior
  psi ~ dunif(0,1)
  #likelihoods (except for s priors)
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1] #extract activity center cell
    #categorical activity center likelihood for this cell, equivalent to zero's trick
    dummy.data[i] ~ dCell(InSS[s.cell[i]])
    d2[i,1:J] <- GetD2(s = s[i,1:2], X = X[1:J,1:2],z=z[i])
    pd.B[i,1:J] <- GetDetectionProb(p0=p0B,sigma=sigma,d2=d2[i,1:J],z=z[i])
    pd.L[i,1:J] <- GetDetectionProb(p0=p0L,sigma=sigma,d2=d2[i,1:J],z=z[i])
    pd.R[i,1:J] <- GetDetectionProb(p0=p0R,sigma=sigma,d2=d2[i,1:J],z=z[i])
    y.B.true[i,1:J] ~ dBernoulliVectorBoth(pd.B[i,1:J],K1D=K1D[1:J],J.cams=J.cams[1:J],z=z[i])
    y.L.true[i,1:J] ~ dBernoulliVectorSingle(pd.L[i,1:J],K1D=K1D[1:J],J.cams=J.cams[1:J],z=z[i])
    y.R.true[i,1:J] ~ dBernoulliVectorSingle(pd.R[i,1:J],K1D=K1D[1:J],J.cams=J.cams[1:J],z=z[i])
  }
  #have to trick Nimble to know ID is part of model by using it somewhere in model statement
  IDdummy <- IDdummyfun(ID.L=ID.L[1:n.L],ID.R=ID.R[1:n.R])
  N <- sum(z[1:M])
})# end model
