NimModel <- nimbleCode({
  #detection function priors
  p0B~dunif(0,1)
  #if you do not want to share p0L and p0R
  p0L~dunif(0,1)
  p0R~dunif(0,1)
  # #if you want to share p0L and R
  # p0S ~ dunif(0,1)
  # p0L <- p0S
  # p0R <- p0S
  sigma~dunif(0,20)
  #data augmentation prior
  psi~dunif(0,1)
  #likelihoods (except for s priors)
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    d2[i,1:J] <- GetD2(s = s[i,1:2], X = X[1:J,1:2],z=z[i])
    pd.B[i,1:J] <- GetDetectionProb(p0=p0B,sigma=sigma,d2=d2[i,1:J],z=z[i])
    pd.L[i,1:J] <- GetDetectionProb(p0=p0L,sigma=sigma,d2=d2[i,1:J],z=z[i])
    pd.R[i,1:J] <- GetDetectionProb(p0=p0R,sigma=sigma,d2=d2[i,1:J],z=z[i])
    y.B.true[i,1:J,1:K] ~ dBernoulliVectorBoth(pd.B[i,1:J],K2D=K2D[1:J,1:K],z=z[i])
    y.L.true[i,1:J,1:K] ~ dBernoulliVectorSingle(pd.L[i,1:J],K2D=K2D[1:J,1:K],z=z[i])
    y.R.true[i,1:J,1:K] ~ dBernoulliVectorSingle(pd.R[i,1:J],K2D=K2D[1:J,1:K],z=z[i])
  }
  #have to trick Nimble to know ID is part of model by using it somewhere in model statement
  IDdummy <- IDdummyfun(ID.L=ID.L[1:n.L])
  N <- sum(z[1:M])
})# end model
