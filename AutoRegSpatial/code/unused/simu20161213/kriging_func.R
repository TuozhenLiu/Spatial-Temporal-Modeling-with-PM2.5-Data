

rho.func<-function(x, alpha)
{
  return(2*pnorm(x+alpha)-1)
}

kriging.typeII<-function(Y_train, rhos_train, dist_train, dist_pred, thetaE, thetaX)
{
  betaE = thetaE[1]; sig2E = thetaE[2]
  alpha = thetaX[1]; betaX = thetaX[2]; sig2X = thetaX[3]
  Time = ncol(Y_train)
  eps_train = Y_train[,-1] - rhos_train*Y_train[,-ncol(Y)]
  SigE = cov.exp(dist_train, beta = betaE, sig2 = sig2E)
  ce0 = cov.exp(dist_pred, beta = betaE, sig2 = sig2E)
  eps_pred = t(ce0)%*%solve(SigE)%*%eps_train
  
  X1_train = qnorm((rhos_train + 1)/2); #alpha = mean(X1_train)
  X_train = X1_train - alpha
  
  SigX = cov.exp(dist_train, beta = betaX, sig2 = sig2X)
  cx0 = cov.exp(dist_pred, beta = betaX, sig2 = sig2X)
  
  cx0SigX = crossprod(cx0, solve(SigX))
  mu_x0 = cx0SigX%*%X_train
  sig2x0 = SigX[1,1] - diag(cx0SigX%*%cx0)
  rand_x0 = apply(cbind(mu_x0, sqrt(sig2x0)), 1, 
                  function(x) return(rnorm(n = 1000, mean = x[1], sd = x[2])))
  rhos_pred = rho.func(rand_x0, alpha)
  rhosk_pred = sapply(1:(Time-2), function(k) colMeans(rhos_pred^k))
  Y_pred = rowSums(eps_pred*cbind(1, rhosk_pred))
  return(Y_pred)
}

kriging.typeI<-function(Y_train, rhos_train, dist_train, dist_pred, thetaE, thetaX)
{
  betaE = thetaE[1]; sig2E = thetaE[2]
  alpha = thetaX[1]; betaX = thetaX[2]; sig2X = thetaX[3]
  
  Time = ncol(Y_train)
  Y1 = Y_train - colMeans(Y_train)
  SigY = tcrossprod(Y1)/Time
  ce0 = cov.exp(dist_pred, beta = betaE, sig2 = sig2E)
  
  X1_train = qnorm((rhos_train + 1)/2); #alpha = mean(X1_train)
  X_train = X1_train - alpha
  
  SigX = cov.exp(dist_train, beta = betaX, sig2 = sig2X)
  cx0 = cov.exp(dist_pred, beta = betaX, sig2 = sig2X)
  
  cx0SigX = crossprod(cx0, solve(SigX))
  mu_x0 = cx0SigX%*%X_train
  sig2x0 = SigX[1,1] - diag(cx0SigX%*%cx0)
  rand_x0 = apply(cbind(mu_x0, sqrt(sig2x0)), 1, 
                  function(x) return(rnorm(n = 1000, mean = x[1], sd = x[2])))
  rhos_pred = rho.func(rand_x0, alpha)
  c1rho = sapply(rhos_train, function(x) colMeans(1/(1-x*rhos_pred)))
  Y_pred = (t(ce0)*c1rho)%*%solve(SigY)%*%Y_train[,1]
  return(Y_pred)
}

