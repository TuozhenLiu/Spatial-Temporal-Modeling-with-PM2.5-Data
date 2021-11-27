
# rho function form
rho.func<-function(x, alpha)
{
  return(2*exp(x+alpha)/(1+exp(x+alpha))-1)
}

# type II kriging: optimal kriging
krigingMu.typeII<-function(Y_train, rhos_train, dist_train, dist_pred, thetaE, thetaX, mu)
{
  # obtain corresponding parameters
  betaE = thetaE[1]; sig2E = thetaE[2]
  alpha = thetaX[1]; betaX = thetaX[2]; sig2X = thetaX[3]
  Time = ncol(Y_train)
  
  eps_train = Y_train[,-1] - rhos_train*Y_train[,-Time]-mu # obtain residuals at each time point
  SigE = cov.exp(dist_train, beta = betaE, sig2 = sig2E) # use exponential covariance model to obtain covariance values
  ce0 = cov.exp(dist_pred, beta = betaE, sig2 = sig2E) # calculate covariance between s1,...,sN and s0 (s0 can be a vector)
  eps_pred = t(ce0)%*%solve(SigE)%*%eps_train # get kriging value for s0
  
  # obtain X by rhos
  y1 = (rhos_train + 1)/2
  X1_train = log(y1/(1-y1))#qnorm((rhos_train + 1)/2); #alpha = mean(X1_train)
  X_train = X1_train - alpha
  
  # similarly obtain covariance for X for s0, s1,..., sN
  SigX = cov.exp(dist_train, beta = betaX, sig2 = sig2X)
  cx0 = cov.exp(dist_pred, beta = betaX, sig2 = sig2X)
  cx0SigX = crossprod(cx0, solve(SigX)) 
  
  mu_x0 = cx0SigX%*%X_train # mu_x at s0
  sig2x0 = SigX[1,1] - diag(cx0SigX%*%cx0) # sigma_x^2 at s0
  sig2x0[sig2x0<0] = 10^(-5) # avoid the negative values
  rand_x0 = apply(cbind(mu_x0, sqrt(sig2x0)), 1, 
                  function(x) return(rnorm(n = 1000, mean = x[1], sd = x[2]))) # generate random samples of X(s0)
  rhos_pred = rho.func(rand_x0, alpha) # then obtain the random samples of rho(X(s0)) (recorded as rho0 in the following)
  rhos_inv_pred = colMeans(1/(1-rhos_pred)) # obtain random samples of 1/(1-rho0)
  
  ncol = Time-2
  rhosk_pred = sapply(ncol:1, function(k) colMeans(rhos_pred^k)) # obtain random samples of rho0^k
  #muY = mean(Y_train[,Time])
  # kriging on unobserved locations
  if (is.null(dim(dist_pred)))
    Y_pred = sum(eps_pred[, 1:(ncol+1)]*c(rhosk_pred, 1)) + rhos_inv_pred*mu # if there is only one kriging location
  else
    Y_pred = rowSums(eps_pred[, 1:(ncol+1)]*cbind(rhosk_pred, 1)) + rhos_inv_pred*mu # for >=2 kriging 
  return(Y_pred)
}

krigingMu.typeI<-function(Y_train, rhos_train, loc_train, loc_pred, thetaE, thetaX, mu, Nrep = 20)
{
  # obtain corresponding parameters
  Time = ncol(Y_train)
  betaE = thetaE[1]; sig2E = thetaE[2]
  alpha = thetaX[1]; betaX = thetaX[2]; sig2X = thetaX[3]
  
  # split training and testing dataset
  n_train = nrow(loc_train)
  n_pred = nrow(loc_pred)
  loc = rbind(loc_train, loc_pred)
  dist_loc2 = as.matrix(dist(loc))
  N = nrow(loc)
  
  # generate Y samples to calculate covariances for Y
  cy1 = matrix(0, nrow = n_train, ncol = n_pred)
  SigY1 = matrix(0, nrow = n_train, ncol = n_train)
  for (i in 1:Nrep)
  {
    #show(i)
    rhos = simu.rho(beta1 = betaX, sig2 = sig2X, N, alpha = alpha, dist_loc2)
    Y = simu.Y(beta1 = betaE, sigy2 = sig2E, N = N, Time = Time, 
               rhos = rhos, dist_loc2, mu)
    Y0 = Y - rowMeans(Y)
    SigY1 = SigY1 + tcrossprod(Y0[1:n_train,])/Time
    if (n_train+1==N)
      cy1 = cy1 + colSums(Y0[(n_train+1):N,]*t(Y0[1:n_train,]))/Time
    else
      cy1 = cy1 + tcrossprod(Y0[1:n_train,], Y0[(n_train+1):N,])/Time
    
  }
  # obtain covariance of Y to obtain Cy0 and SigY
  cy0 = cy1/Nrep
  SigY = SigY1/Nrep
  muY = mean(Y_train[,Time])
  Y_pred = t(cy0)%*%solve(SigY)%*%(Y_train[,Time] - muY) + muY # the best linear predictor
  
  return(Y_pred)
}


### This is a faster version
krigingMu.typeI<-function(Y_train, rhos_train, loc_train, loc_pred, thetaE, thetaX, mu, Nrep = 10)
{
  # obtain corresponding parameters
  Time = ncol(Y_train)
  betaE = thetaE[1]; sig2E = thetaE[2]
  alpha = thetaX[1]; betaX = thetaX[2]; sig2X = thetaX[3]
  
  # split training and testing dataset
  n_train = nrow(loc_train)
  n_pred = nrow(loc_pred)
  loc = rbind(loc_train, loc_pred)
  dist_loc2 = as.matrix(dist(loc))
  N = nrow(loc)
  
  # generate Y samples to calculate covariances for Y
  rhos_R = simu.rho.rep(beta1 = betaX, sig2 = sig2X, N, Nrep, alpha = alpha, dist_loc2)
  Y_R = simu.Y.rep(beta1 = betaE, sigy2 = sig2E, N, Nrep, Time, rhos_R, dist_loc2, mu)
  Y_R0 = Y_R - rowMeans(Y_R)
  rr = rep(1:Nrep, each = N)
  Y_R0_list = split.data.frame(Y_R0, factor(rr))
  SigY = Reduce("+", lapply(Y_R0_list, function(Y0) tcrossprod(Y0[1:n_train,])/Time))/Nrep
  if (n_train+1==N)
    cy0 = Reduce("+", lapply(Y_R0_list, function(Y0) colSums(Y0[(n_train+1):N,]*t(Y0[1:n_train,]))/Time))/Nrep
  else
    cy0 = Reduce("+", lapply(Y_R0_list, function(Y0) tcrossprod(Y0[1:n_train,], Y0[(n_train+1):N,])/Time))/Nrep
  
  muY = mean(Y_train[,Time])
  Yhat_train = mu + rhos_train*Y_train[,Time-1]
  Y_pred = t(cy0)%*%solve(SigY)%*%(Y_train[,Time] - muY) + muY # the best linear predictor
  
  return(Y_pred)
}



krigingMu.typeIII<-function(Y_train, rhos_train, loc_train, loc_pred, thetaE, thetaX, mu, Nrep = 20)
{
  Time = ncol(Y_train)
  betaE = thetaE[1]; sig2E = thetaE[2]
  alpha = thetaX[1]; betaX = thetaX[2]; sig2X = thetaX[3]
  
  n_train = nrow(loc_train)
  n_pred = nrow(loc_pred)
  loc = rbind(loc_train, loc_pred)
  dist_loc2 = as.matrix(dist(loc))
  #dist_train = dist_loc2[1:n_train, 1:n_train]
  #dist_pred = dist_loc2[1:n_train, (n_train+1):nrow(loc)]
  N = nrow(loc)
  train_ind = 1:n_train
  pred_ind = (n_train+1):N
  
  cy1 = matrix(0, nrow = n_train*2, ncol = n_pred)
  SigY1 = matrix(0, nrow = n_train*2, ncol = n_train*2)
  for (i in 1:Nrep)
  {
    #show(i)
    rhos = simu.rho(beta1 = betaX, sig2 = sig2X, N, alpha = alpha, dist_loc2)
    Y = simu.Y(beta1 = betaE, sigy2 = sig2E, N = N, Time = Time, 
               rhos = rhos, dist_loc2, mu)
    Y0 = Y - rowMeans(Y)
    Y1 = rbind(Y0[train_ind,1:(Time-1)], Y0[train_ind,2:Time])
    Y2 = c(Y0[pred_ind,2:Time], Y0[pred_ind,2:Time])
    #tr_ind = c(1:n_train, (N+1):(N+n_train))
    SigY1 = SigY1 + tcrossprod(Y1)/(Time-1)
    
    if (n_train+1==N)
      cy1 = cy1 + colSums(Y2*t(Y1))/(Time - 1)
    else
      cy1 = cy1 + tcrossprod(Y1, Y2)/(Time - 1)
    
  }
  cy0 = cy1/Nrep
  SigY = SigY1/Nrep
  muY = mean( mean(Y_train[,Time]))
  #Y1 = Y_train - rowMeans(Y_train)
  #SigY = tcrossprod(Y1)/Time
  Y_pred = t(cy0)%*%solve(SigY)%*%(c(Y_train[,Time-1],Y_train[,Time]) - muY) + muY
  #hist(Y_pred)
  
  return(Y_pred)
}



