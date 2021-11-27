

### this file is for predictive kriging

predict.typeII<-function(Y_train, rhos_train, mus_train, dist_train, dist_pred, thetaE, thetaX, thetaMu, 
                         cov.type = "exp", mixed  = F)
{
  if (mixed == T)
  {
    cov.func1 = cov.exp
    cov.func2 = cov.quad
  }
  else
  {
    if (cov.type=="exp")
      cov.func1 = cov.exp
    else
      cov.func1 = cov.quad
    cov.func2 = cov.func1
  }
  
  
  betaE = thetaE[1]; sig2E = thetaE[2]
  alpha = thetaX[1]; betaX = thetaX[2]; sig2X = thetaX[3]
  mu = thetaMu[1]; betaMu = thetaMu[2]; sig2Mu = thetaMu[3]
  
  Time = ncol(Y_train)
  eps_train = Y_train[,-1] - rhos_train*Y_train[,-Time]-mus_train
  SigE = cov.func2(dist_train, beta = betaE, sig2 = sig2E)
  ce0 = cov.func2(dist_pred, beta = betaE, sig2 = sig2E)
  eps_pred = t(ce0)%*%solve(SigE)%*%eps_train
  eps_pred[,Time-1] = 0
  
  SigMu = cov.func1(dist_train, beta = betaMu, sig2 = sig2Mu)
  cmu0 = cov.func1(dist_pred, beta = betaMu, sig2 = sig2Mu)
  mus0_train = mus_train - mu
  mu_pred = t(cmu0)%*%solve(SigMu)%*%mus0_train + mu
  
  rhos_train[which(abs(rhos_train)>1)] = sign(rhos_train[which(abs(rhos_train)>1)])*rep(0.999, sum(abs(rhos_train)>1))
  y1 = (rhos_train + 1)/2
  X1_train = log(y1/(1-y1))#qnorm((rhos_train + 1)/2); #alpha = mean(X1_train)
  X_train = X1_train - alpha
  
  SigX = cov.func1(dist_train, beta = betaX, sig2 = sig2X)
  cx0 = cov.func1(dist_pred, beta = betaX, sig2 = sig2X)
  
  cx0SigX = crossprod(cx0, solve(SigX))
  mu_x0 = cx0SigX%*%X_train
  sig2x0 = SigX[1,1] - diag(cx0SigX%*%cx0)
  sig2x0[sig2x0<0] = 10^(-5)
  rand_x0 = apply(cbind(mu_x0, sqrt(sig2x0)), 1, 
                  function(x) return(rnorm(n = 2000, mean = x[1], sd = x[2])))
  rhos_pred = rho.func(rand_x0, alpha)
  rhos_inv_pred = colMeans(1/(1-rhos_pred))
  
  ncol = Time-2
  rhosk_pred = sapply(ncol:1, function(k) colMeans(rhos_pred^k))
  if (is.null(dim(dist_pred)))
    Y_pred = sum(eps_pred[, 1:(ncol+1)]*c(rhosk_pred, 1)[1:(ncol+1)]) + rhos_inv_pred*mu_pred
  else
    Y_pred = rowSums(eps_pred[, 1:(ncol+1)]*cbind(rhosk_pred, 1)[,1:(ncol+1)]) + rhos_inv_pred*mu_pred
  return(Y_pred)
}

# predict.typeI<-function(Y_train, rhos_train, loc_train, loc_pred, thetaE, thetaX, mu, Nrep = 10)
# {
#   Time = ncol(Y_train)
#   betaE = thetaE[1]; sig2E = thetaE[2]
#   alpha = thetaX[1]; betaX = thetaX[2]; sig2X = thetaX[3]
#   
#   n_train = nrow(loc_train)
#   n_pred = nrow(loc_pred)
#   loc = rbind(loc_train, loc_pred)
#   dist_loc2 = as.matrix(dist(loc))
#   #dist_train = dist_loc2[1:n_train, 1:n_train]
#   #dist_pred = dist_loc2[1:n_train, (n_train+1):nrow(loc)]
#   N = nrow(loc)
#   
#   cy1 = matrix(0, nrow = n_train, ncol = n_pred)
#   SigY1 = matrix(0, nrow = n_train, ncol = n_train)
#   for (i in 1:Nrep)
#   {
#     #show(i)
#     rhos = simu.rho(beta = betaX, sig2 = sig2X, N, alpha = alpha, dist_loc2)
#     Y = simu.Y(beta = betaE, sigy2 = sig2E, N = N, Time = Time, 
#                rhos = rhos, dist_loc2, mu)
#     Y0 = Y - rowMeans(Y)
#     SigY1 = SigY1 + tcrossprod(Y0[1:n_train,])/Time
#     if (n_train+1==N)
#       cy1 = cy1 + colSums(Y0[(n_train+1):N,]*t(Y0[1:n_train,]))/Time
#     else
#       cy1 = cy1 + tcrossprod(Y0[1:n_train,], Y0[(n_train+1):N,])/Time
#     
#   }
#   cy0 = cy1/Nrep
#   SigY = SigY1/Nrep
#   muY = mean(Y_train[,Time])
#   #Y1 = Y_train - rowMeans(Y_train)
#   #SigY = tcrossprod(Y1)/Time
#   Yhat_train = mu + rhos_train*Y_train[,Time-1]
#   Y_pred = t(cy0)%*%solve(SigY)%*%(Yhat_train - muY) + muY
#   #hist(Y_pred)
#   
#   return(Y_pred)
# }


### This is a faster version
predict.typeI<-function(Y_train, rhos_train, mus_train, loc_train, loc_pred, thetaE, thetaX, thetaMu, 
                        cov.type = "exp", Nrep = 10, mixed = F)
{
  if (!mixed)
  {
    cov.type1 = cov.type
    cov.type2 = cov.type
  }
  else
  {
    cov.type1 = "exp"
    cov.type2 = "quad"
  }
    
  Time = ncol(Y_train)
  betaE = thetaE[1]; sig2E = thetaE[2]
  alpha = thetaX[1]; betaX = thetaX[2]; sig2X = thetaX[3]
  mu = thetaMu[1]; betaMu = thetaMu[2]; sig2Mu = thetaMu[3]
  
  n_train = nrow(loc_train)
  n_pred = nrow(loc_pred)
  loc = rbind(loc_train, loc_pred)
  dist_loc2 = as.matrix(dist(loc))
  #dist_train = dist_loc2[1:n_train, 1:n_train]
  #dist_pred = dist_loc2[1:n_train, (n_train+1):nrow(loc)]
  N = nrow(loc)
  
  rhos_R = simu.rho.rep(beta = betaX, sig2 = sig2X, N, Nrep, alpha = alpha, dist_loc2, cov.type1)
  mus_R = simu.mu.rep(beta = betaMu, sig2 = sig2Mu, N, Nrep, mu = mu, dist_loc2, cov.type1)
  Y_R = simu.Y.rep(beta = betaE, sigy2 = sig2E, N, Nrep, Time, rhos_R, dist_loc2, mus_R, cov.type2)
  Y_R0 = Y_R - rowMeans(Y_R)
  rr = rep(1:Nrep, each = N)
  Y_R0_list = split.data.frame(Y_R0, factor(rr))
  SigY = Reduce("+", lapply(Y_R0_list, function(Y0) tcrossprod(Y0[1:n_train,])/Time))/Nrep
  if (n_train+1==N)
    cy0 = Reduce("+", lapply(Y_R0_list, function(Y0) colSums(Y0[(n_train+1):N,]*t(Y0[1:n_train,]))/Time))/Nrep
  else
    cy0 = Reduce("+", lapply(Y_R0_list, function(Y0) tcrossprod(Y0[1:n_train,], Y0[(n_train+1):N,])/Time))/Nrep
  
  muY = mean(Y_train[,Time])
  #Y1 = Y_train - rowMeans(Y_train)
  #SigY = tcrossprod(Y1)/Time
  Yhat_train = mus_train + rhos_train*Y_train[,Time-1]
  Y_pred = t(cy0)%*%solve(SigY)%*%(Yhat_train - muY) + muY
  #hist(Y_pred)
  return(Y_pred)
}
