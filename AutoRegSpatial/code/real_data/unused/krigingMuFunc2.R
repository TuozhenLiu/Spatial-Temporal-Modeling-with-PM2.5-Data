

krigingMu.typeIV<-function(Ymat, tr_ind)
{
  Time = ncol(Ymat)
  Y0 = Ymat - rowMeans(Ymat)
  SigY = tcrossprod(Y0[tr_ind,])
  cy0 = Y0[tr_ind,]%*%(Y0[-tr_ind,])
  muY = (mean(Ymat[tr_ind,]))
  Y_pred = t(cy0)%*%solve(SigY)%*%(Ymat[tr_ind,Time] - muY) + muY
  return(Y_pred)
}


krigingMu.typeVI<-function(Ymat, tr_ind)
{
  Time = ncol(Ymat)
  pred_ind = setdiff(1:nrow(Ymat), tr_ind)
  Y0 = Ymat - rowMeans(Ymat)
  
  Y1 = rbind(Y0[tr_ind,1:(Time-1)], Y0[tr_ind,2:Time])
  Y2 = c(Y0[pred_ind,2:Time])
  #tr_ind = c(1:n_train, (N+1):(N+n_train))
  SigY = tcrossprod(Y1)
  cy0 = Y1%*%Y2
  
  #SigY = tcrossprod(Y0[tr_ind,])
  #cy0 = Y0[tr_ind,]%*%(Y0[-tr_ind,])
  muY = (mean(Ymat[tr_ind,]))
  Y_pred = t(cy0)%*%solve(SigY)%*%(c(Ymat[tr_ind,Time-1],Ymat[tr_ind,Time]) - muY) + muY
  return(Y_pred)
}


krigingMu.typeV<-function(Y_train, rhos_train, loc_train, loc_pred, thetaE, thetaX, mu, Nrep = 20)
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
  
  cy1 = matrix(0, nrow = n_train, ncol = n_pred)
  SigY1 = matrix(0, nrow = n_train, ncol = n_train)
  Y1 = Y_train - rowMeans(Y_train)
  for (i in 1:Nrep)
  {
    #show(i)
    rhos = simu.rho(beta1 = betaX, sig2 = sig2X, N, alpha = alpha, dist_loc2)
    Y = simu.Y(beta1 = betaE, sigy2 = sig2E, N = N, Time = Time, 
               rhos = rhos, dist_loc2, mu)
    Y0 = Y - rowMeans(Y)
    #SigY1 = SigY1 + tcrossprod(Y0[1:n_train,])/Time
    if (n_train+1==N)
      cy1 = cy1 + colSums(Y0[(n_train+1):N,]*t(Y1[1:n_train,]))/Time
    else
      cy1 = cy1 + tcrossprod(Y1[1:n_train,], Y0[(n_train+1):N,])/Time
    
  }
  cy0 = cy1/Nrep
  #SigY = SigY1/Nrep
  muY = mean(Y_train)
  
  SigY = tcrossprod(Y1)/Time
  Y_pred = t(cy0)%*%solve(SigY)%*%(Y_train[,Time] - muY) + muY
  #hist(Y_pred)
  
  return(Y_pred)
}



