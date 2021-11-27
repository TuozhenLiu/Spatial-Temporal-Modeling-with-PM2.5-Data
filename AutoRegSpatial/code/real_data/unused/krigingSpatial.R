
### boxplot everyday
for (i in 1:12)
{
  show(i)
  png(paste0("E:/OneDrive/1. Academic/spatial_autoreg/report/boxplot/LogPM25_month",i,".png"), width = 800, height = 480)
  boxplot(log(dat_sub$PM2_5[dat_sub$month1==i])~dat_sub$day1[dat_sub$month1==i])
  dev.off()
}


lse.theta0<-function(Y, dist_loc2)
{
  Y = Y - mean(Y)
  sig_hat = tcrossprod(Y)
  
  theta = c(2, 2)
  iter = 1
  del = 1
  while(mean(abs(del))>10^-3&iter<1000)
  {
    #cat(mean(abs(del)), " ", theta, "\n")
    del = lse.step(sig_hat, beta1 = theta[1], sigy2 = theta[2], dist_loc2)
    theta = theta - del*0.2
    if (any(theta<0|theta>10))
      theta = runif(2, 0.1,5)
    iter = iter+1
  }
  return(abs(theta))
}


krigingS<-function(Y_train, loc_train, loc_pred, thetaY)
{
  betaY = thetaY[1]; sig2Y = thetaY[2]
  n_train = nrow(loc_train)
  loc = rbind(loc_train, loc_pred)
  dist_loc2 = as.matrix(dist(loc))
  N = nrow(loc)
  SigY = cov.exp(dist_loc2, beta = betaY, sig2 = sig2Y)
  cy0 = SigY[(n_train+1):N, 1:n_train]
  SigY0 = SigY[1:n_train, 1:n_train]
  
  muY = mean(Y_train)
  if (is.vector(loc_pred))
  {
    cy0 = matrix(cy0, nrow = 1)
    Y_pred = (cy0)%*%solve(SigY0)%*%(Y_train - muY) + muY 
  }
  else
    Y_pred = (cy0)%*%solve(SigY0)%*%(Y_train - muY) + muY
  return(Y_pred)
}



(thetaY = lse.theta0(Ymat[,ncol(Ymat)], dist_loc2))


station_rho_loc = merge(loc, rhos_df, by = "station")
loc = station_rho_loc[,2:3]*10
dist_loc2 = as.matrix(dist(loc))

#pred_win = 60
begin_ind = 1
pred_win = 30


set.seed(1234)
for (begin_ind in c(1, 60, 90, 120, 150))
{
  end_ind = begin_ind + pred_win -1
  for (est_win in c(90))
  {
    
    Yt_krigs = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    
    Time = ncol(Ymat)
    for (t in begin_ind:end_ind)
    {
      time_win = (Time - t - est_win + 1):(Time - t + 1)
      Ymat_win = Ymat[, time_win]
      for (i in 1:nrow(Ymat))
      {
        cat(t, i, "\r")
        
        Y_train = Ymat_win[-i,ncol(Ymat_win)]
        
        (thetaY = lse.theta0(Y_train, dist_loc2[-i,-i]))
        
        Yt_krigs[i, t-begin_ind+1] = krigingS(Y_train, 
                                              loc_train = loc[-i,], 
                                              loc_pred = loc[i,], 
                                              thetaY)
        
      }
    }
    
    krigs_rmse = sqrt(mean((Yt_krigs - Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    muY = mean(Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])
    sse = sqrt(mean((Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)]- muY)^2))
    
    cat("BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
        "Krig RMSE: ", krigs_rmse,"\n",
        "Krig RMSE R2: ", 1-(krigs_rmse/sse)^2, "\n")
  }
}


