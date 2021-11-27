


lse.theta0<-function(Y, dist_loc2)
{
  Y = Y - mean(Y)
  sig_hat = tcrossprod(Y)
  
  theta = c(0.1, 20)
  iter = 1
  del = 1
  while(mean(abs(del))>10^-3&iter<1000)
  {
    #cat(mean(abs(del)), " ", theta, "\n")
    del = lse.step(sig_hat, beta1 = theta[1], sigy2 = theta[2], dist_loc2)
    theta = theta - del*0.5
    if (any(theta<0|theta>10))
      theta = runif(2, 0.1,10)
    iter = iter+1
  }
  return(abs(theta))
}
krigingST<-function(Y_train, loc_train, loc_pred, thetaY)
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
  if (nrow(loc_pred) ==1)
  {
    cy0 = matrix(cy0, nrow = 1)
    Y_pred = (cy0)%*%solve(SigY0)%*%(Y_train - muY) + muY 
  }
  else
    Y_pred = (cy0)%*%solve(SigY0)%*%(Y_train - muY) + muY
  return(Y_pred)
}


dat_sub$day = as.numeric(dat_sub$date)-min(as.numeric(dat_sub$date))+1
dat_sub1 = merge(dat_sub[,c("station", "PM2_5", "day", "month")], station_rho_loc[,1:3], by = "station")
dat_sub1[,c("lat", "lon")] = 10*dat_sub1[,c("lat", "lon")] 
dat_sub1$day1 = dat_sub1$day/365*10
dat_sub1 = dat_sub1[complete.cases(dat_sub1),]
dat_sub1$Y = log(dat_sub1$PM2_5)
summary(dat_sub1)

###
krigST = matrix(0, nrow = 32, ncol = length(334:364))
for (dd in 334:364)
{
  ind = which(is.element(dat_sub1$day, (dd-30):dd))
  dat_sub2 = dat_sub1[ind,]
  loct = dat_sub2[,5:7]
  Y = dat_sub2$Y
  lastday_ind = which(dat_sub2$day == max(dat_sub2$day))
  
  dist_loc2 = as.matrix(dist(loct))
  (thetaY = lse.theta0(Y, dist_loc2))
  for (i in 1:32)
  {
    cat(dd, i, "\r")
    loc_train = loct[-lastday_ind[i],]
    loc_pred = loct[lastday_ind[i],]
    krigST[i, dd-334+1] = krigingST(Y_train = Y[-lastday_ind[i]], loc_train, loc_pred, thetaY)
  }
}


pred_win = 30
set.seed(1234)
for (begin_ind in c(1, 60, 90, 120, 150))
{
  end_ind = begin_ind + pred_win -1
  for (est_win in c(30))
  {
    
    Yt_krigST = matrix(0, nrow = 32, ncol = pred_win)
    
    Time = max(dat_sub1$day)
    for (t in begin_ind:end_ind)
    {
      time_win = (Time - t - est_win + 1):(Time - t + 1)
      #ind = which(is.element(dat_sub1$day, (dd-30):dd))
      ind = which(is.element(dat_sub1$day, time_win))
      dat_sub2 = dat_sub1[ind,]
      loct = dat_sub2[,5:7]
      Y = dat_sub2$Y
      lastday_ind = which(dat_sub2$day == max(dat_sub2$day))
      lastday_match_ind = match(dat_sub2$station[lastday_ind], station_rho_loc$station)
      
      dist_loc2 = as.matrix(dist(loct))
      (thetaY = lse.theta0(Y, dist_loc2))
      for (i in 1:length(lastday_match_ind))
      {
        ii = lastday_match_ind[i]
        loc_train = loct[-lastday_ind[i],]
        loc_pred = loct[lastday_ind[i],]
        Yt_krigST[ii, t-begin_ind+1] = krigingST(Y_train = Y[-lastday_ind[i]], loc_train, loc_pred, thetaY)
        cat(t, i, Yt_krigST[ii, t-begin_ind+1], "\r")
      }
    }
    Yt_krigST[Yt_krigST==0] = NA
    krigST_rmse = sqrt(mean((Yt_krigST - Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2, na.rm = T))
    muY = mean(Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)], na.rm = T)
    sse = sqrt(mean((Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)]- muY)^2, na.rm = T))
    
    cat("BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
        "Krig RMSE: ", krigST_rmse,"\n",
        "Krig RMSE R2: ", 1-(krigST_rmse/sse)^2, "\n")
  }
}

