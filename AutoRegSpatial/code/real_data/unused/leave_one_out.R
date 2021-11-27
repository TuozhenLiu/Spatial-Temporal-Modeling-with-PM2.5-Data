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
  for (est_win in c(90, 120, 150))
  {
    Yt_pred1 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    Yt_pred2 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    
    Yt_krig1 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    Yt_krig2 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    Yt_krig3 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    Yt_krig4 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    Yt_krig5 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    Yt_krig6 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    
    Time = ncol(Ymat)
    for (t in begin_ind:end_ind)
    {
      time_win = (Time - t - est_win + 1):(Time - t + 1)
      Ymat_win = Ymat[, time_win]
      for (i in 1:nrow(Ymat))
      {
        cat(t, i, "\r")
        pred_ind = i
        train_ind = setdiff(1:nrow(Ymat), pred_ind)
        
        dist_train = dist_loc2[train_ind, train_ind]
        dist_pred = dist_loc2[train_ind, pred_ind]
        Y_train = Ymat_win[train_ind,]
        
        rhosMu = estRhoMuIter(Y_train)
        mu = mean(rhosMu[1,])
        rhos_est = rhosMu[2,]
        
        #X1 = qnorm((rhos_est+1)/2)
        #alpha = mean(X1)
        #(thetaX = c(alpha, lse.X(X1 - alpha, dist_train)))
        thetaX = estThetaRho(rhos_est, dist_train)
        (thetaE = lse.theta(Y_train, dist_train, rhos = rhos_est, mu = mu))
        Yt_pred2[i, t-begin_ind+1] = predict.typeII(Y_train, rhos_est, dist_train, dist_pred, 
                                        thetaE, thetaX, mu = mu)
        Yt_pred1[i, t-begin_ind+1] = predict.typeI(Y_train, rhos_est, loc_train = loc[train_ind,], loc_pred = loc[pred_ind,],
                                       thetaE, thetaX, mu = mu)
        
        Yt_krig2[i, t-begin_ind+1] = krigingMu.typeII(Y_train, rhos_est, dist_train, dist_pred, 
                                        thetaE, thetaX, mu = mu)
        (Yt_krig3[i, t-begin_ind+1] = krigingMu.typeIII(Y_train, rhos_est, loc_train = loc[train_ind,], loc_pred = loc[pred_ind,],
                                                       thetaE, thetaX, mu = mu))
        Yt_krig4[i, t-begin_ind+1] = krigingMu.typeIV(Ymat_win, tr_ind = train_ind)
        Yt_krig6[i, t-begin_ind+1] = krigingMu.typeVI(Ymat_win, tr_ind = train_ind)
        
        Yt_krig5[i, t-begin_ind+1] = krigingMu.typeV(Y_train, rhos_est, loc_train = loc[train_ind,], loc_pred = loc[pred_ind,],
                                                     thetaE, thetaX, mu = mu)
        
        Yt_krig1[i, t-begin_ind+1] = krigingMu.typeI(Y_train, rhos_est, loc_train = loc[train_ind,], loc_pred = loc[pred_ind,],
                                       thetaE, thetaX, mu = mu)
        
      }
    }
    
    pred1_rmse = sqrt(mean((Yt_pred1 - Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    pred2_rmse = sqrt(mean((Yt_pred2 - Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    
    krig1_rmse = sqrt(mean((Yt_krig1 - Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    krig2_rmse = sqrt(mean((Yt_krig2 - Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    krig3_rmse = sqrt(mean((Yt_krig3 - Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    krig4_rmse = sqrt(mean((Yt_krig4 - Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    krig5_rmse = sqrt(mean((Yt_krig5 - Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    krig6_rmse = sqrt(mean((Yt_krig6 - Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    
    muY = mean(Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])
    sse = sqrt(mean((Ymat[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)]- muY)^2))
    
    cat("BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
        "Krig RMSE: ", krig1_rmse, krig2_rmse, krig3_rmse, krig4_rmse, krig6_rmse,"\n",
        "Krig RMSE R2: ", 1-(krig1_rmse/sse)^2, 1-(krig2_rmse/sse)^2, 1-(krig3_rmse/sse)^2, 1-(krig4_rmse/sse)^2, 1-(krig6_rmse/sse)^2, "\n",
        "Pred RMSE: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
  }
}


hist(Yt_pred2)

