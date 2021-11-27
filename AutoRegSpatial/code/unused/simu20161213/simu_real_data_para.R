
### 本文件为比较Type I and II kriging的结果

setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/")

### source the functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/krigingMuFunc.R") # kriging functions 
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/city/realFunc_city.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/predictFunc.R") # prediction functions (similar to kriging functions)



######################################################################################################

set.seed(1234)

mu = 0.85; betaE = 0.5; sig2E = 0.21
betaX = 0.53; sig2X = 0.03; alpha = 2.05

Ns = c(100, 200, 400)
res_rmse = list()
res_r2 = list()
set.seed(2345)
Nrep = 10
begin_ind = 1 # 从最后一期往前数10期
est_win = c(Time - pred_win)
pred_win = 10

#N  = 200; Time = 200
for (k in 1:length(Ns))
{
  N = Ns[k]
  Time = N
  cat("N: ", N, " Time: ", Time, "\n")
  perf_r2 = matrix(0, nrow = Nrep, ncol = 4)
  perf_rmse = matrix(0, nrow = Nrep, ncol = 4)
  for (r in 1:Nrep)
  {
    loc1 = data.frame(simu.loc(N))
    dist_loc2 = as.matrix(dist(loc1))
    rhos = simu.rho(beta1 = betaX, sig2 = sig2X, N = nrow(loc1), alpha = alpha, dist_loc2)
    Ymat = simu.Y(beta1 = betaE, sigy2 = sig2E, N = nrow(loc1), Time = Time, rhos = rhos, dist_loc2, mu = mu)
    #pred_win = 60
    
    #set.seed(2345)
    pred_inds = sort(sample(1:nrow(loc1), 1))
    
    for (begin_ind in c(1))
    {
      end_ind = begin_ind + pred_win -1
      for (est_win in c(Time - pred_win))
      {
        
        Yt_krig1 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
        Yt_krig2 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
        Yt_pred1 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
        Yt_pred2 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
        
        for (t in begin_ind:end_ind)
        {
          time_win = (Time - t - est_win + 1):(Time - t + 1)
          Ymat_win = Ymat[, time_win]
          for (i in pred_inds)
          {
            #cat(t, i, " | ")
            pred_ind = i
            train_ind = setdiff(1:nrow(Ymat), pred_ind)
            
            dist_train = dist_loc2[train_ind, train_ind]
            dist_pred = dist_loc2[train_ind, pred_ind]
            Y_train = Ymat_win[train_ind,]
            
            rhosMu = estRhoMuIter(Y_train)
            mu_est = mean(rhosMu[1,])
            rhos_est = rhosMu[2,]
            
            (thetaX = estThetaRho(rhos_est, dist_train))
            (thetaE = lse.theta(Y_train, dist_train, rhos = rhos_est, mu = mu_est))
            
            Yt_krig2[i, t-begin_ind+1] = krigingMu.typeII(Y_train, rhos_est, dist_train, dist_pred, 
                                                          thetaE, thetaX, mu = mu_est)
            
            Yt_krig1[i, t-begin_ind+1] = krigingMu.typeI(Y_train, rhos_est, loc_train = loc1[train_ind,], loc_pred = loc1[pred_ind,],
                                                         thetaE, thetaX, mu = mu_est)
            
            Yt_pred2[i, t-begin_ind+1] = predict.typeII(Y_train, rhos_est, dist_train, dist_pred, 
                                                        thetaE, thetaX, mu = mu_est)
            Yt_pred1[i, t-begin_ind+1] = predict.typeI(Y_train, rhos_est, loc_train = loc1[train_ind,], loc_pred = loc1[pred_ind,],
                                                       thetaE, thetaX, mu = mu_est)
            
            
          }
          ind1 = 1:(t-begin_ind+1)
          ind2 = (Time -begin_ind+1):(Time - t + 1)
          pred1_rmse = sqrt(mean((Yt_pred1[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
          pred2_rmse = sqrt(mean((Yt_pred2[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
          
          krig1_rmse = sqrt(mean((Yt_krig1[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
          krig2_rmse = sqrt(mean((Yt_krig2[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
          
          muY = mean(Ymat[pred_inds, ind2])
          sse = sqrt(mean((Ymat[pred_inds, ind2]- muY)^2))
          
          #cat(#"BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
          #  "Time point: ", t, 
          #  " Krig RMSE: ", krig1_rmse, krig2_rmse, "\n",
          #  "Krig RMSE R2: ", 1-(krig1_rmse/sse)^2, 1-(krig2_rmse/sse)^2, "\n",
          #  "Pred RMSE R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
          cat(t, "\r")
          
        }
        
        pred1_rmse = sqrt(mean((Yt_pred1[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
        pred2_rmse = sqrt(mean((Yt_pred2[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
        
        krig1_rmse = sqrt(mean((Yt_krig1[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
        krig2_rmse = sqrt(mean((Yt_krig2[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
        #krig3_rmse = sqrt(mean((Yt_krig3[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
        
        muY = mean(Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])
        sse = sqrt(mean((Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)]- muY)^2))
        
        cat("BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
            #"Krig RMSE: ", krig1_rmse, krig2_rmse, "\n",
            "Krig RMSE R2: ", 1-(krig1_rmse/sse)^2, 1-(krig2_rmse/sse)^2,  "\n",
            "Pred RMSE R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
        
      }
    }
    perf_r2[r,] = c(1-(krig1_rmse/sse)^2, 1-(krig2_rmse/sse)^2, 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2)
    perf_rmse[r,] = c(krig1_rmse, krig2_rmse, pred1_rmse, pred2_rmse)
  }
  cat("\n N: ", N, " Time: ", Time, "\n",
      "Krig RMSE R2: ", colMeans(perf_r2[,1:2]),  "\n",
      "Pred RMSE R2: ", colMeans(perf_r2[,3:4]), "\n")
  res_r2[[k]] = perf_r2
  res_rmse[[k]] = perf_rmse
}













dif_mat = Yt_krig2[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)]
sse_mat = scale(Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)], center = T, scale = F)
sd(1 - colMeans(dif_mat^2)/colMeans(sse_mat^2))


sd(sqrt(colMeans(dif_mat^2)))