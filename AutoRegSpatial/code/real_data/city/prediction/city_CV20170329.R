
### 本文件为比较Type I and II kriging的结果

setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/")

### source the functions
#source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/krigingMuFunc.R") # kriging functions 
source("E:/OneDrive/1. Academic/spatial_autoreg/code/simulator.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/est_infer_func.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/predictFunc2.R") # prediction functions (similar to kriging functions)

load("rda/loc1.rda") # 读入经纬度数据
load("rda/Ymat.rda")
load("rda/Ymat_detrend.rda")

Ymat = Ymat_detrend
loc = loc1

### set the covariance function type
cov.type = "quad"
loc1 = loc[,3:4]*1
if (cov.type== "quad")
  loc1 = loc1*4
dist_loc2 = as.matrix(dist(loc1))


### predictive kriging using sliding window
begin_ind = 1
pred_win = ncol(Ymat) - 60 # use the first 60 days for training
set.seed(2345)
#set.seed(12)
Nloc = 100 # randomly select 100 cities for PK

for (begin_ind in c(1))
{
  end_ind = begin_ind + pred_win -1
  for (est_win in c(60))
  {
    Yt_pred1 = matrix(0, nrow = Nloc, ncol = pred_win)
    Yt_pred2 = matrix(0, nrow = Nloc, ncol = pred_win)
    Ymat_sub = matrix(0, nrow = Nloc, ncol = ncol(Ymat))
    
    Time = ncol(Ymat)
    for (t in begin_ind:end_ind)
    {
      pred_inds = sort(sample(1:nrow(loc1), Nloc))
      
      time_win = (Time - t - est_win + 1):(Time - t + 1)
      Ymat_win = Ymat[, time_win]
      
      cat(t, " | ")
      pred_ind = pred_inds
      i = pred_inds
      train_ind = setdiff(1:nrow(Ymat), pred_ind)
      
      dist_train = dist_loc2[train_ind, train_ind]
      dist_pred = dist_loc2[train_ind, pred_ind]
      Y_train = Ymat_win[train_ind,]
      
      rhosMu = estRhoMu(Y_train)
      mus_est = rhosMu[1,]
      rhos_est = rhosMu[2,]
      
      (thetaMu = estThetaMu(mus_est, dist_train, infer = F,  cov.type = cov.type))
      (thetaX = estThetaRho(rhos_est, dist_train, infer = F, cov.type = cov.type))
      (thetaE = lse.theta(Y_train, dist_train, rhos = rhos_est, mu = mus_est, infer = F, cov.type = cov.type))
      
      
      Ymat_sub[, Time - t + 1] = Ymat[pred_inds, Time - t + 1]
      Yt_pred2[, t-begin_ind+1] = predict.typeII(Y_train, rhos_est, mus_est, dist_train, dist_pred, 
                                                  thetaE, thetaX, thetaMu, cov.type = cov.type) 
      Yt_pred1[, t-begin_ind+1] = predict.typeI(Y_train, rhos_est, mus_est, loc_train = loc1[train_ind,], loc_pred = loc1[pred_ind,],
                                                 thetaE, thetaX, thetaMu, cov.type = cov.type)
      
      ind1 = 1:(t-begin_ind+1)
      ind2 = (Time -begin_ind+1):(Time - t + 1)
      cat("time win:")
      pred1_rmse = sqrt(median((Yt_pred1[, ind1] - Ymat_sub[, ind2])^2))
      pred2_rmse = sqrt(median((Yt_pred2[, ind1] - Ymat_sub[, ind2])^2))
      
      # krig1_rmse = sqrt(mean((Yt_krig1[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
      # krig2_rmse = sqrt(mean((Yt_krig2[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
      # krig3_rmse = sqrt(mean((Yt_krig3[pred_inds, t-begin_ind+1] - Ymat[pred_inds,Time - t + 1])^2))
      
      muY = mean(Ymat_sub[, ind2])
      sse = sqrt(mean((Ymat_sub[, ind2]- muY)^2))
      
      cat(#"BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
        "Time point: ", t, 
        " Pred RMSE: ", pred1_rmse, pred2_rmse, "\n",
        #"Krig RMSE R2: ", 1-(krig1_rmse/sse)^2, 1-(krig2_rmse/sse)^2, "\n",
        "Pred RMSE R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
      
      pred1_rmse = sqrt(mean((Yt_pred1[, ind1] - Ymat_sub[, ind2])^2))
      pred2_rmse = sqrt(mean((Yt_pred2[, ind1] - Ymat_sub[, ind2])^2))
      
      
      muY = mean(Ymat_sub[, ind2])
      sse = sqrt(mean((Ymat_sub[, ind2]- muY)^2))
      
      cat("Pred Mean R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
      
    }
    
    pred1_rmse = sqrt(mean((Yt_pred1 - Ymat_sub[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    pred2_rmse = sqrt(mean((Yt_pred2 - Ymat_sub[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    
    # krig1_rmse = sqrt(mean((Yt_krig1[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    # krig2_rmse = sqrt(mean((Yt_krig2[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    # #krig3_rmse = sqrt(mean((Yt_krig3[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    
    muY = mean(Ymat_sub[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])
    sse = sqrt(mean((Ymat_sub[,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)]- muY)^2))
    
    cat("BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
        #"Krig RMSE: ", krig1_rmse, krig2_rmse, "\n",
        #"Krig RMSE R2: ", 1-(krig1_rmse/sse)^2, 1-(krig2_rmse/sse)^2,  "\n",
        "Pred RMSE R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
    
    
  }
}

