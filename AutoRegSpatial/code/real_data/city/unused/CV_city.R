
### 本文件为比较Type I and II kriging的结果

setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/")

### source the functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/krigingMuFunc.R") # kriging functions 
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/city/realFunc_city.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/predictFunc.R") # prediction functions (similar to kriging functions)

load("rda/loc.rda") # 读入经纬度数据
load("rda/Ymat.rda")
load("rda/rhos_df.rda")
set.seed(1234)
city_ind = (1:nrow(rhos_df))
Ymat = Ymat[city_ind,]
rhos_df = rhos_df[city_ind,]

station_rho_loc = merge(rhos_df,loc, by = "citycode")
Ymat = Ymat[match(station_rho_loc$citycode, rhos_df$citycode),]

loc1 = station_rho_loc[,c(4:5)]
dist_loc2 = as.matrix(dist(loc1))


#pred_win = 60
begin_ind = 1
pred_win = 300
#set.seed(2345)
pred_inds = sort(sample(1:nrow(loc1), 60))

set.seed(2345)
for (begin_ind in c(1))
{
  end_ind = begin_ind + pred_win -1
  for (est_win in c(60))
  {
    
    Yt_krig1 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    Yt_krig2 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    #Yt_krig3 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    Yt_pred1 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    Yt_pred2 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    
    Time = ncol(Ymat)
    for (t in begin_ind:end_ind)
    {
      time_win = (Time - t - est_win + 1):(Time - t + 1)
      Ymat_win = Ymat[, time_win]
      
        cat(t, " | ")
        pred_ind = pred_inds
        i = pred_inds
        train_ind = setdiff(1:nrow(Ymat), pred_ind)
        
        dist_train = dist_loc2[train_ind, train_ind]
        dist_pred = dist_loc2[train_ind, pred_ind]
        Y_train = Ymat_win[train_ind,]
        
        rhosMu = estRhoMuIter(Y_train)
        mu = mean(rhosMu[1,])
        rhos_est = rhosMu[2,]
        
        (thetaX = estThetaRho(rhos_est, dist_train))
        (thetaE = lse.theta(Y_train, dist_train, rhos = rhos_est, mu = mu))
        
        Yt_krig2[i, t-begin_ind+1] = krigingMu.typeII(Y_train, rhos_est, dist_train, dist_pred, 
                                                      thetaE, thetaX, mu = mu)
        #(Yt_krig3[i, t-begin_ind+1] = krigingMu.typeIII(Y_train, rhos_est, loc_train = loc1[train_ind,], loc_pred = loc1[pred_ind,],
        #                                                thetaE, thetaX, mu = mu))
        
        Yt_krig1[i, t-begin_ind+1] = krigingMu.typeI(Y_train, rhos_est, loc_train = loc1[train_ind,], loc_pred = loc1[pred_ind,],
                                                     thetaE, thetaX, mu = mu)
        
        Yt_pred2[i, t-begin_ind+1] = predict.typeII(Y_train, rhos_est, dist_train, dist_pred, 
                                                    thetaE, thetaX, mu = mu)
        Yt_pred1[i, t-begin_ind+1] = predict.typeI(Y_train, rhos_est, loc_train = loc1[train_ind,], loc_pred = loc1[pred_ind,],
                                                   thetaE, thetaX, mu = mu)
        
        
      
      ind1 = 1:(t-begin_ind+1)
      ind2 = (Time -begin_ind+1):(Time - t + 1)
      pred1_rmse = sqrt(median((Yt_pred1[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
      pred2_rmse = sqrt(median((Yt_pred2[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
      
      krig1_rmse = sqrt(mean((Yt_krig1[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
      krig2_rmse = sqrt(mean((Yt_krig2[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
      #krig3_rmse = sqrt(mean((Yt_krig3[pred_inds, t-begin_ind+1] - Ymat[pred_inds,Time - t + 1])^2))
      
      muY = mean(Ymat[pred_inds, ind2])
      sse = sqrt(mean((Ymat[pred_inds, ind2]- muY)^2))
      
      cat(#"BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
        "Time point: ", t, 
        " Pred RMSE: ", pred1_rmse, pred2_rmse, "\n",
        "Krig RMSE R2: ", 1-(krig1_rmse/sse)^2, 1-(krig2_rmse/sse)^2, "\n",
        "Pred RMSE R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
      
    }
    
    pred1_rmse = sqrt(mean((Yt_pred1[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    pred2_rmse = sqrt(mean((Yt_pred2[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    
    krig1_rmse = sqrt(mean((Yt_krig1[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    krig2_rmse = sqrt(mean((Yt_krig2[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    #krig3_rmse = sqrt(mean((Yt_krig3[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    
    muY = mean(Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])
    sse = sqrt(mean((Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)]- muY)^2))
    
    cat("BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
        "Krig RMSE: ", krig1_rmse, krig2_rmse, "\n",
        "Krig RMSE R2: ", 1-(krig1_rmse/sse)^2, 1-(krig2_rmse/sse)^2,  "\n",
        "Pred RMSE R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
  }
}


dif_mat = Yt_krig2[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)]
sse_mat = scale(Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)], center = T, scale = F)
sd(1 - colMeans(dif_mat^2)/colMeans(sse_mat^2))


sd(sqrt(colMeans(dif_mat^2)))