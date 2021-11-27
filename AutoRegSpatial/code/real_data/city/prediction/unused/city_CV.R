
### 本文件为比较Type I and II kriging的结果

setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/")

### source the functions
#source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/krigingMuFunc.R") # kriging functions 
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/city/realFunc_city2.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/exp/est_infer_func.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/predictFunc2.R") # prediction functions (similar to kriging functions)

load("rda/loc1.rda") # 读入经纬度数据
load("rda/Ymat.rda")

### date cut
begin_date = as.Date("2015-01-01") 
end_date = as.Date("2015-12-31") 
dd = as.Date(begin_date:end_date, origin = "1970-01-01")
Ymat0 = Ymat
gr = cut(dd, as.Date(c("2015-01-01", "2015-04-01", "2015-07-01", "2015-10-01", "2015-12-31")),include.lowest = T)
gr_num = split(1:length(dd), gr)

###
#Ymat = Ymat[,(ncol(Ymat)-(91)):ncol(Ymat)]


pred_inds = grep("(北京)|(天津)|(上海)|(南京)|(成都)|(广州)|(深圳)", loc1$cityname)
# set.seed(1234)
# pred_inds = sort(sample(1:nrow(Ymat), 50))


loc = loc1
loc1 = loc[,3:4]
dist_loc2 = as.matrix(dist(loc1))
cov.type = "exp"
if (cov.type == "quad")
  dist_loc2 = dist_loc2



est_win = 60
begin_ind = est_win + 1
pred_win = ncol(Ymat) - est_win
end_ind = begin_ind + pred_win - 1

set.seed(2345)
Yt_pred1 = matrix(0, nrow = nrow(Ymat), ncol = ncol(Ymat))
Yt_pred2 = matrix(0, nrow = nrow(Ymat), ncol = ncol(Ymat))


Time = ncol(Ymat)
for (t in begin_ind:end_ind)
{
  time_win = (t-est_win):t  #(Time - t - est_win + 1):(Time - t + 1)
  Ymat_win = Ymat[, time_win]
  
  cat(t, " | ")
  pred_ind = sample(1:nrow(loc1), 50)#pred_inds
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
  
  Yt_pred2[pred_ind, t] = predict.typeII(Y_train, rhos_est, mus_est, dist_train, dist_pred, 
                                              thetaE, thetaX, thetaMu, cov.type = cov.type)
  Yt_pred1[pred_ind, t] = predict.typeI(Y_train, rhos_est, mus_est, loc_train = loc1[train_ind,], loc_pred = loc1[pred_ind,],
                                             thetaE, thetaX, thetaMu, cov.type = cov.type)
  
  ind1 = begin_ind:t
  pred1_rmse = sqrt(median((Yt_pred1[pred_inds, ind1] - Ymat[pred_inds, ind1])^2))
  pred2_rmse = sqrt(median((Yt_pred2[pred_inds, ind1] - Ymat[pred_inds, ind1])^2))
  
  
  muY = mean(Ymat[pred_inds, ind1])
  sse = sqrt(mean((Ymat[pred_inds, ind1]- muY)^2))
  
  cat("Time point: ", t, 
      "Pred RMSE: ", pred1_rmse, pred2_rmse, "\n",
      "Pred RMSE R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
  
  pred1_rmse = sqrt(mean((Yt_pred1[pred_inds, ind1] - Ymat[pred_inds, ind1])^2))
  pred2_rmse = sqrt(mean((Yt_pred2[pred_inds, ind1] - Ymat[pred_inds, ind1])^2))
  
  muY = mean(Ymat[pred_inds, ind1])
  sse = sqrt(mean((Ymat[pred_inds, ind1]- muY)^2))
  
  cat("Pred Mean R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
  
}

pred1_rmse = sqrt(mean((Yt_pred1[pred_inds, ind1] - Ymat[pred_inds, ind1])^2))
pred2_rmse = sqrt(mean((Yt_pred2[pred_inds, ind1] - Ymat[pred_inds, ind1])^2))

muY = mean(Ymat[pred_inds, ind1])
sse = sqrt(mean((Ymat[pred_inds, ind1]- muY)^2))

cat("BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
    #"Krig RMSE: ", krig1_rmse, krig2_rmse, "\n",
    #"Krig RMSE R2: ", 1-(krig1_rmse/sse)^2, 1-(krig2_rmse/sse)^2,  "\n",
    "Pred RMSE R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")

