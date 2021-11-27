


### 本文件为比较Type I and II predictive kriging的结果


setwd("E:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/exp/")

### source the functions
#source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/krigingMuFunc.R") # kriging functions 
source("../../real_data/city/realFunc_city2.R") # real data functions
source("../../real_data/predictFunc2.R") # prediction functions (similar to kriging functions)
source("infer_thetaE_func20170311.R")
#source("E:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/infer_thetaE_func20170226.R")

######################################################################################################

set.seed(1234)

Ns = c(200, 300)
Times = c(300, 500)

### strong spatial correlation
betaE = 1; sig2E = 2
betaX = 1; sig2X = 1.2; alpha = 0.8
betaMu = 1; sig2Mu = 1.5; mu = 1


# ### weak spatial correlation
# betaE = 2; sig2E = 2
# betaX = 1.5; sig2X = 1.2; alpha = 0.8
# betaMu = 2; sig2Mu = 1.5; mu = 1
# 

set.seed(2345)
Nrep = 200
begin_ind = 1 # 从最后一期往前数1期
pred_win = 1

res_rmse_all = rep(list(), length(Ns))
res_r2_all = rep(list(), length(Ns))

#N  = 200; Time = 200
for (k in 1:length(Ns))
{
  N = Ns[k]
  res_rmse = matrix(0, nrow = length(Times), 2)
  res_r2 = matrix(0, nrow = length(Times), 2)
  
  for (tt in 1:length(Times))
  {
    Time = Times[tt]
    est_win = Time - pred_win #c(Time - pred_win)
    n_pred = floor(N/5)
    
    cat("N: ", N, " Time: ", Time, "\n")
    perf_r2 = matrix(0, nrow = Nrep, ncol = 2)
    perf_rmse = matrix(0, nrow = Nrep, ncol = 2)
    
    for (r in 1:Nrep)
    {
      loc1 = data.frame(simu.loc(N))
      dist_loc2 = as.matrix(dist(loc1))
      rhos = simu.rho(beta1 = betaX, sig2 = sig2X, N, alpha = alpha, dist_loc2)
      mus = simu.mu(beta1 = betaMu, sig2 = sig2Mu, N, mu = mu, dist_loc2)
      Ymat = simu.Y(beta1 = betaE, sigy2 = sig2E, N = N, Time = Time, rhos = rhos, dist_loc2, mu = mus)
      
      pred_inds = sort(sample(1:nrow(loc1), n_pred))
      
      for (begin_ind in c(1))
      {
        end_ind = begin_ind + pred_win -1
        
        Yt_krig1 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
        Yt_krig2 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
        Yt_pred1 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
        Yt_pred2 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
        
        time_win = (Time - end_ind - est_win + 1):(Time - end_ind + 1)
        Ymat_win = Ymat[, time_win]
        
        pred_ind = pred_inds
        train_ind = setdiff(1:nrow(Ymat), pred_ind)
        
        dist_train = dist_loc2[train_ind, train_ind]
        dist_pred = dist_loc2[train_ind, pred_ind]
        Y_train = Ymat_win[train_ind,]
        
        rhos_mu = estRhoMu(Y_train)
        mus_est = rhos_mu[1,]
        rhos_est = rhos_mu[2,]
        
        thetaMu = estThetaMu(mus_est, dist_train, infer = F)
        (thetaX = estThetaRho(rhos_est, dist_train, infer = F))
        (thetaE = lse.theta(Y_train, dist_train, rhos = rhos_est, mu = mus_est, infer = F))
        
        
        for (t in begin_ind:end_ind)
        {
          
          Yt_pred2[pred_ind, t-begin_ind+1] = predict.typeII(Y_train, rhos_est, mus_est, dist_train, dist_pred, 
                                                             thetaE, thetaX, thetaMu)
          Yt_pred1[pred_ind, t-begin_ind+1] = predict.typeI(Y_train, rhos_est, mus_est, loc_train = loc1[train_ind,], loc_pred = loc1[pred_ind,],
                                                            thetaE, thetaX, thetaMu)
          
          ind1 = 1:(t-begin_ind+1)
          ind2 = (Time -begin_ind+1):(Time - t + 1)
          pred1_rmse = sqrt(mean((Yt_pred1[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
          pred2_rmse = sqrt(mean((Yt_pred2[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
          
          muY = mean(Ymat[pred_inds, ind2])
          sse = sqrt(mean((Ymat[pred_inds, ind2]- muY)^2))
          #cat(t, " | ")
        }
        
        pred1_rmse = sqrt(median((Yt_pred1[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
        pred2_rmse = sqrt(median((Yt_pred2[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
        
        muY = mean(Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])
        sse = sqrt(mean((Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)]- muY)^2))
      }
      perf_r2[r,] = c(1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2)
      perf_rmse[r,] = c(pred1_rmse, pred2_rmse)
      cat(r, perf_r2[r,], perf_rmse[r,], "\r")
    }
    cat("\n N: ", N, " Time: ", Time, "\n",
        #"Krig RMSE R2: ", colMeans(perf_r2[,1:2]),  "\n",
        "Pred RMSE R2: ", colMeans(perf_r2, na.rm = T), "\n")
    res_r2[tt,] = (colMeans(perf_r2, na.rm = T))
    res_rmse[tt,] = (colMeans(perf_rmse, na.rm = T))
  }
  res_r2_all[[k]] = res_r2
  res_rmse_all[[k]] = res_rmse
}



# save(res_r2_all, file = "../../../data/simu/res_r2_all_weak.rda")
# save(res_rmse_all, file = "../../../data/simu/res_rmse_all_weak.rda")


save(res_r2_all, file = "../../../data/simu/res_r2_all_strong.rda")
save(res_rmse_all, file = "../../../data/simu/res_rmse_all_strong.rda")





for (k in 1:2)
{
  for (t in 1:2)
  {
    cat((res_r2_all[[k]][t,]), "\n")
  }
}









