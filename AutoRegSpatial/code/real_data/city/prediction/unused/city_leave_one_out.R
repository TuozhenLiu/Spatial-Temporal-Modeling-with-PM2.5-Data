
rm(list = ls())


### 本文件为比较Type I and II kriging的结果

setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/")

### source the functions
#source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/krigingMuFunc.R") # kriging functions 
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/city/realFunc_city2.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/exp/est_infer_func.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/predictFunc2.R") # prediction functions (similar to kriging functions)

load("rda/loc1.rda") # 读入经纬度数据
load("rda/Ymat.rda")
# estimate the parameters in rhos (theta_X)
lse.X<-function(Y, dist_loc2, infer = T, theta_true, cov.type = "exp", p.value = F) # Y here is X-alpha, which is centered
{
  Time = ncol(Y)
  if(is.null(Time))
    Time = 1
  if (is.null(dim(Y)))
    sig_hat = tcrossprod(Y)
  else
    sig_hat = tcrossprod(Y)/ncol(Y)
  
  #sig_hat = tcrossprod(Y) 
  sig2_hat = mean(Y^2)#var(as.vector(Y))
  
  theta = c(0.5)
  iter = 1; del = 1
  if (cov.type == "quad")
  {
    lse.step = lse.step.quad
    infer.cov = infer.cov.quad
  }
  if (infer&p.value)
    theta_true[2] = sig2_hat
  
  
  while(mean(abs(del))>10^-3&iter<1000)
  {
    #cat(mean(abs(del)), " ", theta, "\n")
    del = lse.step(sig_hat, beta = theta, sigy2 = theta_true[2], dist_loc2)
    theta = theta - del*0.5
    if (any(theta<0|theta>10))
      theta = runif(1, 0.01,10)
    iter = iter+1
  }
  theta = abs(theta)
  if (!infer)
    return(c(theta, sig2_hat))
  
  ### for spatial parameter inference
  est_cov = infer.cov(sig_hat, beta = theta, sigy2 = theta_true[2], dist_loc2, Time = Time)
  CI = (theta-1.96*sqrt((est_cov))<theta_true[1])&(theta_true[1]<theta+1.96*sqrt((est_cov)))
  
  ### for sigma inference
  
  if (cov.type == "quad")
    sig = cov.quad(dist_loc2, theta, sig2_hat)
  else
    sig = cov.exp(dist_loc2, theta, sig2_hat)
  
  sig2_cov = 2*tr(sig%*%sig)/nrow(dist_loc2)^2/Time
  sig2_CI = (sig2_hat-1.96*sqrt((sig2_cov))<theta_true[2])&(theta_true[2]<sig2_hat+1.96*sqrt((sig2_cov)))
  
  return(list(c(abs(theta), sig2_hat), c(sqrt(est_cov), sqrt(sig2_cov)), c(CI, sig2_CI)))
  
  #return(list(abs(theta), est_cov, CI)) # since exponential covariance model is symmetric of beta, we restrict theta to be positive
}




# ii = loc1[,3]>100
# loc1 = loc1[ii,]
# Ymat = Ymat[ii,]

#load("rda/rhos_df.rda")
set.seed(1234)
# city_ind = (1:nrow(rhos_df))
# Ymat = Ymat[city_ind,]
# rhos_df = rhos_df[city_ind,]

loc = loc1
#pred_inds = grep("(北京)|(天津)|(上海)|(南京)|(成都)|(广州)|(深圳)", loc1$cityname)

# station_rho_loc = merge(rhos_df,loc, by = "citycode")
# Ymat = Ymat[match(station_rho_loc$citycode, rhos_df$citycode),]

cov.type = "quad"
loc1 = loc[,3:4]*4
dist_loc2 = as.matrix(dist(loc1))
# if (cov.type == 'quad')
#   dist_loc2 = dist_loc2*4
#Ymat = Ymat0[,ind1]


begin_ind = 1
pred_win = 305
set.seed(2345)
set.seed(12)
pred_inds = sort(sample(1:nrow(loc1), 50)) #1:nrow(loc1)#sort(sample(1:nrow(loc1), 50))#[-c(9,49)]


for (begin_ind in c(1))
{
  end_ind = begin_ind + pred_win -1
  for (est_win in c(60))
  {
    
    # Yt_krig1 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    # Yt_krig2 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    # Yt_krig3 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    Yt_pred1 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    Yt_pred2 = matrix(0, nrow = nrow(Ymat), ncol = pred_win)
    
    Time = ncol(Ymat)
    for (t in begin_ind:end_ind)
    {
      time_win = (Time - t - est_win + 1):(Time - t + 1)
      Ymat_win = Ymat[, time_win]
      
      for (i in pred_inds)
      {
        cat(t, i,  "\r")
        train_ind = setdiff(1:nrow(Ymat), i)
        
        dist_train = dist_loc2[train_ind, train_ind]
        dist_pred = dist_loc2[train_ind, i]
        Y_train = Ymat_win[train_ind,]
        
        rhosMu = estRhoMu(Y_train)
        mus_est = rhosMu[1,]
        rhos_est = rhosMu[2,]
        
        thetaMu = estThetaMu(mus_est, dist_train, infer = F,  cov.type = cov.type)
        (thetaX = estThetaRho(rhos_est, dist_train, infer = F, cov.type = cov.type))
        (thetaE = lse.theta(Y_train, dist_train, rhos = rhos_est, 
                            mu = mus_est, infer = F, cov.type = cov.type))
        
        
        Yt_pred2[i, t-begin_ind+1] = predict.typeII(Y_train, rhos_est, mus_est, dist_train, dist_pred, 
                                                    thetaE, thetaX, thetaMu, cov.type = cov.type)
        Yt_pred1[i, t-begin_ind+1] = predict.typeI(Y_train, rhos_est, mus_est, loc_train = loc1[train_ind,], loc_pred = loc1[i,],
                                                   thetaE, thetaX, thetaMu, cov.type = cov.type)
      }
      ind1 = 1:(t-begin_ind+1)
      ind2 = (Time -begin_ind+1):(Time - t + 1)
      pred1_rmse = sqrt(median((Yt_pred1[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
      pred2_rmse = sqrt(median((Yt_pred2[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
      
      muY = mean(Ymat[pred_inds, ind2])
      sse = sqrt(mean((Ymat[pred_inds, ind2]- muY)^2))
      
      cat(#"BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
        "Time point: ", t, 
        " Pred RMSE: ", pred1_rmse, pred2_rmse, "\n",
        #"Krig RMSE R2: ", 1-(krig1_rmse/sse)^2, 1-(krig2_rmse/sse)^2, "\n",
        "Pred RMSE R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
      
      pred1_rmse = sqrt(mean((Yt_pred1[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
      pred2_rmse = sqrt(mean((Yt_pred2[pred_inds, ind1] - Ymat[pred_inds, ind2])^2))
      
      
      muY = mean(Ymat[pred_inds, ind2])
      sse = sqrt(mean((Ymat[pred_inds, ind2]- muY)^2))
      
      cat("Pred Mean R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
      
    }
    
    pred1_rmse = sqrt(mean((Yt_pred1[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    pred2_rmse = sqrt(mean((Yt_pred2[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    
    # krig1_rmse = sqrt(mean((Yt_krig1[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    # krig2_rmse = sqrt(mean((Yt_krig2[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    # #krig3_rmse = sqrt(mean((Yt_krig3[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])^2))
    
    muY = mean(Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)])
    sse = sqrt(mean((Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)]- muY)^2))
    
    cat("BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
        #"Krig RMSE: ", krig1_rmse, krig2_rmse, "\n",
        #"Krig RMSE R2: ", 1-(krig1_rmse/sse)^2, 1-(krig2_rmse/sse)^2,  "\n",
        "Pred RMSE R2: ", 1- (pred1_rmse/sse)^2, 1- (pred2_rmse/sse)^2, "\n")
    
    
  }
}

# save(Yt_pred1, file = "E:/OneDrive/1. Academic/spatial_autoreg/data/realdata/Yt_pred1_quad.rda")
# 
# save(Yt_pred2, file = "E:/OneDrive/1. Academic/spatial_autoreg/data/realdata/Yt_pred2_quad.rda")
