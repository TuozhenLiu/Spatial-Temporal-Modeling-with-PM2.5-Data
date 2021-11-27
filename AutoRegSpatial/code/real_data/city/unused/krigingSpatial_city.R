
### 本文件用于比较单纯用spatial的结果(相当于T = 1)

lse.theta0<-function(Y, dist_loc2)
{
  Y = Y - mean(Y)
  sig_hat = tcrossprod(Y)
  theta = lse.X(Y, dist_loc2, infer = F, theta_true = c(1,1))
  return(theta)
  
  theta = c(1, var(as.vector(Y)))
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


setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/")

### source the functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/unused/krigingMuFunc.R") # kriging functions 
source("E:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/exp/est_infer_func.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/city/realFunc_city2.R") # real data functions
#source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/predictFunc.R") # prediction functions (similar to kriging functions)

load("rda/loc1.rda") # 读入经纬度数据
load("rda/Ymat.rda")
load("rda/rhos_df.rda")

#(thetaY = lse.theta0(Ymat[,ncol(Ymat)], dist_loc2))


# station_rho_loc = merge(rhos_df,loc, by = "citycode")
# Ymat = Ymat[match(station_rho_loc$citycode, rhos_df$citycode),]
# loc1 = station_rho_loc[,4:5]
dist_loc2 = as.matrix(dist(loc1[,3:4]))

begin_ind = 1
set.seed(2345)
#pred_inds = sort(sample(1:nrow(loc1), 30))

pred_inds = grep("(北京)|(天津)|(上海)|(南京)|(保定)|(石家庄)|(成都)|(广州)|(深圳)", loc1$cityname)

set.seed(1234)

Yt_krigs = matrix(0, nrow = nrow(Ymat), ncol = length(begin_ind:end_ind))
Time = ncol(Ymat)
for (t in begin_ind:end_ind)
{
  for (i in pred_inds)
  {
    cat(t, i, "\r")
    
    Y_train = Ymat[-i,t]
    
    (thetaY = lse.theta0(Y_train, dist_loc2[-i,-i]))
    
    Yt_krigs[i, t-begin_ind+1] = krigingS(Y_train, 
                                          loc_train = loc1[-i,3:4], 
                                          loc_pred = loc1[i,3:4], 
                                          thetaY)
    
  }
  
  krigs_rmse = sqrt(mean((Yt_krigs[pred_inds, t] - Ymat[pred_inds,t])^2))
  
  muY = mean(Ymat[pred_inds,t])
  sse = sqrt(mean((Ymat[pred_inds, t]- muY)^2))
  
  cat("BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
      "Krig RMSE: ", krigs_rmse,"\n",
      "Krig RMSE R2: ", 1-(krigs_rmse/sse)^2, "\n")
}

krigs_rmse = sqrt(mean((Yt_krigs[pred_inds,] - 
                          Ymat[pred_inds, begin_ind:end_ind])^2))
muY = mean(Ymat[pred_inds,begin_ind:end_ind])
sse = sqrt(mean((Ymat[pred_inds,begin_ind:end_ind]- muY)^2))

cat("BackPred Start Point: ", begin_ind, "Est_win: ", est_win, "Pred_win: ", pred_win, "\n",
    "Krig RMSE: ", krigs_rmse,"\n",
    "Krig RMSE R2: ", 1-(krigs_rmse/sse)^2, "\n")


resi = Ymat[pred_inds, begin_ind:end_ind] - Yt_krigs[pred_inds,]


#(北京)|(天津)|(上海)|(南京)|(保定)|(石家庄)|(成都)|(广州)|(深圳)
cities = as.character(loc1$cityname[pred_inds]) #c('北京', '天津', '上海', '南京', '保定', '石家庄', '成都', '广州', '深圳')
par(mfrow = c(1,2))

for (i in 1:9)
{
  plot(Yt_krigs[pred_inds,][i,], resi[i,], xlab = "Kriging 值", 
       ylab = "Kriging Residuals", main = cities[i])
  
  cat(cities[i], " : ",  mean(resi[i,]),  # residual的均值
      mean(resi[i,]<0), # 预报值低于kriging值的比例
      1-mean(resi[i,]^2)/var(Ymat[pred_inds, ][i,]), "\n ") 
      
}
save(Yt_krigs, file = "rda/Yt_krigs.rda")







dif_mat = Yt_krigs[pred_inds,] - Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)]
sse_mat = scale(Ymat[pred_inds,(Time -begin_ind+1):(Time -begin_ind+1 - pred_win +1)], center = T, scale = F)
sd(1 - colMeans(dif_mat^2)/colMeans(sse_mat^2))



near_city = apply(dist_loc2, 1, function(x){
  order(x, decreasing = F)[2]
})

Ymat_near = Ymat[near_city,]
sqrt(mean((Ymat - Ymat_near)^2))

plot(exp(Ymat_near[,3]), exp(Ymat[,3]))

