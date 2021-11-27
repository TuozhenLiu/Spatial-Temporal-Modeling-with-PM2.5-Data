### This file is for:
### fix sigma to estimate spatial dependence para & mean para
rm(list=ls())
setwd("E:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/exp/")

### format specify
specify_decimal = function(x, k) gsub('\\s+','',format(round(x, k), nsmall=k))                                   ### format function for keeping 2 digits after


#source("../../real_data/krigingMuFunc.R")
source("../../real_data/city/realFunc_city2.R")
source("../../real_data/predictFunc2.R")
source("infer_thetaE_func20170311.R")

Ns = c(200, 300)
Times = c(300, 500)


### show spatial correlation levels
set.seed(1234)
loc = simu.loc(100)
dist_loc2 = as.matrix(dist(loc))
para_cov = cov.exp(dist_loc2, 2,1)
hist(para_cov[upper.tri(para_cov, diag = F)])

betaE = 1.5; sig2E = 1
betaX = 1.5; sig2X = 1; alpha = 0.8
betaMu = 2; sig2Mu = 1; mu = 1


thetaX_true = c(alpha, betaX, sig2X)
theta = c(alpha, betaX, sig2X, betaE, sig2E, mu, betaMu, sig2Mu)

Nrep = 200
N = 100
set.seed(2345)
estRes0 = rep(list(), length(Ns))


i = 1; t = 1; r = 1
for (i in 1:length(Ns))
{
  estRes = rep(list(), length(Times))
  for (t in 1:length(Times))
  {
    N = Ns[i]; Time = Times[t]
    
    ### for the estimated values
    theta_est = matrix(0, nrow = 8, ncol = Nrep)
    thetaECI = matrix(T, nrow = 2, ncol = Nrep); thetaE_sd = matrix(0, nrow = 2, ncol = Nrep)
    thetaXCI = matrix(T, nrow = 3, ncol = Nrep); thetaX_sd = matrix(0, nrow = 3, ncol = Nrep)
    thetaMuCI = matrix(T, nrow = 3, ncol = Nrep); thetaMu_sd = matrix(0, nrow = 3, ncol = Nrep)
    
    ### for real true values (rhos, mus to be true)
    theta_true_est = matrix(0, nrow = 8, ncol = Nrep)
    thetaE_trueCI = matrix(T, nrow = 2, ncol = Nrep); thetaE_true_sd = matrix(0, nrow = 2, ncol = Nrep)
    thetaX_trueCI = matrix(T, nrow = 3, ncol = Nrep); thetaX_true_sd = matrix(0, nrow = 3, ncol = Nrep)
    thetaMu_trueCI = matrix(T, nrow = 3, ncol = Nrep); thetaMu_true_sd = matrix(0, nrow = 3, ncol = Nrep)
    
    for (r in 1:Nrep)
    {
      loc = simu.loc(N, cov.type = "exp")
      dist_loc2 = as.matrix(dist(loc))
      rhos = simu.rho(beta = betaX, sig2 = sig2X, N, alpha = alpha, dist_loc2, cov.type = "exp")
      mus = simu.mu(beta = betaMu, sig2 = sig2Mu, N, mu = mu, dist_loc2, cov.type = "exp")
      Y = simu.Y(beta = betaE, sigy2 = sig2E, N = N, Time = Time, rhos = rhos, 
                 dist_loc2, mu = mus, cov.type = "exp")
      
      train_ind = 1:N #1:floor(2/3*N)
      dist_train = dist_loc2[train_ind, train_ind]
      Y_train = Y[train_ind,]
      rhos_train = rhos[train_ind]
      mus_train = mus[train_ind]
      
      ### estimate rho
      rhos_mu = estRhoMu(Y_train)
      mus_est = rhos_mu[1,]
      rhos_est = rhos_mu[2,]
      
      thetaE_list = lse.theta(Y_train, dist_train, rhos = rhos_est, mu = mus_est, 
                              theta_true = c(betaE, sig2E), cov.type = "exp")
      thetaE = thetaE_list[[1]]
      thetaECI[,r] = thetaE_list[[3]]
      thetaE_sd[,r] = (thetaE_list[[2]])
      
      
      
      thetaX_list = estThetaRho(rhos_est, dist_train, theta_true = c(alpha, betaX, sig2X), cov.type = "exp")
      thetaX = thetaX_list[[1]]
      thetaXCI[,r] = thetaX_list[[3]]
      thetaX_sd[,r] = (thetaX_list[[2]])
      
      thetaMu_list = estThetaMu(mus_est, dist_train, theta_true = c(mu, betaMu, sig2Mu), cov.type = "exp")
      thetaMu = thetaMu_list[[1]]
      thetaMuCI[,r] = thetaMu_list[[3]]
      thetaMu_sd[,r] = (thetaMu_list[[2]])
      
      
      
      ### true values
      
      thetaE_true_list = lse.theta(Y_train, dist_train, rhos = rhos_train, mu = mus_train, theta_true = c(betaE, sig2E), cov.type = "exp")
      thetaE_true = thetaE_true_list[[1]]
      thetaE_trueCI[,r] = thetaE_true_list[[3]]
      thetaE_true_sd[,r] = (thetaE_true_list[[2]])
      
      thetaX_true_list = estThetaRho(rhos_train, dist_train, theta_true = c(alpha, betaX, sig2X), cov.type = "exp")
      thetaX_true = thetaX_true_list[[1]]
      thetaX_trueCI[,r] = thetaX_true_list[[3]]
      thetaX_true_sd[,r] = (thetaX_true_list[[2]])
      
      thetaMu_true_list = estThetaMu(mus_train, dist_train, theta_true = c(mu, betaMu, sig2Mu), cov.type = "exp")
      thetaMu_true = thetaMu_true_list[[1]]
      thetaMu_trueCI[,r] = thetaMu_true_list[[3]]
      thetaMu_true_sd[,r] = (thetaMu_true_list[[2]])
      
      
      theta_est[,r] = c(thetaX, thetaE, thetaMu)
      theta_true_est[,r] = c(thetaX_true, thetaE_true, thetaMu_true)
      
      
      cat(r, "\r")
    }
    estRes[[t]] = list(theta_est[c(4,5,6,7,8,1,2,3),], rbind(thetaECI,thetaMuCI, thetaXCI))
    
    cat("\n N: ", N, " Time:", Time, "\n",
        "Estimation Bias: ", rowMeans(theta_est-theta), "\n",
        "Estimation Bias: ", sqrt(rowMeans((theta_est-theta)^2)), "\n",
        "Coverage Prob: ", rowMeans(thetaXCI), rowMeans(thetaECI), rowMeans(thetaMuCI), "\n",
        "Coverage Prob (True): ", rowMeans(thetaX_trueCI), rowMeans(thetaE_trueCI), rowMeans(thetaMu_trueCI),  "\n",
        "Estimation SD: ", sqrt(rowMeans((theta_est-rowMeans(theta_est))^2)), "\n",
        "Estimation SD (True): ", sqrt(rowMeans((theta_true_est-rowMeans(theta_true_est))^2)), "\n",
        "Difference (||Est-True||): ", sqrt(rowMeans((theta_est-theta_true_est)^2)), "\n",
        "||Est-True||/SE: ", 
        sqrt(rowMeans((theta_est-theta_true_est)^2))/sqrt(rowMeans((theta_true_est-rowMeans(theta_true_est))^2)), "\n", 
        "Theta Asym SD: ", rowMeans(thetaX_sd), rowMeans(thetaE_sd),rowMeans(thetaMu_sd),"\n",
        "Theta_true Asym SD: ", rowMeans(thetaX_true_sd), rowMeans(thetaE_true_sd),  rowMeans(thetaMu_true_sd),
        "\n")
  }
  estRes0[[i]] = estRes
}
