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

Ns = c(100, 200, 300)
Times = c(300, 500)


### for quad function
betaE = 1; sig2E = 1
betaX = 1; sig2X = 1; alpha = 0.8
betaMu = 1; sig2Mu = 1; mu = 1

### show spatial correlation levels
set.seed(1234)
loc = simu.loc(100, cov.type = "quad")
dist_loc2 = as.matrix(dist(loc))
para_cov = cov.quad(dist_loc2, 2,1)
hist(para_cov[upper.tri(para_cov, diag = F)])


thetaX_true = c(alpha, betaX, sig2X)
theta = c(alpha, betaX, sig2X, betaE, sig2E, mu, betaMu, sig2Mu)

Nrep = 300
N = 100
Time = 100 #Times[2]
set.seed(2345)
estRes0 = rep(list(), length(Ns))

for (i in 1:length(Ns))
{
  estRes = rep(list(), length(Times))
  for (t in 1:1)
  {
    N = Ns[i]; Time = Times[t]
    
    ### for the estimated values
    theta_est = matrix(0, nrow = 8, ncol = Nrep)
    thetaECI = matrix(T, nrow = 2, ncol = Nrep); thetaE_sd = matrix(0, nrow = 2, ncol = Nrep)
    thetaXCI = matrix(T, nrow = 2, ncol = Nrep); thetaX_sd = matrix(0, nrow = 2, ncol = Nrep)
    thetaMuCI = matrix(T, nrow = 2, ncol = Nrep); thetaMu_sd = matrix(0, nrow = 2, ncol = Nrep)
    muCI = rep(T, Nrep); mu_sd = rep(0, Nrep)
    alphaCI = rep(T, Nrep); alpha_sd = rep(0, Nrep)
    
    ### for real true values (rhos, mus to be true)
    theta_true_est = matrix(0, nrow = 8, ncol = Nrep)
    thetaE_trueCI = matrix(T, nrow = 2, ncol = Nrep); thetaE_true_sd = matrix(0, nrow = 2, ncol = Nrep)
    thetaX_trueCI = matrix(T, nrow = 2, ncol = Nrep); thetaX_true_sd = matrix(0, nrow = 2, ncol = Nrep)
    thetaMu_trueCI = matrix(T, nrow = 2, ncol = Nrep); thetaMu_true_sd = matrix(0, nrow = 2, ncol = Nrep)
    mu_trueCI = rep(T, Nrep); mu_true_sd = rep(0, Nrep)
    alpha_trueCI = rep(T, Nrep); alpha_true_sd = rep(0, Nrep)
    
    muCI = matrix(T, ncol = Nrep, nrow = 3)
    mu_est = matrix(0, ncol = Nrep, nrow = 3)
    mu_sd = matrix(0, ncol = Nrep, nrow = 3)
    
    
    for (r in 1:Nrep)
    {
      loc = simu.loc(N, cov.type = "quad")*2
      dist_loc2 = as.matrix(dist(loc))
      
      mus = simu.mu(beta = betaMu, sig2 = sig2Mu, N, mu = mu, dist_loc2)
      
      train_ind = 1:N #1:floor(2/3*N)
      dist_train = dist_loc2[train_ind, train_ind]
      mus_train = mus[train_ind]
      
      thetaMu_true_list = estThetaMu(mus, dist_loc2, infer = T, 
                                     theta_true = c(mu, betaMu, sig2Mu), cov.type = "quad")
      #thetaMu_true_list = lse.X(mus_train, dist_loc2 = dist_train, theta_true = c(mu, betaMu, sig2Mu))
      # estThetaMu(mus_train, dist_train, theta_true = c(betaMu, sig2Mu), 
      #                              cov.type = "quad")
      # thetaMu_true = thetaMu_true_list[[1]]
      # thetaMu_trueCI[,r] = thetaMu_true_list[[2]]
      # thetaMu_true_sd[,r] = (thetaMu_true_list[[3]])
      # 
      # mu_true_infer_res = infer.mu(thetaMu_true, mu_true = mu, dist_loc2, Time, 
      #                              cov.type = "quad")
      # mu_true_sd[r] = sqrt(mu_true_infer_res[[1]])
      # mu_trueCI[r] = mu_true_infer_res[[2]]
      # 
      # mu_est[,r] = c(thetaMu_true)
      # mu_sd[,r] = c(mu_true_infer_res[[1]], thetaMu_true_list[[2]])
      # muCI[,r] = c(thetaMu_true_list[[3]], mu_true_infer_res[[2]])
      
      mu_est[,r] = thetaMu_true_list[[1]]
      mu_sd[,r] = thetaMu_true_list[[2]]
      muCI[,r] = (thetaMu_true_list[[3]])
      
      cat(r, "\r")
    }
    cat("N: ", N,"\n",
        "Bias: ", rowMeans(mu_est - c(mu, betaMu, sig2Mu)), 
        "SD: ", apply(mu_est,1, sd),
        "Asym SD: ", rowMeans(mu_sd), 
        "CP: ", rowMeans(muCI), "\n")
  }
}





