setwd("E:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/")

source("../real_data/krigingMuFunc.R")
source("../real_data/city/realFunc_city.R")
source("../real_data/predictFunc.R")
source("infer_thetaE_func20170217.R")

Ns = c(100, 200, 400)
Times = c(200, 400, 800)

mu = 2; betaE = 1; sig2E = 3
#betaX = 1.5; sig2X = 2; alpha = 1
betaX = 1; sig2X = 0.8; alpha = 0.1


#betaE = 1; sig2E = 2
#betaX = 2; sig2X = 2; alpha = 0
thetaX_true = c(alpha, betaX, sig2X)
theta = c(alpha, betaX, sig2X, betaE, sig2E, mu)

Nrep = 500
N = 100
Time = 100 #Times[2]
set.seed(2345)
estRes0 = rep(list(), 3)

for (i in 1:length(Ns))
{
  estRes = rep(list(), length(Times))
  for (t in 1:length(Times))
  {
    N = Ns[i]; Time = Times[t]
    res0 = matrix(0, nrow = 4, ncol = Nrep)
    res = matrix(0, nrow = 4, ncol = Nrep)
    res_pred = matrix(0, nrow = 4, ncol = Nrep)
    theta_est = matrix(0, nrow = 6, ncol = Nrep)
    thetaECI = matrix(T, nrow = 2, ncol = Nrep)
    thetaE_sd = matrix(0, nrow = 2, ncol = Nrep)
    
    thetaXCI = matrix(T, nrow = 2, ncol = Nrep)
    thetaX_sd = matrix(0, nrow = 2, ncol = Nrep)
    
    muCI = rep(T, Nrep)
    mu_sd = rep(0, Nrep)
    
    theta_true_est = matrix(0, nrow = 6, ncol = Nrep)
    thetaE_trueCI = matrix(T, nrow = 2, ncol = Nrep)
    thetaE_true_sd = matrix(0, nrow = 2, ncol = Nrep)
    
    thetaX_trueCI = matrix(T, nrow = 2, ncol = Nrep)
    thetaX_true_sd = matrix(0, nrow = 2, ncol = Nrep)
    
    mu_trueCI = rep(T, Nrep)
    mu_true_sd = rep(0, Nrep)
    
    #thetaX_true_est = matrix(0, nrow = 3, ncol = Nrep)
    for (r in 1:Nrep)
    {
      loc = simu.loc(N)
      dist_loc2 = as.matrix(dist(loc))
      rhos = simu.rho(beta1 = betaX, sig2 = sig2X, N, alpha = alpha, dist_loc2)
      Y = simu.Y(beta1 = betaE, sigy2 = sig2E, N = N, Time = Time, rhos = rhos, dist_loc2, mu = mu)
      
      
      train_ind = 1:N #1:floor(2/3*N)
      krig_ind = setdiff(1:N, train_ind)
      dist_train = dist_loc2[train_ind, train_ind]
      dist_krig = dist_loc2[train_ind, krig_ind]
      Y_train = Y[train_ind,]
      rhos_train = rhos[train_ind]
      
      rhos_mu_list = estRhoMuIter(Y_train)
      mu_est = rhos_mu_list[[1]]
      
      rhos_est = rhos_mu_list[[2]]
      
      thetaE_list = lse.theta(Y_train, dist_train, rhos = rhos_est, mu = mu_est, theta_true = c(betaE, sig2E))
      thetaE = thetaE_list[[1]]
      thetaECI[,r] = thetaE_list[[3]]
      thetaE_sd[,r] = sqrt(diag(thetaE_list[[2]]))
      
      mu_infer_res = infer.mu(mu_est, mu_true = mu, thetaE, dist_loc2 = dist_train, Time, Y_train, rhos_est)
      mu_sd[r] = sqrt(mu_infer_res[[1]])
      muCI[r] = mu_infer_res[[2]]
      
      thetaX_list = estThetaRho(rhos_est, dist_train, theta_true = c(betaX, sig2X))
      thetaX = thetaX_list[[1]]
      thetaXCI[,r] = thetaX_list[[3]]
      thetaX_sd[,r] = sqrt(diag(thetaX_list[[2]]))
      
      ### true rhos
      mu_true = sum(Y_train[,-1] - rhos_train*Y_train[,-Time])/(length(train_ind)*Time)
      mu_true_infer_res = infer.mu(mu_true, mu_true = mu, thetaE, dist_loc2 = dist_train, Time, Y_train, rhos_train)
      mu_true_sd[r] = sqrt(mu_true_infer_res[[1]])
      mu_trueCI[r] = mu_true_infer_res[[2]]
      
      thetaE_true_list = lse.theta(Y_train, dist_train, rhos = rhos_train, mu = mu, theta_true = c(betaE, sig2E))
      thetaE_true = thetaE_true_list[[1]]
      thetaE_trueCI[,r] = thetaE_true_list[[3]]
      thetaE_true_sd[,r] = sqrt(diag(thetaE_true_list[[2]]))
      
      thetaX_true_list = estThetaRho(rhos_train, dist_train, theta_true = c(betaX, sig2X))
      thetaX_true = thetaX_true_list[[1]]
      thetaX_trueCI[,r] = thetaX_true_list[[3]]
      thetaX_true_sd[,r] = sqrt(diag(thetaX_true_list[[2]]))
      
      #thetaX_true_est[,r] = estThetaRho(rhos_train, dist_train)
      theta_est[,r] = c(thetaX, thetaE, mu_est)
      theta_true_est[,r] = c(thetaX_true, thetaE_true, mu_true)
      
      
      cat(r, "\r")
    }
    #estRes[[t]] = list(res0, res, theta_est, thetaX_true_est)
    cat("\n N: ", N, " Time:", Time, "\n",
        "Estimation Bias: ", rowMeans(theta_est-theta), "\n",
        "Coverage Prob: ", rowMeans(thetaXCI), rowMeans(thetaECI), mean(muCI), "\n",
        "Coverage Prob (True): ", rowMeans(thetaX_trueCI), rowMeans(thetaE_trueCI), mean(mu_trueCI),  "\n",
        #"thetaXTrueRho Bias: ", rowMeans(thetaX_true_est-thetaX_true), "\n",
        "Estimation SD: ", sqrt(rowMeans((theta_est-rowMeans(theta_est))^2)), "\n",
        "Estimation SD (True): ", sqrt(rowMeans((theta_true_est-rowMeans(theta_true_est))^2)), "\n",
        "Difference (||Est-True||): ", sqrt(rowMeans((theta_est-theta_true_est)^2)), "\n",
        "||Est-True||/SE: ", 
        sqrt(rowMeans((theta_est-theta_true_est)^2))/sqrt(rowMeans((theta_true_est-rowMeans(theta_true_est))^2)), "\n", 
        "Theta Asym SD: ", rowMeans(thetaX_sd), rowMeans(thetaE_sd), mean(mu_sd),
        "Theta_true Asym SD: ", rowMeans(thetaX_true_sd), rowMeans(thetaE_true_sd), mean(mu_true_sd),
        #"thetaXTrueRho SD: ", sqrt(rowMeans((thetaX_true_est-rowMeans(thetaX_true_est))^2)), 
        "\n")
  }
  estRes0[[i]] = estRes
}
save(estRes0, file = "../data/estRes0.rda")


for (i in 1:length(Ns))
{
  cat(Ns[i], " & ", 2*Ns[i], " & ")
  rr = estRes[[i]]
  theta_est = rr[[3]]
  est = round((rowMeans(theta_est-theta))*10,2)
  SD = round(sqrt(rowMeans((theta_est-rowMeans(theta_est))^2))*10,2)
  cat(paste0(paste(est, SD, sep = " (", collapse = ") & "), ")"))
  kbias = round(rowMeans(rr[[1]][3:4,])*10, 2)
  ksd = round(sqrt(rowMeans(rr[[2]][3:4,])), 2)
  cat(" & ", paste0(paste(kbias, ksd, sep = " (", collapse = ") & "), ")"), "\\\\ \n")
}

