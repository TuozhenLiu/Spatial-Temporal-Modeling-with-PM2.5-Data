setwd("E:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/")

source("../real_data/krigingMuFunc.R")
source("../real_data/city/realFunc_city_wave.R")
#source("../real_data/predictFunc2.R")
source("infer_wave_func20170308.R")

Ns = c(100, 200, 400)
Times = 5#c(200, 400, 800)

phiE = 1; sig2E = 3
phiX = 2; sig2X = 1; alpha = 0
phiMu = 3; sig2Mu = 1; mu = 1

#phiE = 1; sig2E = 2
#phiX = 2; sig2X = 2; alpha = 0
thetaX_true = c(alpha, phiX, sig2X)
theta = c(alpha, phiX, sig2X, phiE, sig2E, mu, phiMu, sig2Mu)

Nrep = 100
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
    
    theta_est = matrix(0, nrow = 8, ncol = Nrep)
    thetaECI = matrix(T, nrow = 2, ncol = Nrep)
    thetaE_sd = matrix(0, nrow = 2, ncol = Nrep)
    
    thetaXCI = matrix(T, nrow = 2, ncol = Nrep)
    thetaX_sd = matrix(0, nrow = 2, ncol = Nrep)
    
    thetaMuCI = matrix(T, nrow = 2, ncol = Nrep)
    thetaMu_sd = matrix(0, nrow = 2, ncol = Nrep)
    
    muCI = rep(T, Nrep)
    mu_sd = rep(0, Nrep)
    
    alphaCI = rep(T, Nrep)
    alpha_sd = rep(0, Nrep)
    
    theta_true_est = matrix(0, nrow = 8, ncol = Nrep)
    thetaE_trueCI = matrix(T, nrow = 2, ncol = Nrep)
    thetaE_true_sd = matrix(0, nrow = 2, ncol = Nrep)
    
    thetaX_trueCI = matrix(T, nrow = 2, ncol = Nrep)
    thetaX_true_sd = matrix(0, nrow = 2, ncol = Nrep)
    
    thetaMu_trueCI = matrix(T, nrow = 2, ncol = Nrep)
    thetaMu_true_sd = matrix(0, nrow = 2, ncol = Nrep)
    
    mu_trueCI = rep(T, Nrep)
    mu_true_sd = rep(0, Nrep)
    
    alpha_trueCI = rep(T, Nrep)
    alpha_true_sd = rep(0, Nrep)
    
    #thetaX_true_est = matrix(0, nrow = 3, ncol = Nrep)
    for (r in 1:Nrep)
    {
      loc = simu.loc(N)
      dist_loc2 = as.matrix(dist(loc))
      rhos = simu.rho(phi = phiX, sig2 = sig2X, N, alpha = alpha, dist_loc2)
      mus = simu.mu(phi = phiMu, sig2 = sig2Mu, N, mu = mu, dist_loc2)
      Y = simu.Y(phi = phiE, sigy2 = sig2E, N = N, Time = Time, rhos = rhos, dist_loc2, mu = mus)
      
      
      train_ind = 1:N #1:floor(2/3*N)
      dist_train = dist_loc2[train_ind, train_ind]
      #dist_krig = dist_loc2[train_ind, krig_ind]
      Y_train = Y[train_ind,]
      rhos_train = rhos[train_ind]
      mus_train = mus[train_ind]
      
      rhos_mu = estRhoMu(Y_train)
      mus_est = rhos_mu[1,]
      rhos_est = rhos_mu[2,]
      
      thetaE_list = lse.theta(Y_train, dist_train, rhos = rhos_est, mu = mus_est, theta_true = c(phiE, sig2E))
      thetaE = thetaE_list[[1]]
      (thetaECI[,r] = thetaE_list[[3]])
      thetaE_sd[,r] = sqrt(diag(thetaE_list[[2]]))
      
      
      
      (thetaX_list = estThetaRho(rhos_est, dist_train, theta_true = c(phiX, sig2X)))
      (thetaX = thetaX_list[[1]])
      thetaXCI[,r] = thetaX_list[[3]]
      thetaX_sd[,r] = sqrt(diag(thetaX_list[[2]]))
      
      (thetaMu_list = estThetaMu(mus_est, dist_train, theta_true = c(phiMu, sig2Mu)))
      thetaMu = thetaMu_list[[1]]
      thetaMuCI[,r] = thetaMu_list[[3]]
      thetaMu_sd[,r] = sqrt(diag(thetaMu_list[[2]]))
      
      mu_infer_res = infer.mu(thetaMu, mu_true = mu, dist_loc2, Time)
      mu_sd[r] = sqrt(mu_infer_res[[1]])
      muCI[r] = mu_infer_res[[2]]
      
      alpha_infer_res = infer.alpha(thetaX, alpha_true = alpha, dist_loc2, Time)
      alpha_sd[r] = sqrt(alpha_infer_res[[1]])
      alphaCI[r] = alpha_infer_res[[2]]
      
      
      ### true values
      
      thetaE_true_list = lse.theta(Y_train, dist_train, rhos = rhos_train, mu = mus_train, theta_true = c(phiE, sig2E))
      thetaE_true = thetaE_true_list[[1]]
      thetaE_trueCI[,r] = thetaE_true_list[[3]]
      thetaE_true_sd[,r] = sqrt(diag(thetaE_true_list[[2]]))
      
      (thetaX_true_list = estThetaRho(rhos_train, dist_train, theta_true = c(phiX, sig2X)))
      thetaX_true = thetaX_true_list[[1]]
      thetaX_trueCI[,r] = thetaX_true_list[[3]]
      thetaX_true_sd[,r] = sqrt(diag(thetaX_true_list[[2]]))
      
      thetaMu_true_list = estThetaMu(mus_train, dist_train, theta_true = c(phiMu, sig2Mu))
      thetaMu_true = thetaMu_true_list[[1]]
      thetaMu_trueCI[,r] = thetaMu_true_list[[3]]
      thetaMu_true_sd[,r] = sqrt(diag(thetaMu_true_list[[2]]))
      
      mu_true_infer_res = infer.mu(thetaMu_true, mu_true = mu, dist_loc2, Time)
      mu_true_sd[r] = sqrt(mu_true_infer_res[[1]])
      mu_trueCI[r] = mu_true_infer_res[[2]]
      
      alpha_true_infer_res = infer.alpha(thetaX_true, alpha_true = alpha, dist_loc2, Time)
      alpha_true_sd[r] = sqrt(alpha_true_infer_res[[1]])
      alpha_trueCI[r] = alpha_true_infer_res[[2]]
      
      #thetaX_true_est[,r] = estThetaRho(rhos_train, dist_train)
      theta_est[,r] = c(thetaX, thetaE, thetaMu)
      theta_true_est[,r] = c(thetaX_true, thetaE_true, thetaMu_true)
      
      
      cat(r, "\r")
    }
    #estRes[[t]] = list(res0, res, theta_est, thetaX_true_est)
    cat("\n N: ", N, " Time:", Time, "\n",
        "Estimation Bias: ", rowMeans(theta_est-theta), "\n",
        "Coverage Prob: ", mean(alphaCI), rowMeans(thetaXCI), rowMeans(thetaECI), mean(muCI), rowMeans(thetaMuCI), "\n",
        "Coverage Prob (True): ", mean(alpha_trueCI), rowMeans(thetaX_trueCI), rowMeans(thetaE_trueCI), mean(mu_trueCI), rowMeans(thetaMu_trueCI),  "\n",
        #"thetaXTrueRho Bias: ", rowMeans(thetaX_true_est-thetaX_true), "\n",
        "Estimation SD: ", sqrt(rowMeans((theta_est-rowMeans(theta_est))^2)), "\n",
        "Estimation SD (True): ", sqrt(rowMeans((theta_true_est-rowMeans(theta_true_est))^2)), "\n",
        "Difference (||Est-True||): ", sqrt(rowMeans((theta_est-theta_true_est)^2)), "\n",
        "||Est-True||/SE: ", 
        sqrt(rowMeans((theta_est-theta_true_est)^2))/sqrt(rowMeans((theta_true_est-rowMeans(theta_true_est))^2)), "\n", 
        "Theta Asym SD: ", mean(alpha_sd), rowMeans(thetaX_sd), rowMeans(thetaE_sd), mean(mu_sd), rowMeans(thetaMu_sd),"\n",
        "Theta_true Asym SD: ", mean(alpha_true_sd), rowMeans(thetaX_true_sd), rowMeans(thetaE_true_sd), mean(mu_true_sd), rowMeans(thetaMu_true_sd),
        #"thetaXTrueRho SD: ", sqrt(rowMeans((thetaX_true_est-rowMeans(thetaX_true_est))^2)), 
        "\n")
  }
  estRes0[[i]] = estRes
}
# save(estRes0, file = "../data/estRes0.rda")
# 
# 
# for (i in 1:length(Ns))
# {
#   cat(Ns[i], " & ", 2*Ns[i], " & ")
#   rr = estRes[[i]]
#   theta_est = rr[[3]]
#   est = round((rowMeans(theta_est-theta))*10,2)
#   SD = round(sqrt(rowMeans((theta_est-rowMeans(theta_est))^2))*10,2)
#   cat(paste0(paste(est, SD, sep = " (", collapse = ") & "), ")"))
#   kbias = round(rowMeans(rr[[1]][3:4,])*10, 2)
#   ksd = round(sqrt(rowMeans(rr[[2]][3:4,])), 2)
#   cat(" & ", paste0(paste(kbias, ksd, sep = " (", collapse = ") & "), ")"), "\\\\ \n")
# }
# 
