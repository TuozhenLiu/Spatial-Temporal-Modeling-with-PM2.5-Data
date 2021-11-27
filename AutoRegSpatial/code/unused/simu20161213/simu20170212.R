setwd("E:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/")

source("../real_data/krigingMuFunc.R")
source("../real_data/city/realFunc_city.R")
source("../real_data/predictFunc.R")

Ns = c(100, 200, 400)
Times = c(200, 400, 800)

mu = 1; betaE = 2; sig2E = 3
betaX = 1.5; sig2X = 1; alpha = 2

#betaE = 1; sig2E = 2
#betaX = 2; sig2X = 2; alpha = 0
thetaX_true = c(alpha, betaX, sig2X)
theta = c(alpha, betaX, sig2X, betaE, sig2E, mu)

Nrep = 200
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
    thetaX_true_est = matrix(0, nrow = 3, ncol = Nrep)
    for (r in 1:Nrep)
    {
      loc = simu.loc(N)
      dist_loc2 = as.matrix(dist(loc))
      rhos = simu.rho(beta1 = betaX, sig2 = sig2X, N, alpha = alpha, dist_loc2)
      Y = simu.Y(beta1 = betaE, sigy2 = sig2E, N = N, Time = Time, rhos = rhos, dist_loc2, mu = mu)
      
      
      train_ind = 1:floor(2/3*N)
      krig_ind = setdiff(1:N, train_ind)
      dist_train = dist_loc2[train_ind, train_ind]
      dist_krig = dist_loc2[train_ind, krig_ind]
      Y_train = Y[train_ind,]
      rhos_train = rhos[train_ind]
      
      rhos_mu_est = estRhoMuIter(Y_train)
      mu_est = rhos_mu_est[1,1]
      rhos_est = rhos_mu_est[2,]
      thetaE = lse.theta(Y_train, dist_train, rhos = rhos_est, mu = mu_est)
      thetaX = estThetaRho(rhos_est, dist_train)
      #thetaX_true_est[,r] = estThetaRho(rhos_train, dist_train)
      theta_est[,r] = c(thetaX, thetaE, mu_est)
      
      cat(r, "\r")
    }
    estRes[[t]] = list(res0, res, theta_est, thetaX_true_est)
    cat("\n N: ", N, " Time:", Time, "\n",
        "Estimation Bias: ", rowMeans(theta_est-theta), "\n",
        #"thetaXTrueRho Bias: ", rowMeans(thetaX_true_est-thetaX_true), "\n",
        "Estimation SD: ", sqrt(rowMeans((theta_est-rowMeans(theta_est))^2)), "\n",
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

