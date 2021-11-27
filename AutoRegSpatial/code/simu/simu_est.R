

set.seed(123)
#hist(cov.quad(dist_loc2, 3, 1))
Ns = c(200, 300)
Times = c(300, 500)


muCI = matrix(T, ncol = Nrep, nrow = 3)
mu_est = matrix(0, ncol = Nrep, nrow = 3)
mu_sd = matrix(0, ncol = Nrep, nrow = 3)
rhoCI = matrix(T, ncol = Nrep, nrow = 3)
rho_est = matrix(0, ncol = Nrep, nrow = 3)
rho_sd = matrix(0, ncol = Nrep, nrow = 3)
ECI = matrix(T, ncol = Nrep, nrow = 2)
E_est = matrix(0, ncol = Nrep, nrow = 2)
E_sd = matrix(0, ncol = Nrep, nrow = 2)

muTCI = matrix(T, ncol = Nrep, nrow = 3)
muT_est = matrix(0, ncol = Nrep, nrow = 3)
muT_sd = matrix(0, ncol = Nrep, nrow = 3)
rhoTCI = matrix(T, ncol = Nrep, nrow = 3)
rhoT_est = matrix(0, ncol = Nrep, nrow = 3)
rhoT_sd = matrix(0, ncol = Nrep, nrow = 3)
ETCI = matrix(T, ncol = Nrep, nrow = 2)
ET_est = matrix(0, ncol = Nrep, nrow = 2)
ET_sd = matrix(0, ncol = Nrep, nrow = 2)

estRes0 =list()

for (i in 1:2)
{
  N = Ns[i]
  estRes = list()
  for (t in 1:2)
  {
    Time = Times[t]
    for (r in 1:Nrep)
    {
      cat(r, " | ")
      loc = simu.loc(N,  cov.type)*2
      dist_loc2 = as.matrix(dist(loc))
      rhos = simu.rho(beta = betaX, sig2 = sig2X, N, alpha = alpha, dist_loc2, cov.type = cov.type)
      mus = simu.mu(beta = betaMu, sig2 = sig2Mu, N, mu = mu, dist_loc2, cov.type = cov.type)
      Y = simu.Y(beta = betaE, sigy2 = sig2E, N = N, Time = Time, rhos = rhos, 
                 dist_loc2, mu = mus, cov.type = cov.type)
      rhosMu = estRhoMu(Y = Y)
      mus_est = rhosMu[1,]
      rhos_est = rhosMu[2,]
      
      estMu_list = estThetaMu(mus_est, dist_loc2, infer = T, 
                              theta_true = c(mu, betaMu, sig2Mu), cov.type, est_sig = F)
      estrho_list = estThetaRho(rhos_est, dist_loc2, infer = T, 
                                theta_true = c(alpha, betaX, sig2X), cov.type, est_sig = F)
      
      estE_list = lse.theta(Y, dist_loc2, rhos = as.vector(rhos_est), mu = as.vector(mus_est), 
                            theta_true = c(betaE, sig2E), cov.type = cov.type, est_sig = F)
      
      mu_est[,r] = estMu_list[[1]]
      mu_sd[,r] = (estMu_list[[2]])
      muCI[,r] = estMu_list[[3]]
      
      rho_est[,r] = estrho_list[[1]]
      rho_sd[,r] = (estrho_list[[2]])
      rhoCI[,r] = estrho_list[[3]]
      
      E_est[,r] = estE_list[[1]]
      E_sd[,r] = (estE_list[[2]])
      ECI[,r] = estE_list[[3]]
      
      ###
      
      estMu_list = estThetaMu(mus, dist_loc2, infer = T, 
                              theta_true = c(mu, betaMu, sig2Mu), cov.type = cov.type, est_sig = F)
      estrho_list = estThetaRho(rhos, dist_loc2, infer = T, 
                                theta_true = c(alpha, betaX, sig2X), cov.type = cov.type, est_sig = F)
      
      estE_list = lse.theta(Y, dist_loc2, rhos = as.vector(rhos), mu = as.vector(mus), 
                            theta_true = c(betaE, sig2E), cov.type = cov.type, est_sig = F)
      
      muT_est[,r] = estMu_list[[1]]
      muT_sd[,r] = (estMu_list[[2]])
      muTCI[,r] = estMu_list[[3]]
      
      rhoT_est[,r] = estrho_list[[1]]
      rhoT_sd[,r] = (estrho_list[[2]])
      rhoTCI[,r] = estrho_list[[3]]
      
      ET_est[,r] = estE_list[[1]]
      ET_sd[,r] = (estE_list[[2]])
      ETCI[,r] = estE_list[[3]]
      
    }
    estRes[[t]] = list(rbind(E_est, mu_est, rho_est), rbind(ECI,muCI, rhoCI))
    cat("N: ", N, " Time: ", Time, "\n",
        "Bias: ", rowMeans(estRes[[t]][[1]]) - c(betaE, sig2E, mu, betaMu, sig2Mu, alpha, betaX, sig2X), "\n",
        "True Bias: ", rowMeans(ET_est - c(betaE, sig2E)),
        rowMeans(muT_est - c(mu, betaMu, sig2Mu)), 
        rowMeans(rhoT_est - c(alpha, betaX, sig2X)), "\n",
        "SD: ", apply(estRes[[t]][[1]],1, sd), "\n",
        "True SD: ",  apply(ET_est,1, sd),
        apply(muT_est,1, sd), apply(rhoT_est,1, sd), "\n",
        #"Asym SD: ", rowMeans(mu_sd), 
        "CP: ",  rowMeans(estRes[[t]][[2]]), "\n",
        "True CP: ",  rowMeans(ETCI), rowMeans(muTCI), rowMeans(rhoTCI), "\n")
  }
  estRes0[[i]] = estRes
}
