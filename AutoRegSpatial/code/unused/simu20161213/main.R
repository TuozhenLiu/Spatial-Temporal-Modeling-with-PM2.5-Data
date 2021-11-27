
### C_e(s) = sigy2*exp(-beta1^2|s|)
### objective function: \sum_{ij}(\hat sigma_{ij} - sigma_{ij})^2; sigma_{ij} = cov(Y_i, Y_j)
### for newton raphson algorithm see the notebook



library(MASS)
library(Matrix)
source("func.R")


N = 50; Time = 200
rhos = rep(0.5, N)
simu_res = simu.Y(beta1 = 4, sigy2 = 1, N, Time, rhos = rhos)
Y = simu_res$Y; dist_loc2 = simu_res$dist_loc2
rhos = estRho(Y)
(aa = lse.theta(Y, dist_loc2, rhos = rhos))


### plot the objective function
beta1 = 4
bs = seq(2,7,0.1); obj = rep(0, length(bs))
sig_hat = tcrossprod(Y)/ncol(Y)
for (i in 1:length(bs))
{
  b1 = bs[i]
  sig = exp(-b1^2*dist_loc2)
  sig_del = sig_hat - sig
  obj[i] = sum(sig_del^2)
}

plot(bs, obj, type = "l")
bs[which.min(obj)]

S = seq(0.5, 4, 0.1); obj = rep(0, length(S))
for (i in 1:length(S))
{
  ss = S[i]
  sig = ss*exp(-beta1^2*dist_loc2)
  sig_del = sig_hat - sig
  obj[i] = sum(sig_del^2)
}

plot(S, obj, type = "l")
S[which.min(obj)]

### 
set.seed(1234)
Ns = c(50, 100, 200, 400)
Times = c(200, 500, 1000, 2000)*sqrt(Ns)/10
Nrep = 500
betaE = 4; sig2E = 1
betaX = 2; sig2X = 0.5; alpha = 1

thetaEst = rep(list(matrix(0, nrow = 2, ncol = Nrep)), 4)
thetaRho = rep(list(matrix(0, nrow = 3, ncol = Nrep)), 4)
thetaTRho = rep(list(matrix(0, nrow = 3, ncol = Nrep)), 4)
#rhosEst = rep(list(matrix(0, nrow = N, ncol = Nrep)), 4)

for (k in 1:4)
{
  Time = Times[k]
  N = Ns[k]
  for (i in 1:Nrep)
  {
    cat(i, "\r")
    dist_loc2 = simu.dist(N)
    rhos = simu.rho(beta1 = betaX, sig2 = sig2X, N, alpha = alpha, dist_loc2)
    Y = simu.Y(beta1 = betaE, sigy2 = sig2E, N = N, Time = Time, rhos = rhos, dist_loc2)
    #rhosEst[[k]][,i] = estRho(Y)
    rhos_est = estRho(Y)
    (thetaRho[[k]][,i] = estThetaRho(rhos_est, dist_loc2))
    (thetaTRho[[k]][,i] = estThetaRho(rhos, dist_loc2))
    thetaEst[[k]][,i] = lse.theta(Y, dist_loc2,rhos = rhos_est)
    #etaEst[k,i] = (lse.theta(Y, dist_loc2))
  }
  cat("N: ", N, " Time: ", Time, #"\n Rho Bias: ", mean(rhosEst[[k]] - rhos), 
      #"\n Rho SD: ", mean((rhosEst[[k]] - rowMeans(rhosEst[[k]]))^2),
      "\n ThetaRho Bias: ", rowMeans(thetaRho[[k]]) - c(alpha, betaX, sig2X), 
      "\n ThetaRho SD: ", rowMeans((thetaRho[[k]] - rowMeans(thetaRho[[k]]))^2),
      "\n ThetaTRho Bias: ", rowMeans(thetaTRho[[k]]) - c(alpha, betaX, sig2X), 
      "\n ThetaTRho SD: ", rowMeans((thetaTRho[[k]] - rowMeans(thetaTRho[[k]]))^2),
      "\n Theta Bias: ", rowMeans(thetaEst[[k]]) - c(betaE, sig2E), 
      "\n Theta SD", rowMeans((thetaEst[[k]] - rowMeans(thetaEst[[k]]))^2) , "\n")
}


### Kriging



