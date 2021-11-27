
### C_e(s) = sigy2*exp(-beta1^2|s|)
### objective function: \sum_{ij}(\hat sigma_{ij} - sigma_{ij})^2; sigma_{ij} = cov(Y_i, Y_j)
### for newton raphson algorithm see the notebook



library(MASS)
library(Matrix)

cov.exp <- function(x, beta, sig2)
{
  return(sig2*exp(-beta^2*x))
}

simu.Y <- function(beta1, sigy2, N, Time, rhos)
{
  X_loc = runif(N, 0, sqrt(400)/10)
  Y_loc = runif(N, 0, sqrt(400)/10)
  loc = cbind(X_loc, Y_loc)
  dist_loc2 = as.matrix(dist(loc))
  
  #### simulate epsilon
  cov_e = cov.exp(dist_loc2, beta1, sigy2)
  
  #### simulate epsilon
  eig_e = eigen(cov_e)
  sqrt_value = sqrt(eig_e$values)
  eps_mat = sapply(1:Time, function(t){
    eps0 = rnorm(N, 0, 1)
    eps = eig_e$vectors%*%(sqrt_value*eps0)
  })
  Y = matrix(0, nrow = N, ncol = Time)
  Y[,1] = eps_mat[,1]
  for (t in 2:Time)
  {
    Y[,t] = Y[,t-1]*rhos + eps_mat[,t]
  }
  
  return(list(Y = Y, dist_loc2 = dist_loc2))
}

estRho<-function(Y)
{
  Time = ncol(Y)
  rhos = apply(Y, 1, function(y){
    sum(y[-Time]*y[-Time])^(-1)*sum(y[-Time]*y[-1])
  })
  return(rhos)
}



lse.step <- function(sig_hat, beta1, sigy2, dist_loc2)
{
  exp_loc = exp(-beta1^2*dist_loc2)
  sig = sigy2*exp_loc
  sig_del = sig_hat - sig
  sig_del2 = sig_hat - 2*sig
  
  grad_beta1 = 4*beta1*sum(sig_del*sig*dist_loc2)
  grad_sigy2 = -2*sum(sig_del*exp_loc)
  grad_para = c(grad_beta1, grad_sigy2)
  
  ### hessian matrix
  hmat = matrix(0,2,2)
  hmat[1,1] = (-8*sum(sig_del2*sig*beta1^2*dist_loc2^2) + 4*sum(sig_del*sig*dist_loc2))
  hmat[1,2] = 4*beta1*sum(sig_del2*exp_loc*dist_loc2)
  hmat[2,1] = hmat[1,2]
  hmat[2,2] = 2*sum(exp_loc^2)
  eig_hmat = eigen(hmat)
  if (any(eig_hmat$values<0))
    hmat = (eig_hmat$vectors)%*%(abs(eig_hmat$values)*t(eig_hmat$vectors))
  
  del = solve(hmat)%*%grad_para #grad_beta1/hmat
  return(del)
}


lse.theta<-function(Y, dist_loc2, rhos)
{
  Time = ncol(Y)
  Y0 = Y
  Y = Y[,-1] - rhos*Y[,-Time]
  if (is.null(dim(Y)))
    sig_hat = tcrossprod(Y)
  else
    sig_hat = tcrossprod(Y)/ncol(Y)
  
  #yy = log(sig_hat^2)[upper.tri(dist_loc2)]
  #xx = cbind(-dist_loc2[upper.tri(dist_loc2)], 1)
  #theta = as.vector(solve(crossprod(xx))%*%colSums(xx*yy))
  #theta[1] = sqrt(theta[1]);  theta[2] = exp(theta[2])
  #theta = sqrt(abs(sum(xx*yy)/(4*sum(xx^2))))
  theta = c(2, 2)
  iter = 1
  del = 1
  while(mean(abs(del))>10^-3&iter<1000)
  {
    #cat(mean(abs(del)), " ", theta, "\n")
    del = lse.step(sig_hat, beta1 = theta[1], sigy2 = theta[2], dist_loc2)
    theta = theta - del*0.5
    if (any(theta<0|theta>10))
      theta = runif(2, 2,5)
    iter = iter+1
  }
  return(abs(theta))
}



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
Ns = c(200, 500, 1000, 2000)
Times = c(200, 500, 1000, 2000)
Nrep = 1000
beta1 = 4; sigy2 = 1
rhos = runif(N, 0.2, 0.5)

thetaEst = rep(list(matrix(0, nrow = 2, ncol = Nrep)), 4)
rhosEst = rep(list(matrix(0, nrow = N, ncol = Nrep)), 4)
for (k in 1:4)
{
  Time = Times[k]
  N = 50
  for (i in 1:Nrep)
  {
    cat(i, "\r")
    simu_res = simu.Y(beta1 = 4, sigy2 = 1, N = N, Time = Time, rhos = rhos)
    Y = simu_res$Y; dist_loc2 = simu_res$dist_loc2
    rhosEst[[k]][,i] = estRho(Y)
    thetaEst[[k]][,i] = (lse.theta(Y, dist_loc2,rhos = rhosEst[[k]][,i]))
    #etaEst[k,i] = (lse.theta(Y, dist_loc2))
  }
  cat("Time: ", Time, "\n Rho Bias: ", mean(rhosEst[[k]] - rhos), 
      "\n Rho SD: ", mean((rhosEst[[k]] - rowMeans(rhosEst[[k]]))^2),
      "\n Theta Bias: ", rowMeans(thetaEst[[k]]) - c(beta1, sigy2), 
      "\n Theta SD", rowMeans((thetaEst[[k]] - rowMeans(thetaEst[[k]]))^2) , "\n")
}


### Kriging



