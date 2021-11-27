
cov.exp <- function(x, beta, sig2)
{
  return(sig2*exp(-beta^2*x))
}
simu.Y <- function(beta1, sigy2, N, Time, rhos, dist_loc2, mu)
{
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
    Y[,t] = Y[,t-1]*rhos + eps_mat[,t] + mu
  }
  return(Y = Y)
}

simu.loc<-function(N)
{
  X_loc = runif(N, 0, sqrt(N)/10)
  Y_loc = runif(N, 0, sqrt(N)/10)
  loc = cbind(X_loc, Y_loc)
  return(loc)
}

simu.dist<-function(N)
{
  X_loc = runif(N, 0, sqrt(N)/10)
  Y_loc = runif(N, 0, sqrt(N)/10)
  loc = cbind(X_loc, Y_loc)
  dist_loc2 = as.matrix(dist(loc))
  return(dist_loc2)
}

simu.rho <- function(beta1, sig2, N, alpha, dist_loc2)
{
  cov_X = cov.exp(dist_loc2, beta1, sig2)
  
  #### simulate epsilon
  eig_X = eigen(cov_X)
  sqrt_value = sqrt(eig_X$values)
  X = eig_X$vectors%*%(sqrt_value*rnorm(N, 0, 1))
  #rhos = 2*pnorm(X+alpha)-1
  rhos = 2*exp(X+alpha)/(1+exp(X+alpha))-1
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
      theta = runif(2, 0.3,5)
    iter = iter+1
  }
  return(abs(theta))
}


lse.X<-function(Y, dist_loc2)
{
  sig_hat = tcrossprod(Y)
  
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

estThetaRho<-function(rhos, dist_loc2)
{
  #X1 = qnorm((rhos+1)/2)
  y1 = (rhos+1)/2
  X1 = log(y1/(1-y1))
  alpha = mean(X1)
  thetaX = lse.X(X1 - alpha, dist_loc2)
  return(c(alpha, thetaX))
}

estRho<-function(Y)
{
  Time = ncol(Y)
  rhos = apply(Y, 1, function(y){
    sum(y[-Time]*y[-Time])^(-1)*sum(y[-Time]*y[-1])
  })
  rhos[rhos> 0.9999999] = 0.9999999
  rhos[rhos< -0.9999999] = -0.9999999
  return(rhos)
}
