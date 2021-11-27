
tr<-function(M)
{
  return(sum(diag(M)))
}



infer.mu<-function(thetaMu, mu_true, dist_loc2, Time)
{
  mu = thetaMu[1]
  N = nrow(dist_loc2); Time = ncol(Y)
  #mu_cov = sum(tcrossprod(eps))/(N^2*Time^2)
  
  sigE = cov.quad(dist_loc2, thetaMu[2], sig2 = thetaMu[3])
  mu_cov = sum(sigE)/N^2
  CI = (mu-1.96*sqrt(mu_cov)<mu_true)&(mu_true<mu+1.96*sqrt(mu_cov))
  return(list(mu_cov, CI))
}



infer.alpha<-function(thetaX, alpha_true, dist_loc2, Time)
{
  N = nrow(dist_loc2); Time = ncol(Y)
  alpha = thetaX[1]
  sigE = cov.quad(dist_loc2, thetaX[2], sig2 = thetaX[3])
  alpha_cov = sum(sigE)/N^2
  CI = (alpha-1.96*sqrt(alpha_cov)<alpha_true)&(alpha_true<alpha+1.96*sqrt(alpha_cov))
  return(list(alpha_cov, CI))
}


## estimate rho parameters
# estimate theta_X
estThetaRho<-function(rhos, dist_loc2, theta_true)
{
  #X1 = qnorm((rhos+1)/2)
  rhos[which(abs(rhos)>1)] = sign(rhos[which(abs(rhos)>1)])*rep(0.999, sum(abs(rhos)>1))
  y1 = (rhos+1)/2
  X1 = log(y1/(1-y1))
  alpha = mean(X1)
  thetaX = lse.X(X1 - alpha, dist_loc2, theta_true)
  thetaX[[1]] = c(alpha, thetaX[[1]])
  return(thetaX)
}


estThetaMu<-function(mus, dist_loc2, theta_true)
{
  mu = mean(mus)
  thetaX = lse.X(mus-mu, dist_loc2, theta_true)
  thetaX[[1]] = c(mu, thetaX[[1]])
  return(thetaX)
}



# estimate the parameters in rhos (theta_X)
lse.X<-function(Y, dist_loc2, theta_true) # Y here is X-alpha, which is centered
{
  sig_hat = tcrossprod(Y) 
  
  theta = c(1, var(as.vector(Y)))
  iter = 1; del = 1
  while(mean(abs(del))>10^-3&iter<1000)
  {
    #cat(mean(abs(del)), " ", theta, "\n")
    del = lse.step(sig_hat, beta = theta[1], sigy2 = theta[2], dist_loc2)
    theta = theta - del*0.5
    if (any(theta<0|theta>10))
      theta = runif(2, 0.1,10)
    iter = iter+1
  }
  theta = abs(theta)
  est_cov = infer.cov(sig_hat, beta = theta[1], sigy2 = theta[2], dist_loc2, Time = 1)
  CI = (theta-1.96*sqrt(diag(est_cov))<theta_true)&(theta_true<theta+1.96*sqrt(diag(est_cov)))
  
  return(list(abs(theta), est_cov, CI)) # since exponential covariance model is symmetric of beta, we restrict theta to be positive
}


lse.theta<-function(Y, dist_loc2, rhos, mu, theta_true)
{
  Time = ncol(Y)
  Y0 = Y
  Y = Y[,-1] - rhos*Y[,-Time] - mu # residuals
  if (is.null(dim(Y)))
    sig_hat = tcrossprod(Y)
  else
    sig_hat = tcrossprod(Y)/ncol(Y)
  
  theta = c(1, 0.5) # the initial value
  iter = 1; del = 1
  while(mean(abs(del))>10^-3&iter<1000)
  {
    #cat(mean(abs(del)), " ", theta, "\n")
    del = lse.step(sig_hat, beta = theta[1], sigy2 = theta[2], dist_loc2) # each step
    theta = theta - del*0.5
    if (any(theta<0|theta>5)) # if the value is beyond this, it might not converge due to the initial values
      theta = runif(2, 0.1,5)
    iter = iter+1
  }
  theta = abs(theta)
  est_cov = infer.cov(sig_hat, beta = theta[1], sigy2 = theta[2], dist_loc2, Time)
  CI = (theta-1.96*sqrt(diag(est_cov))<theta_true)&(theta_true<theta+1.96*sqrt(diag(est_cov)))
  return(list(abs(theta), est_cov, CI))
}


# for each newton-raphson iteration, the following function gives each step
# only works for rational quadratic model
lse.step <- function(sig_hat, beta, sigy2, dist_loc2)
{
  beta_s = beta^2+dist_loc2^2
  sig = sigy2*beta^2/beta_s
  
  sig_del = sig_hat - sig
  
  ### first order derivative: gradient
  first_beta = 2*beta*sigy2*dist_loc2^2/beta_s^2
  first_sig2 = beta^2/beta_s
  grad_beta = - 2*sum(sig_del*first_beta)
  grad_sigy2 = -2*sum(sig_del*first_sig2)
  
  
  grad_para = c(grad_beta, grad_sigy2)
  
  ### hessian matrix
  sec_beta = (-8*sigy2*beta^2*dist_loc2^2 + 2*sigy2*dist_loc2^2*beta_s)/beta_s^3
  sec_betasig2 = 2*beta*dist_loc2^2/beta_s^2
  
  
  hmat = matrix(0,2,2)
  hmat[1,1] = -sum(sig_del*sec_beta) + sum(first_beta^2)
  hmat[1,2] = -sum(sig_del*sec_betasig2)+ sum(first_beta*first_sig2)
  hmat[2,1] = hmat[1,2]
  hmat[2,2] = sum(first_sig2^2)
  if (any(!is.finite(hmat)))
    return(grad_para)
  (eig_hmat = eigen(hmat))
  if (any(eig_hmat$values<0))
    hmat = (eig_hmat$vectors)%*%(abs(eig_hmat$values)*t(eig_hmat$vectors))
  
  del = solve(hmat)%*%grad_para #grad_beta/hmat
  
  return(del)
}


infer.cov<-function(sig_hat, beta, sigy2, dist_loc2, Time)
{
  beta_s = beta^2+dist_loc2^2
  sig = sigy2*beta^2/beta_s
  
  sig_del = sig_hat - sig
  
  ### first order derivative: gradient
  first_beta = 2*beta*sigy2*dist_loc2^2/beta_s^2
  first_sig2 = beta^2/beta_s
  grad_beta = - 2*sum(sig_del*first_beta)
  grad_sigy2 = -2*sum(sig_del*first_sig2)
  
  
  grad_para = c(grad_beta, grad_sigy2)
  
  ### hessian matrix
  sec_beta = (-8*sigy2*beta^2*dist_loc2^2 + 2*sigy2*dist_loc2^2*beta_s)/beta_s^3
  sec_betasig2 = 2*beta*dist_loc2^2/beta_s^2
  
  
  hmat = matrix(0,2,2)
  hmat[1,1] = -2*sum(sig_del*sec_beta) + 2*sum(first_beta^2)
  hmat[1,2] = -2*sum(sig_del*sec_betasig2)+ 2*sum(first_beta*first_sig2)
  hmat[2,1] = hmat[1,2]
  hmat[2,2] = 2*sum(first_sig2^2)
  
  ### first order covariance
  W1 = first_beta
  W2 = first_sig2
  var_grad = matrix(0, 2,2)
  var_grad[1,1] = 2*tr(W1%*%sig%*%W1%*%sig)
  var_grad[1,2] = 2*tr(W1%*%sig%*%W2%*%sig)
  var_grad[2,1] = var_grad[1,2]
  var_grad[2,2] = 2*tr(W2%*%sig%*%W2%*%sig)
  var_grad = var_grad*4/Time
  
  return(solve(hmat)%*%var_grad%*%solve(hmat))
}

