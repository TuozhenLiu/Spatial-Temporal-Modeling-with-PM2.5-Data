
### 1. Basic functions
# generate the locations
simu.loc<-function(N, cov.type = "exp")
{
  X_loc = runif(N, 0, sqrt(N))
  Y_loc = runif(N, 0, sqrt(N))
  loc = cbind(X_loc, Y_loc)
  if (cov.type=="quad")
    loc = loc*sqrt(10)
  return(loc)
}

# get trace of a matrix
tr<-function(M)
{
  return(sum(diag(M)))
}

# rho transformation
rho.func<-function(x, alpha)
{
  return(2*exp(x+alpha)/(1+exp(x+alpha))-1)
}


# rational quadratic cov
cov.quad <- function(x, beta, sig2)
{
  cc = sig2/(1+x^2/beta^2)
  return(cc)
}

# exponential covariance model
cov.exp <- function(x, beta, sig2)
{
  return(sig2*exp(-beta^2*x))
}


### 2. Parameter estimation and inference


# # estimate rho and mu by iterations (for mu is a constant, not used in this paper)
# estRhoMuIter<-function(Y)
# {
#   Time = ncol(Y); N = nrow(Y)
#   rhos = rep(0.5, nrow(Y))
#   mu = 0
#   del = 1
#   while (del>10^{-4})
#   {
#     #show(del)
#     mu1 = sum(Y[,-1] - rhos*Y[,-Time])/(N*Time)
#     rhos1 = rowSums(Y[,-Time]*(Y[,-1]-mu1))/rowSums(Y[,-Time]^2)
#     del = sum(abs(c(mu1-mu, rhos1-rhos)))/(N+1)
#     mu = mu1; rhos = rhos1
#   }
#   
#   return(list(mu, rhos))
# }

# Step I: estimate rho mu for each location

estRhoMu<-function(Y)
{
  Time = ncol(Y)
  rhosmu = apply(Y, 1, function(y){
    x = cbind(1, y[-Time])
    return(solve(crossprod(x))%*%crossprod(x, y[-1]))  # estimate rho for each location
  })
  
  return(rhosmu = rhosmu)
}

# inference for mu (mean)
infer.mu<-function(thetaMu, mu_true, dist_loc2, Time, cov.type = "exp")
{
  mu = thetaMu[1]
  N = nrow(dist_loc2)
  #mu_cov = sum(tcrossprod(eps))/(N^2*Time^2)
  if (cov.type=="exp")
    sigE = cov.exp(dist_loc2, thetaMu[2], sig2 = thetaMu[3])
  else
    sigE = cov.quad(dist_loc2, thetaMu[2], sig2 = thetaMu[3])
  mu_cov = sum(sigE)/N^2
  CI = (mu-1.96*sqrt(mu_cov)<mu_true)&(mu_true<mu+1.96*sqrt(mu_cov))
  return(list(sqrt(mu_cov), CI))
}


# Step II: estimate rho and mu parameters
# estimate theta_X
estThetaRho<-function(rhos, dist_loc2, infer = T, theta_true = c(1,1,1), cov.type = "exp", p.value = F, est_sig = T)
{
  rhos[which(abs(rhos)>1)] = sign(rhos[which(abs(rhos)>1)])*rep(0.9999, sum(abs(rhos)>1))
  y1 = (rhos+1)/2
  X1 = log(y1/(1-y1))
  thetaX = estThetaMu(X1, dist_loc2, infer, theta_true, cov.type, p.value, est_sig = est_sig)
  return(thetaX)
}


# estimate theta_Mu
estThetaMu<-function(mus, dist_loc2, infer = T, theta_true = c(1,1,1), cov.type = "exp", p.value = F, est_sig = T)
{
  mu_hat = mean(mus)
  thetaX = lse.X(mus-mu_hat, dist_loc2, infer, theta_true[2:3], cov.type, p.value, est_sig = est_sig)
  if (!infer)
    return(c(mu_hat, thetaX))
  
  theta = thetaX[[1]][1]
  if (!infer)
    return(c(mu_hat, thetaX[[1]]))
  
  mu_list = infer.mu(thetaMu = c(mu_hat, theta, theta_true[3]), 
                     mu_true = theta_true[1], dist_loc2, Time = 1, cov.type)
  
  thetaX[[1]] = c(mu_hat, thetaX[[1]])
  thetaX[[2]] = c(mu_list[[1]], thetaX[[2]])
  thetaX[[3]] = c(mu_list[[2]], thetaX[[3]])
  if (p.value)
    thetaX[[4]] = (1-pnorm(abs(thetaX[[1]]/thetaX[[2]])))*2
  return(thetaX)
}


# estimate theta_e (need to filter first)
lse.theta<-function(Y, dist_loc2, rhos, mu, infer = T, theta_true = c(1,1), cov.type = "exp", p.value = F, est_sig = T)
{
  Time = ncol(Y)
  Y0 = Y
  Y = Y[,-1] - rhos*Y[,-Time] - mu # filter to obtain residuals
  
  thetaE_list = lse.X(Y, dist_loc2, infer, theta_true, cov.type, p.value, est_sig = est_sig)
  if (p.value)
    thetaE_list[[4]] = (1-pnorm(abs(thetaE_list[[1]]/thetaE_list[[2]])))*2
  
  return(thetaE_list)
  
}

### 3. general functions need to be called
# estimate decentered processes (no replicates in time)
lse.X<-function(Y, dist_loc2, infer = T, theta_true, cov.type = "exp", p.value = F, est_sig = T) # Y here is X-alpha, which is centered
{
  Time = ncol(Y)
  if(is.null(Time))
    Time = 1
  if (is.null(dim(Y)))
    sig_hat = tcrossprod(Y)
  else
    sig_hat = tcrossprod(Y)/ncol(Y)
  
  sig2_hat = mean(Y^2)#var(as.vector(Y))
  
  theta = c(1.5)
  iter = 1; del = 1
  if (cov.type == "quad")
  {
    lse.step = lse.step.quad
    infer.cov = infer.cov.quad
  }
  if ((infer&p.value)|est_sig)
    theta_true[2] = sig2_hat
  
  while(mean(abs(del))>10^-3&iter<1000)
  {
    del = lse.step(sig_hat, beta = theta, sigy2 = theta_true[2], dist_loc2)
    theta = theta - del*0.5
    if (any(theta<0|theta>5))
      theta = runif(1, 0.1,5)
    iter = iter+1
  }
  theta = abs(theta)
  if (!infer)
    return(c(theta, sig2_hat))
  
  # for spatial parameter inference
  est_cov = infer.cov(sig_hat, beta = theta, sigy2 = theta_true[2], dist_loc2, Time = Time)
  CI = (theta-1.96*sqrt((est_cov))<theta_true[1])&(theta_true[1]<theta+1.96*sqrt((est_cov)))
  
  # for sigma inference
  
  if (cov.type == "quad")
    sig = cov.quad(dist_loc2, theta, sig2_hat)
  else
    sig = cov.exp(dist_loc2, theta, sig2_hat)
  
  sig2_cov = 2*tr(sig%*%sig)/nrow(dist_loc2)^2/Time
  sig2_CI = (sig2_hat-1.96*sqrt((sig2_cov))<theta_true[2])&(theta_true[2]<sig2_hat+1.96*sqrt((sig2_cov)))
  
  return(list(c(abs(theta), sig2_hat), c(sqrt(est_cov), sqrt(sig2_cov)), c(CI, sig2_CI)))
  
  #return(list(abs(theta), est_cov, CI)) # since exponential covariance model is symmetric of beta, we restrict theta to be positive
}


# newton raphson for exp function

lse.step <- function(sig_hat, beta, sigy2, dist_loc2)
{
  exp_loc = exp(-beta^2*dist_loc2)
  sig = sigy2*exp_loc
  sig_del = sig_hat - sig
  sig_del2 = sig_hat - 2*sig
  
  # first order derivative: gradient
  grad_beta = 4*beta*sum(sig_del*sig*dist_loc2)
  #grad_sigy2 = -2*sum(sig_del*exp_loc)
  grad_para = c(grad_beta)
  
  # hessian matrix
  hmat = (-8*sum(sig_del2*sig*beta^2*dist_loc2^2) + 4*sum(sig_del*sig*dist_loc2))
  
  del = grad_para/hmat #grad_beta/hmat
  return(del)
}

# asymptotic covariance for exp function
infer.cov<-function(sig_hat, beta, sigy2, dist_loc2, Time)
{
  exp_loc = exp(-beta^2*dist_loc2)
  sig = sigy2*exp_loc
  sig_del = sig_hat - sig
  sig_del2 = sig_hat - 2*sig
  
  W1 = -2*beta*sig*dist_loc2
  #W2 = exp_loc
  var_grad = 2*tr(W1%*%sig%*%W1%*%sig)
  var_grad = var_grad*4/Time
  
  # hessian matrix
  hmat = (-8*sum(sig_del2*sig*beta^2*dist_loc2^2) + 4*sum(sig_del*sig*dist_loc2))
  return(var_grad/hmat^2)
}

# newton raphson for exp function
lse.step.quad <- function(sig_hat, beta, sigy2, dist_loc2)
{
  dist_loc2 = dist_loc2^2
  beta_s = beta^2+dist_loc2
  sig = sigy2*beta^2/beta_s
  
  sig_del = sig_hat - sig
  
  # first order derivative: gradient
  first_beta = 2*beta*sigy2*dist_loc2/beta_s^2
  first_sig2 = beta^2/beta_s
  grad_beta = - 2*sum(sig_del*first_beta)
  grad_sigy2 = -2*sum(sig_del*first_sig2)
  
  
  grad_para = grad_beta #c(grad_beta, grad_sigy2)
  
  # hessian matrix
  sec_beta = (-8*sigy2*beta^2*dist_loc2 + 2*sigy2*dist_loc2*beta_s)/beta_s^3
  #sec_betasig2 = 2*beta*dist_loc2/beta_s^2
  
  
  hmat = -2*sum(sig_del*sec_beta) + 2*sum(first_beta^2)
  
  del = grad_para/hmat #grad_beta/hmat
  
  return(del)
}


# asymptotic covariance for quad function
infer.cov.quad <- function(sig_hat, beta, sigy2, dist_loc2, Time)
{
  dist_loc2 = dist_loc2^2
  beta_s = beta^2+dist_loc2
  sig = sigy2*beta^2/beta_s
  
  sig_del = sig_hat - sig
  
  # first order derivative: gradient
  first_beta = 2*beta*sigy2*dist_loc2/beta_s^2
  first_sig2 = beta^2/beta_s
  grad_beta = - 2*sum(sig_del*first_beta)
  grad_sigy2 = -2*sum(sig_del*first_sig2)
  
  
  grad_para = grad_beta#c(grad_beta, grad_sigy2)
  
  # hessian matrix
  sec_beta = (-8*sigy2*beta^2*dist_loc2 + 2*sigy2*dist_loc2*beta_s)/beta_s^3
  sec_betasig2 = 2*beta*dist_loc2/beta_s^2
  
  
  hmat = -2*sum(sig_del*sec_beta) + 2*sum(first_beta^2)
  
  # first order covariance
  W1 = first_beta
  var_grad = 2*tr(W1%*%sig%*%W1%*%sig)
  var_grad = var_grad*4/Time
  
  return(var_grad/hmat^2)
}



