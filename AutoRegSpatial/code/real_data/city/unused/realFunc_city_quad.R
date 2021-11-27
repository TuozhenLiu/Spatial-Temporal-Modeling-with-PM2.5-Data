
### Simulation： for bootstrap used in estimation


simu.loc<-function(N)
{
  X_loc = runif(N, 0, sqrt(N))
  Y_loc = runif(N, 0, sqrt(N))
  loc = cbind(X_loc, Y_loc)
  return(loc)
}


# rational quadratic covariance model

cov.quad <- function(x, beta, sig2)
{
  cc = sig2/(1+x^2/beta^2)
  return(cc)
}


# for simulating autocoefficients rho
simu.rho <- function(beta, sig2, N, alpha, dist_loc2)
{
  cov_X = cov.quad(dist_loc2, beta, sig2)
  
  #### simulate X
  eig_X = eigen(cov_X)
  eig_X$values[eig_X$values<0] = 0
  sqrt_value = sqrt(eig_X$values)
  X = eig_X$vectors%*%(sqrt_value*rnorm(N, 0, 1)) # generate X following cov_X
  #rhos = 2*pnorm(X+alpha)-1
  rhos = 2*exp(X+alpha)/(1+exp(X+alpha))-1 # generate rhos across logistic function
  return(rhos)
}



# for simulating autocoefficients rho N*R matrix
simu.rho.rep <- function(beta, sig2, N, R, alpha, dist_loc2)
{
  cov_X = cov.quad(dist_loc2, beta, sig2)
  
  #### simulate X
  eig_X = eigen(cov_X)
  eig_X$values[eig_X$values<0] = 0
  sqrt_value = sqrt(eig_X$values)
  X = eig_X$vectors%*%(sqrt_value*matrix(rnorm(N*R, 0, 1), nrow = N)) # generate X following cov_X
  #rhos = 2*pnorm(X+alpha)-1
  rhos = 2*exp(X+alpha)/(1+exp(X+alpha))-1 # generate rhos across logistic function
  return(as.vector(rhos))
}

simu.mu <- function(beta, sig2, N, mu, dist_loc2)
{
  cov_X = cov.quad(dist_loc2, beta, sig2)
  
  #### simulate X
  eig_X = eigen(cov_X)
  eig_X$values[eig_X$values<0] = 0
  sqrt_value = sqrt(eig_X$values)
  mus = eig_X$vectors%*%(sqrt_value*rnorm(N, 0, 1)) + mu # generate X following cov_X
  return(mus)
}

simu.mu.rep <- function(beta, sig2, N, R, mu, dist_loc2)
{
  cov_X = cov.quad(dist_loc2, beta, sig2)
  
  #### simulate X
  eig_X = eigen(cov_X)
  eig_X$values[eig_X$values<0] = 0
  sqrt_value = sqrt(eig_X$values)
  X = eig_X$vectors%*%(sqrt_value*matrix(rnorm(N*R, 0, 1), nrow = N)) # generate X following cov_X
  #rhos = 2*pnorm(X+alpha)-1
  mus = X + mu # generate rhos across logistic function
  return(as.vector(mus))
}

# simulate Y_t matrix following autoregression model with spatial dependence
simu.Y <- function(beta, sigy2, N, Time, rhos, dist_loc2, mu)
{
  #### simulate covariance of epsilon
  cov_e = cov.quad(dist_loc2, beta, sigy2)
  
  #### simulate epsilon
  eig_e = eigen(cov_e)
  eig_e$values[eig_e$values<0] = 0
  sqrt_value = sqrt(eig_e$values)
  # eps_mat = sapply(1:Time, function(t){ # epsilon follows same distribution across time
  #   eps0 = rnorm(N, 0, 1)
  #   eps = eig_e$vectors%*%(sqrt_value*eps0)
  # })
  eps_mat = eig_e$vectors%*%(sqrt_value*matrix(rnorm(N*Time, 0, 1), nrow = N))
  Y = matrix(0, nrow = N, ncol = Time)
  Y[,1] = eps_mat[,1]
  for (t in 2:Time)
  {
    Y[,t] = Y[,t-1]*rhos + eps_mat[,t] + mu # generate Yt according to AMSD model
  }
  return(Y = Y)
}

simu.Y.rep <- function(beta, sigy2, N, R, Time, rhos, dist_loc2, mus)
{
  #### simulate covariance of epsilon
  cov_e = cov.quad(dist_loc2, beta, sigy2)
  
  #### simulate epsilon
  eig_e = eigen(cov_e)
  eig_e$values[eig_e$values<0] = 0
  sqrt_value = sqrt(eig_e$values)
  eig_e$values[eig_e$values<0] = 0
  eps_mat_list = lapply(1:R, function(r) eig_e$vectors%*%(sqrt_value*matrix(rnorm(N*Time, 0, 1), nrow = N)))
  eps_mat_R = do.call(rbind, eps_mat_list)
  
  #eps_mat = eig_e$vectors%*%(sqrt_value*matrix(rnorm(N*Time, 0, 1), nrow = N))
  
  # eps_mat = sapply(1:Time, function(t){ # epsilon follows same distribution across time
  #   eps0 = rnorm(N, 0, 1)
  #   eps = eig_e$vectors%*%(sqrt_value*eps0)
  # })
  Y = matrix(0, nrow = N*R, ncol = Time)
  Y[,1] = eps_mat_R[,1]
  for (t in 2:Time)
  {
    Y[,t] = Y[,t-1]*rhos + eps_mat_R[,t] + mus # generate Yt according to AMSD model
  }
  return(Y = Y)
}
### Estimation functions

## estimate rho parameters
# estimate theta_X
estThetaRho<-function(rhos, dist_loc2)
{
  #X1 = qnorm((rhos+1)/2)
  rhos[which(abs(rhos)>1)] = sign(rhos[which(abs(rhos)>1)])*rep(0.999, sum(abs(rhos)>1))
  y1 = (rhos+1)/2
  X1 = log(y1/(1-y1))
  alpha = mean(X1)
  thetaX = lse.X(X1 - alpha, dist_loc2)
  return(c(alpha, thetaX))
}

estThetaMu<-function(mus, dist_loc2)
{
  #X1 = qnorm((rhos+1)/2)
  mu = mean(mus)
  thetaMu = lse.X(mus-mu, dist_loc2)
  return(c(mu, thetaMu))
}

# estimate the parameters in rhos (theta_X)
lse.X<-function(Y, dist_loc2) # Y here is X-alpha, which is centered
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
  return(abs(theta)) # since exponential covariance model is symmetric of beta, we restrict theta to be positive
}

# estimate rho and mu by iterations
### Estimation functions

estRhoMu<-function(Y)
{
  Time = ncol(Y)
  rhosmu = apply(Y, 1, function(y){
    x = cbind(1, y[-Time])
    return(solve(crossprod(x))%*%crossprod(x, y[-1]))  # estimate rho for each location
  })
  #mu = mean(rhosmu[1,])
  return(rhosmu = rhosmu)
}

## estimate epsilon parameters
# filter to obtain residuals
filter<-function(Ymat, rhos,mu)
{
  Time = ncol(Ymat)
  #Ymat1 = cbind(mu, rhos, Ymat)
  #eps = apply(Ymat1, 1, function(x) x[-(1:3)]-x[2]*x[-c(1,2,ncol(Ymat1))] - x[1])
  eps = Ymat[,-1] - rhos*Ymat[,-Time] - mu
  return(eps)
}

# estimate the parameters epsilon theta_e
lse.theta<-function(Y, dist_loc2, rhos, mu)
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
    cat(mean(abs(del)), " ", theta, "\n")
    del = lse.step(sig_hat, beta = theta[1], sigy2 = theta[2], dist_loc2) # each step
    theta = theta - del*0.5
    if (any(theta<0|theta>5)) # if the value is beyond this, it might not converge due to the initial values
      theta = runif(2, 0.1,5)
    iter = iter+1
  }
  return(abs(theta))
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
  hmat[1,1] = -2*sum(sig_del*sec_beta) + 2*sum(first_beta^2)
  hmat[1,2] = -2*sum(sig_del*sec_betasig2)+ 2*sum(first_beta*first_sig2)
  hmat[2,1] = hmat[1,2]
  hmat[2,2] = 2*sum(first_sig2^2)
  if (any(!is.finite(hmat)))
    return(grad_para)
  (eig_hmat = eigen(hmat))
  if (any(eig_hmat$values<0))
    hmat = (eig_hmat$vectors)%*%(abs(eig_hmat$values)*t(eig_hmat$vectors))
  
  del = solve(hmat)%*%grad_para #grad_beta/hmat
  
  return(del)
}





# to obtain the kriging surface
trans2mat<-function(Yt_pred, lattice_num)
{
  Yt_pred_mat = t(matrix(Yt_pred, nrow = lattice_num))
  colnames(Yt_pred_mat) = round(seq(116, 117.1, length.out = lattice_num), 3)
  rownames(Yt_pred_mat) = round(seq(39.52, 40.53, length.out = lattice_num), 3)
  return(Yt_pred_mat)
}

