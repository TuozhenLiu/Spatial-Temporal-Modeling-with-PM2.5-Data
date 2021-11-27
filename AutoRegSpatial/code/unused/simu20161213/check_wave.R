
simu.loc<-function(N)
{
  X_loc = runif(N, 0, sqrt(N))
  Y_loc = runif(N, 0, sqrt(N))
  loc = cbind(X_loc, Y_loc)
  return(loc)
}



simu.mu <- function(phi, sig2, N, mu, dist_loc2)
{
  cov_X = cov.wave(dist_loc2, phi, sig2)
  
  #### simulate X
  eig_X = eigen(cov_X)
  eig_X$values[eig_X$values<0] = 0
  sqrt_value = sqrt(eig_X$values)
  mus = eig_X$vectors%*%(sqrt_value*rnorm(N, 0, 1)) + mu # generate X following cov_X
  return(mus)
}

lse.X<-function(Y, dist_loc2) # Y here is X-alpha, which is centered
{
  sig_hat = tcrossprod(Y) 
  
  theta = 1.5
  iter = 1; del = 1
  while(mean(abs(del))>10^-3&iter<1000)
  {
    cat(mean(abs(del)), " ", theta, "\n")
    del = lse.step(sig_hat, phi = theta[1], sigy2 = 1, dist_loc2)
    theta = theta - del*0.5
    if (any(theta<0|theta>10))
      theta = runif(1, 0.1,10)
    iter = iter+1
  }
  return(abs(theta)) # since exponential covariance model is symmetric of beta, we restrict theta to be positive
}

object.func<-function(phi, sig2, sig_hat)
{
  cov_X = cov.wave(dist_loc2, phi, sig2)
  return(mean((cov_X - sig_hat)^2))
}

grad<-function(phi, sigy2, sig_hat)
{
  #exp_loc = exp(-phi^2*dist_loc2)
  #sig = sigy2*exp_loc
  cc = 6.5*pi
  phis = phi/(cc*dist_loc2); diag(phis) = phi/cc
  sin_phis = sin(1/phis); diag(sin_phis) = cc/phi
  cos_phis = cos(1/phis)
  sig = sigy2*phis*sin_phis
  
  #sig = cov.wave(x = dist_loc2, phi, sigy2)
  sig_del = sig_hat - sig
  #sig_del2 = sigy2*phis/phi*sin_phis - sigy2/phi*cos_phis
  ch = 1/(cc*dist_loc2); diag(ch) = 1/cc
  
  ### first order derivative: gradient
  first_phi = sigy2*(ch*sin_phis - 1/phi*cos_phis)
  # first_sig2 = phis*sin_phis
  # grad_sig2 = - sum(sig_del*first_sig2)
  grad_phi = - sum(sig_del*first_phi)/nrow(sig_hat)
  return(grad_phi)
}


# for each newton-raphson iteration, the following function gives each step
# only works for exponential model
lse.step <- function(sig_hat, phi, sigy2, dist_loc2)
{
  #exp_loc = exp(-phi^2*dist_loc2)
  #sig = sigy2*exp_loc
  cc = 6.5*pi
  phis = phi/(cc*dist_loc2); diag(phis) = phi/cc
  sin_phis = sin(1/phis); diag(sin_phis) = cc/phi
  cos_phis = cos(1/phis)
  sig = sigy2*phis*sin_phis
  
  #sig = cov.wave(x = dist_loc2, phi, sigy2)
  sig_del = sig_hat - sig
  #sig_del2 = sigy2*phis/phi*sin_phis - sigy2/phi*cos_phis
  ch = 1/(cc*dist_loc2); diag(ch) = 1/cc
  
  ### first order derivative: gradient
  first_phi = sigy2*(ch*sin_phis - 1/phi*cos_phis)
  # first_sig2 = phis*sin_phis
  # grad_sig2 = - sum(sig_del*first_sig2)
  grad_phi = - sum(sig_del*first_phi)
  
  # grad_phi = 4*phi*sum(sig_del*sig*dist_loc2)
  # grad_sigy2 = -2*sum(sig_del*exp_loc)
  grad_para = c(grad_phi)
  
  ### hessian matrix
  sec_phi = - sigy2/phis/phi*sin_phis
  diag(sec_phi) = -sigy2/phi^2
  #sec_phisig = first_phi/sigy2
  
  hmat = -sum(sig_del*sec_phi) + sum(first_phi^2)
  
  del = grad_phi/hmat
  return(del)
}

hist(cov.wave(dist_loc2, 3, 1))
loc = simu.loc(800); dist_loc2 = as.matrix(dist(loc))
mus = simu.mu(phi = 2, sig2 = 1, N = 800, mu = 0, dist_loc2)
sd(mus)
#lse.X(Y = mus, dist_loc2)

### first order check
xx = seq(0.1, 5, 0.5)
yy = sapply(xx, grad, sigy2 = 1, sig_hat = tcrossprod(mus))
plot(xx, yy, type = "l")

#grad(phi = 2, sigy2 = 1, sig_hat = tcrossprod(mus))

### test the obj func
xx = seq(0.1, 5, 0.5)
yy = sapply(xx, object.func, sig2 = 1, sig_hat = tcrossprod(mus))
plot(xx, yy, type = "l")



