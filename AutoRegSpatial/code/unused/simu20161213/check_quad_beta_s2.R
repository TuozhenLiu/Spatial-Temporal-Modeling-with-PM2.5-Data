
### distance takes square in this case


tr<-function(M)
{
  return(sum(diag(M)))
}



simu.loc<-function(N)
{
  X_loc = runif(N, 0, sqrt(N))*sqrt(10)
  Y_loc = runif(N, 0, sqrt(N))*sqrt(10)
  loc = cbind(X_loc, Y_loc)
  return(loc)
}

cov.quad <- function(x, beta, sig2)
{
  cc = sig2/(1+x^2/beta^2)
  return(cc)
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


object.func<-function(beta, sig2, sig_hat,dist_loc2)
{
  cov_X = cov.quad(dist_loc2, beta, sig2)
  return(mean((cov_X - sig_hat)^2))
}


# for each newton-raphson iteration, the following function gives each step
# only works for rational quadratic model
lse.step <- function(sig_hat, beta, sigy2, dist_loc2)
{
  beta_s = beta^2+dist_loc2
  sig = sigy2*beta^2/beta_s
  
  sig_del = sig_hat - sig
  
  ### first order derivative: gradient
  first_beta = 2*beta*sigy2*dist_loc2/beta_s^2
  first_sig2 = beta^2/beta_s
  grad_beta = - 2*sum(sig_del*first_beta)
  grad_sigy2 = -2*sum(sig_del*first_sig2)
  
  
  grad_para = c(grad_beta, grad_sigy2)
  
  ### hessian matrix
  sec_beta = (-8*sigy2*beta^2*dist_loc2 + 2*sigy2*dist_loc2*beta_s)/beta_s^3
  sec_betasig2 = 2*beta*dist_loc2/beta_s^2
  
  
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
###########################################################################


lse.X<-function(Y, dist_loc2, theta_true) # Y here is X-alpha, which is centered
{
  sig_hat = tcrossprod(Y) 
  
  theta = c(1.5)
  iter = 1; del = 1
  while(mean(abs(del))>10^-3&iter<1000)
  {
    #cat(mean(abs(del)), " ", theta, "\n")
    del = lse.step(sig_hat, beta = theta[1], sigy2 = var(as.vector(Y)), dist_loc2^2)
    theta = theta - del*0.5
    if (any(theta<0|theta>10))
      theta = runif(1, 0.1,10)
    iter = iter+1
  }
  est_cov = infer.cov(sig_hat, beta = theta[1], sigy2 =  var(as.vector(Y)), dist_loc2^2, Time = 1)
  CI = (theta-1.96*sqrt((est_cov))<theta_true)&(theta_true<theta+1.96*sqrt((est_cov)))
  
  sig2_hat = var(as.vector(Y))
  sig = cov.quad(dist_loc2, theta[1], sig2_hat)
  sig2_cov = 2*tr(sig%*%sig)/nrow(dist_loc2)^2
  sig2_CI = (sig2_hat-1.96*sqrt((sig2_cov))<1)&(1<sig2_hat+1.96*sqrt((sig2_cov)))
  return(list(c(abs(theta), sig2_hat), c(sqrt(est_cov), sqrt(sig2_cov)), c(CI, sig2_CI)))
}

# only beta
lse.step <- function(sig_hat, beta, sigy2, dist_loc2)
{
  beta_s = beta^2+dist_loc2
  sig = sigy2*beta^2/beta_s
  
  sig_del = sig_hat - sig
  
  ### first order derivative: gradient
  first_beta = 2*beta*sigy2*dist_loc2/beta_s^2
  first_sig2 = beta^2/beta_s
  grad_beta = - 2*sum(sig_del*first_beta)
  grad_sigy2 = -2*sum(sig_del*first_sig2)
  
  
  grad_para = grad_beta #c(grad_beta, grad_sigy2)
  
  ### hessian matrix
  sec_beta = (-8*sigy2*beta^2*dist_loc2 + 2*sigy2*dist_loc2*beta_s)/beta_s^3
  sec_betasig2 = 2*beta*dist_loc2/beta_s^2
  
  
  hmat = -2*sum(sig_del*sec_beta) + 2*sum(first_beta^2)
  
  del = grad_para/hmat #grad_beta/hmat
  
  return(del)
}


infer.cov<-function(sig_hat, beta, sigy2, dist_loc2, Time)
{
  beta_s = beta^2+dist_loc2
  sig = sigy2*beta^2/beta_s
  
  sig_del = sig_hat - sig
  
  ### first order derivative: gradient
  first_beta = 2*beta*sigy2*dist_loc2/beta_s^2
  first_sig2 = beta^2/beta_s
  grad_beta = - 2*sum(sig_del*first_beta)
  grad_sigy2 = -2*sum(sig_del*first_sig2)
  
  
  grad_para = grad_beta#c(grad_beta, grad_sigy2)
  
  ### hessian matrix
  sec_beta = (-8*sigy2*beta^2*dist_loc2 + 2*sigy2*dist_loc2*beta_s)/beta_s^3
  sec_betasig2 = 2*beta*dist_loc2/beta_s^2
  
  
  hmat = -2*sum(sig_del*sec_beta) + 2*sum(first_beta^2)
  
  ### first order covariance
  W1 = first_beta
  var_grad = 2*tr(W1%*%sig%*%W1%*%sig)
  var_grad = var_grad*4/Time
  
  return(var_grad/hmat^2)
}




set.seed(1234)
#hist(cov.quad(dist_loc2, 3, 1))
Ns = c(100, 200, 400, 800)
Nrep = 500
muCI = rep(T, Nrep)
mu_est = rep(0, Nrep)
mu_sd = rep(0, Nrep)


muCI = matrix(T, ncol = Nrep, nrow = 2)
mu_est = matrix(0, ncol = Nrep, nrow = 2)
mu_sd = matrix(0, ncol = Nrep, nrow = 2)


for (i in 1:4)
{
  N = Ns[i]
  for (r in 1:Nrep)
  {
    cat(r, "\r")
    loc = simu.loc(N); dist_loc2 = as.matrix(dist(loc))
    mus = simu.mu(beta = 1, sig2 = 1, N, mu = 0, dist_loc2)
    estMu_list = lse.X(Y = mus, dist_loc2, theta_true = c(1))
    mu_est[,r] = estMu_list[[1]]
    mu_sd[,r] = (estMu_list[[2]])
    muCI[,r] = estMu_list[[3]]
  }
  cat("N: ", N,"\n",
      "Bias: ", rowMeans(mu_est - 1), 
      "SD: ", apply(mu_est,1, sd),
      "Asym SD: ", rowMeans(mu_sd), 
      "CP: ", rowMeans(muCI), "\n")
}

sd(mus)


### test the obj func

xx = seq(0.1, 5, 0.5)
yy = sapply(xx, object.func, sig2 = 1, sig_hat = tcrossprod(mus))
plot(xx, yy, type = "l")



