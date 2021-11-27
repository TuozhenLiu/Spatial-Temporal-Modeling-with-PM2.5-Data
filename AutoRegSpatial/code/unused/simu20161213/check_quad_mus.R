
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


simu.loc<-function(N, cov.type = "exp")
{
  X_loc = runif(N, 0, sqrt(N))
  Y_loc = runif(N, 0, sqrt(N))
  loc = cbind(X_loc, Y_loc)
  if (cov.type=="quad")
    loc = loc*sqrt(10)
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
  return(list(mu_cov, CI))
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
  mu_hat = mean(Y)
  Y = Y - mu_hat
  sig_hat = tcrossprod(Y) 
  
  theta = c(1.5)
  iter = 1; del = 1
  while(mean(abs(del))>10^-3&iter<1000)
  {
    #cat(mean(abs(del)), " ", theta, "\n")
    del = lse.step(sig_hat, beta = theta, sigy2 = theta_true[3], dist_loc2^2)
    theta = theta - del*0.5
    if (any(theta<0|theta>10))
      theta = runif(1, 0.1,10)
    iter = iter+1
  }
  est_cov = infer.cov(sig_hat, beta = theta, sigy2 =  theta_true[3], dist_loc2^2, Time = 1)
  CI = (theta-1.96*sqrt((est_cov))<theta_true[2])&(theta_true[2]<theta+1.96*sqrt((est_cov)))
  
  mu_list = infer.mu(thetaMu = c(mu_hat, theta, theta_true[2]), 
                     mu_true = theta_true[1], dist_loc2, Time = 1, cov.type = "quad")
  
  sig2_hat = var(as.vector(Y))
  sig = cov.quad(dist_loc2, theta, sig2_hat)
  sig2_cov = 2*tr(sig%*%sig)/nrow(dist_loc2)^2
  sig2_CI = (sig2_hat-1.96*sqrt((sig2_cov))<1)&(1<sig2_hat+1.96*sqrt((sig2_cov)))
  return(list(c(mu_hat, abs(theta), sig2_hat), c(sqrt(mu_list[[1]]), sqrt(est_cov), sqrt(sig2_cov)), 
              c(mu_list[[2]], CI, sig2_CI)))
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


mu = 1; betaMu = 1; sig2Mu = 1

set.seed(1234)
#hist(cov.quad(dist_loc2, 3, 1))
Ns = c(100, 200, 400, 800)
Nrep = 500
muCI = rep(T, Nrep)
mu_est = rep(0, Nrep)
mu_sd = rep(0, Nrep)


muCI = matrix(T, ncol = Nrep, nrow = 3)
mu_est = matrix(0, ncol = Nrep, nrow = 3)
mu_sd = matrix(0, ncol = Nrep, nrow = 3)


for (i in 1:4)
{
  N = Ns[i]
  for (r in 1:Nrep)
  {
    cat(r, "\r")
    loc = simu.loc(N,  cov.type = "quad")*2; dist_loc2 = as.matrix(dist(loc))
    mus = simu.mu(beta = betaMu, sig2 = sig2Mu, N, mu = mu, dist_loc2)
    estMu_list = lse.X(Y = mus, dist_loc2, theta_true = c(mu, betaMu, sig2Mu))
    mu_est[,r] = estMu_list[[1]]
    mu_sd[,r] = (estMu_list[[2]])
    muCI[,r] = estMu_list[[3]]
  }
  cat("N: ", N,"\n",
      "Bias: ", rowMeans(mu_est - c(mu, betaMu, sig2Mu)), 
      "SD: ", apply(mu_est,1, sd),
      "Asym SD: ", rowMeans(mu_sd), 
      "CP: ", rowMeans(muCI), "\n")
}

sd(mus)


### test the obj func

xx = seq(0.1, 5, 0.5)
yy = sapply(xx, object.func, sig2 = 1, sig_hat = tcrossprod(mus))
plot(xx, yy, type = "l")

######################################################################



setwd("E:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/exp/")

### source the functions
source("../../real_data/city/realFunc_city2.R") # real data functions
source("../../real_data/predictFunc2.R") # prediction functions (similar to kriging functions)
source("est_infer_func.R")

### format specify
specify_decimal = function(x, k) gsub('\\s+','',format(round(x, k), nsmall=k))                                   ### format function for keeping 2 digits after


betaE = 1; sig2E = 1
betaX = 1; sig2X = 1; alpha = 0.8
betaMu = 1; sig2Mu = 1; mu = 1
cov.type = "quad"

#mu = 0.8; betaMu = 1; sig2Mu = 1

set.seed(123)
#hist(cov.quad(dist_loc2, 3, 1))
Ns = c(200, 300)
Times = c(300, 500)
Nrep = 200
muCI = rep(T, Nrep)
mu_est = rep(0, Nrep)
mu_sd = rep(0, Nrep)


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


for (i in 1:2)
{
  N = Ns[i]
  for (t in 1:2)
  {
    Time = Times[t]
    for (r in 1:Nrep)
    {
      cat(r, "\r")
      loc = simu.loc(N,  cov.type = "quad")*2; dist_loc2 = as.matrix(dist(loc))
      rhos = simu.rho(beta = betaX, sig2 = sig2X, N, alpha = alpha, dist_loc2, cov.type = cov.type)
      mus = simu.mu(beta = betaMu, sig2 = sig2Mu, N, mu = mu, dist_loc2, cov.type = cov.type)
      Y = simu.Y(beta = betaE, sigy2 = sig2E, N = N, Time = Time, rhos = rhos, 
                 dist_loc2, mu = mus, cov.type = cov.type)
      rhosMu = estRhoMu(Y = Y)
      mus_est = rhosMu[1,]
      rhos_est = rhosMu[2,]
      
      estMu_list = estThetaMu(mus_est, dist_loc2, infer = T, 
                              theta_true = c(mu, betaMu, sig2Mu), cov.type = "quad")
      estrho_list = estThetaRho(rhos_est, dist_loc2, infer = T, 
                                theta_true = c(alpha, betaX, sig2X), cov.type = "quad")
      
      estE_list = lse.theta(Y, dist_loc2, rhos = as.vector(rhos_est), mu = as.vector(mus_est), 
                            theta_true = c(betaE, sig2E), cov.type = "quad")
      
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
                              theta_true = c(mu, betaMu, sig2Mu), cov.type = "quad")
      estrho_list = estThetaRho(rhos, dist_loc2, infer = T, 
                                theta_true = c(alpha, betaX, sig2X), cov.type = "quad")
      
      estE_list = lse.theta(Y, dist_loc2, rhos = as.vector(rhos), mu = as.vector(mus), 
                            theta_true = c(betaE, sig2E), cov.type = "quad")
      
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
    cat("N: ", N, " Time: ", Time, "\n",
        "Bias: ", rowMeans(mu_est - c(mu, betaMu, sig2Mu)), 
        rowMeans(rho_est - c(alpha, betaX, sig2X)),
        rowMeans(E_est - c(betaE, sig2E)), "\n",
        "True Bias: ", rowMeans(muT_est - c(mu, betaMu, sig2Mu)), 
        rowMeans(rhoT_est - c(alpha, betaX, sig2X)),
        rowMeans(ET_est - c(betaE, sig2E)), "\n",
        "SD: ", apply(mu_est,1, sd),
        apply(rho_est,1, sd),
        apply(E_est,1, sd), "\n",
        "True SD: ", apply(muT_est,1, sd),
        apply(rhoT_est,1, sd),
        apply(ET_est,1, sd), "\n",
        #"Asym SD: ", rowMeans(mu_sd), 
        "CP: ", rowMeans(muCI), rowMeans(rhoCI), rowMeans(ECI), "\n",
        "True CP: ", rowMeans(muTCI), rowMeans(rhoTCI), rowMeans(ETCI), "\n")
  }
}
