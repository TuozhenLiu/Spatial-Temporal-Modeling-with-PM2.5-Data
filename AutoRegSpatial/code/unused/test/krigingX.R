


N = 500

for (N in c(50, 100, 200, 400))
{
  dist_loc2 = simu.dist(N)
  mm0 = rep(0, 500)
  mm = rep(0, 500)
  for (r in 1:500)
  {
    cov_X = cov.exp(dist_loc2, betaX, sig2X)
    
    #### simulate epsilon
    eig_X = eigen(cov_X)
    sqrt_value = sqrt(eig_X$values)
    X = eig_X$vectors%*%(sqrt_value*rnorm(N, 0, 1))
    rhos = 2*pnorm(X+alpha)-1
    
    train_ind = 1:(N-2)#sample(1:N, floor(2/3*N))
    pred_ind = setdiff(1:N, train_ind)
    dist_train = dist_loc2[train_ind, train_ind]
    dist_pred = dist_loc2[train_ind, pred_ind]
    rhos_train = rhos[train_ind]
    
    X1_train = qnorm((rhos_train + 1)/2); #alpha = mean(X1_train)
    X_train = X1_train - alpha
    SigX = cov.exp(dist_train, beta = betaX, sig2 = sig2X)
    cx0 = cov.exp(dist_pred, beta = betaX, sig2 = sig2X)
    
    cx0SigX = crossprod(cx0, solve(SigX))
    mu_x0 = cx0SigX%*%X_train
    sig2x0 = SigX[1,1] - diag(cx0SigX%*%cx0)
    rand_x0 = apply(cbind(mu_x0, sqrt(sig2x0)), 1, 
                    function(x) return(rnorm(n = 1000, mean = x[1], sd = x[2])))
    rhos_pred = rho.func(rand_x0, alpha)
    rhos1_pred = colMeans(rhos_pred)
    mm0 = mean(X[pred_ind]-mu_x0)
    mm[r] = mean(rhos[pred_ind] - rhos1_pred)
    cat(r, "\r")
  }
  cat("\n", mean(mm0), mean(mm), "\n")
  cat("\n", mean(mm0^2), mean(mm^2), "\n")
}



rhosk_pred = sapply(1:(Time-2), function(k) colMeans(rhos_pred^k))

mean(X[pred_ind]-mu_x0)
mean(colMeans(rand_x0) - X[pred_ind])


