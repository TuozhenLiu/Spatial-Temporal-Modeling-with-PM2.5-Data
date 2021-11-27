ce0 = cov.exp(dist_pred, beta = betaE, sig2 = sig2E)

X1_train = qnorm((rhos_train + 1)/2); #alpha = mean(X1_train)
X_train = X1_train - alpha

SigX = cov.exp(dist_train, beta = betaX, sig2 = sig2X)
cx0 = cov.exp(dist_pred, beta = betaX, sig2 = sig2X)

cx0SigX = crossprod(cx0, solve(SigX))
mu_x0 = cx0SigX%*%X_train
sig2x0 = SigX[1,1] - diag(cx0SigX%*%cx0)
sig2x0[sig2x0<0] = 10^(-5)
rand_x0 = apply(cbind(mu_x0, sqrt(sig2x0)), 1, 
                function(x) return(rnorm(n = 1000, mean = x[1], sd = x[2])))
rhos_pred = rho.func(rand_x0, alpha)
c1rho = sapply(rhos_train, function(x) colMeans(1/(1-x*rhos_pred)))
c2rho = sapply(rhos_train, function(x) colMeans(1/(1-x)/(1-rhos_pred)))
Y_pred = (t(ce0)*c1rho + mu^2*c2rho)%*%solve(SigY)%*%(Y_train[,Time] - muY) + muY
