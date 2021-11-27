### check the objective function

obj.exp<-function(beta, dist_loc2, sigy2, sig_hat)
{
  exp_loc = exp(-beta^2*dist_loc2)
  sig = sigy2*exp_loc
  sig_del = sig_hat - sig
  return(sum(sig_del^2))
}
obj.quad<-function(beta, dist_loc2, sigy2, sig_hat)
{
  sig = cov.quad(x = dist_loc2, sig2 = sigy2, beta = beta)
  sig_del = sig_hat - sig
  return(sum(sig_del^2))
}

### mus
vec = rhosMu[1,]

### rhos
y1 = (rhosMu[2,]+1)/2
vec = log(y1/(1-y1))

### for residuals
vec = Ymat[,-1] - rhosMu[2,]*Ymat[,-Time] - rhosMu[1,] # residuals

mu_hat = mean(vec)
Y = vec - mu_hat
sig2_hat = mean(Y^2)
sig_hat = tcrossprod(Y)
xx = seq(0, 10, 0.05)
#yy = sapply(xx, obj.exp, dist_loc2 = dist_loc2, sigy2 = sig2_hat, sig_hat = sig_hat)

yy2 = sapply(xx, obj.quad, dist_loc2 = dist_loc2, sigy2 = sig2_hat, sig_hat = sig_hat)
plot(xx, yy2, type  ="l")
xx[which.min(yy2)]


