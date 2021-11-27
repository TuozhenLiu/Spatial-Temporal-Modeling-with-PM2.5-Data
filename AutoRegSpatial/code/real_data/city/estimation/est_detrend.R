
library(showtext)
library(plyr)


rm(list = ls())

setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/")
load("rda/Ymat.rda")
load("rda/Ymat_detrend.rda")
load("rda/loc1.rda")


### source the functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/simulator.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/est_infer_func.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/predictFunc2.R") # prediction functions (similar to kriging functions)


### distance matrix
#Ymat0 = Ymat
dist_loc2 = as.matrix(dist(loc1[,c("lon", "lat")]))

Ymat0 = Ymat
Ymat = Ymat_detrend
Time = ncol(Ymat)
### estimate rhos and mu
rhosMu = estRhoMu(Ymat)
show(mean(rhosMu[1,])/(1-mean(rhosMu[2,])))
model_name = c("{\\sc Exponential}", "{\\sc Quadratic}")
cov.types = c("exp", "quad")
lam = 1

for (i in 1:2)
{
  if (i==2)
  {
    cat(" & $\\mu$ & $\\phi_u$ & $\\alpha$ & $\\phi_x$ & $\\phi_e$ \\\\ \n")
    cat("\\hline \n")
  }
  cov.type = cov.types[i]
  ### for estimation theta_x
  (thetaMu = estThetaMu(rhosMu[1,], dist_loc2*lam, infer = T, p.value = T, cov.type = cov.type))
  (thetaX = estThetaRho(rhosMu[2,], dist_loc2*lam, infer = T, p.value = T, cov.type = cov.type))
  
  ### for estimation theta_e
  (thetaE = lse.theta(Ymat, dist_loc2*lam, rhos = rhosMu[2,], mu = rhosMu[1,], infer = T, p.value = T, cov.type = cov.type))
  
  ### output in latex format
  theta = c(thetaMu[[1]][1:2], thetaX[[1]][1:2], thetaE[[1]][1])
  theta_p = c(thetaMu[[4]][1:2], thetaX[[4]][1:2], thetaE[[4]][1])
  pv_c = specify_decimal(theta_p, 3)
  pv_c = paste0("(", pv_c, ")")
  pv_c[theta_p<0.001] = "($<0.001$)"
  
  est = specify_decimal(theta, 2)
  cat(model_name[i], "&", paste(est,  collapse = " & "), "\\\\ \n")
  cat("& ", paste(pv_c,  collapse = " & "), "\\\\ \n")
  cat("\\hline \n")
}

### draw mus and rhos

png("E:/OneDrive/1. Academic/spatial_autoreg/tex/hist_murho.png", width = 800, height = 550)    
par(mfrow = c(1,2))
hist(rhosMu[1,], xlab = "Estimated Intercept Values", main = "", col = "grey", cex.axis = 1.2, cex = 1.2, cex.lab = 1.2)
hist(rhosMu[2,], xlab = "Estimated Autoregression Coefficients", main = "", col = "grey", cex.axis = 1.2, cex = 1.2, cex.lab = 1.2)
dev.off()

rowMeans(rhosMu) # show mean values



