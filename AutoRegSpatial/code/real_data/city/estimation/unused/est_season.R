
library(showtext)
library(plyr)


rm(list = ls())

setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/")
load("rda/Ymat.rda")
load("rda/loc1.rda")


### source the functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/city/realFunc_city2.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/exp/est_infer_func.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/real_data/predictFunc2.R") # prediction functions (similar to kriging functions)

### distance matrix
Ymat0 = Ymat
dist_loc2 = as.matrix(dist(loc1[,c("lon", "lat")]))


### cut into 4 seasons

begin_date = as.Date("2015-01-01") 
end_date = as.Date("2015-12-31") 
dd = as.Date(begin_date:end_date, origin = "1970-01-01")
Ymat0 = Ymat
gr = cut(dd, as.Date(c("2015-01-01", "2015-04-01", "2015-07-01", "2015-10-01", "2015-12-31")),include.lowest = T)
gr_num = split(1:length(dd), gr)


#### model fitting

est_list = list()
p_val_list = list()
cov_types = c("exp", "quad")
lams = c(4,4)

for (k in 1:2)
{
  est = matrix(0, nrow = 5, ncol = 4)
  p_val = matrix(0, nrow = 5, ncol = 4)
  cov.type = cov_types[k]
  lam = lams[k]
  for (i in 1:4)
  {
    ind1 = gr_num[[i]]
    Ymat = Ymat0[,ind1]
    Time = ncol(Ymat)
    ### estimate rhos and mu
    rhosMu = estRhoMu(Ymat)
    # rhos_df = data.frame(citycode = gsub("PM2_5.", "", colnames(rhosMu)), rhos = rhosMu[2,])
    # hist(rhosMu[2,], xlab = "Autoregressive Coefficients", col = "dodgerblue")
    # dim(rhos_df)
    # save(rhos_df, file = "rda/rhos_df.rda")
    show(mean(rhosMu[1,])/(1-mean(rhosMu[2,])))
    
    ### for estimation theta_x
    (thetaMu = estThetaMu(rhosMu[1,], dist_loc2*lam, infer = T, p.value = T, cov.type = cov.type))
    
    (thetaX = estThetaRho(rhosMu[2,], dist_loc2*lam, infer = T, p.value = T, cov.type = cov.type))
    
    ### for estimation theta_e
    (thetaE = lse.theta(Ymat, dist_loc2*lam, rhos = rhosMu[2,], mu = rhosMu[1,], infer = T, p.value = T, cov.type = cov.type))
    
    est[,i] = c(thetaMu[[1]][1:2], thetaX[[1]][1:2], thetaE[[1]][1])
     if (k==2)
       est[c(2,4,5),i] = est[c(2,4,5),i]*1000
    p_val[,i] = c(thetaMu[[4]][1:2], thetaX[[4]][1:2], thetaE[[4]][1])
  }
  est_list[[k]] = est
  p_val_list[[k]] = p_val
}


###
para_list = list(c("$\\mu$", "$\\beta_u$", "$\\alpha$", "$\\beta_x$", "$\\beta_e$"),
                 c("$\\mu$", "$\\phi_u$", "$\\alpha$", "$\\phi_x$", "$\\phi_e$"))

for (m in 1:2)
{
  est = est_list[[m]]
  p_val = p_val_list[[m]]
  for (k in 1:5)
  {
    est_c = specify_decimal(est[k,], 2)
    pv_c = specify_decimal(p_val[k,], 3)
    
    pv_c = paste0("(", pv_c, ")")
    pv_c[p_val[k,]<0.001] = "***"
    
    #cat(para_list[[m]][k], "&", paste(est_c,  collapse = " & "), "\\\\ \n")
    cat(para_list[[m]][k], "&", paste(paste(est_c, pv_c),  collapse = " & "), "\\\\ \n")
    
  }
  cat("\\hline \n")
  cat("\\multicolumn{5}{c}{Quadratic Spatial Covariance}\\\\ \n")
}




