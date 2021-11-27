### This is the main file

# set the file path as the work path


### source the functions
source("../simulator.R") # simulator functions
source("../predictFunc2.R") # prediction functions (similar to kriging functions)
source("../est_infer_func.R")


### format specify
specify_decimal = function(x, k) gsub('\\s+','',format(round(x, k), nsmall=k))                                   ### format function for keeping 2 digits after


#Nrep = 200
cov.type1 = "exp"
cov.type2 = "quad"

### for the exponential 
### strong spatial correlation
betaE = 2; sig2E = 1
betaX = 1; sig2X = 1; alpha = 0.8
betaMu = 1; sig2Mu = 1; mu = 1


Nrep = 1000
cat("\n  For the mixedratic function, Strong spatial correlation, Estimation repeat: ", Nrep, "\n ")


source("simu_est_mixed.R")


save(estRes0, file = "../../data/simu/estmixed_strong.rda")


Nrep = 200
cat("\n  For the mixed function, Strong spatial correlation, Prediction repeat: ", Nrep, "\n ")


source("simu_predCV_mixed.R")

save(res_r2_all, file = "../../data/simu/res_r2_mixed_strong.rda")
save(res_rmse_all, file = "../../data/simu/res_rmse_mixed_strong.rda")


### weak spatial correlation

betaE = 1; sig2E = 1
betaX = 1.5; sig2X = 1; alpha = 0.8
betaMu = 1.5; sig2Mu = 1; mu = 1

Nrep = 1000
cat("\n  For the mixed function, Weak spatial correlation, Estimation repeat: ", Nrep, "\n ")


source("simu_est_mixed.R")

save(estRes0, file = "../../data/simu/estmixed_weak.rda")

Nrep = 200
cat("\n  For the mixedratic function, Weak spatial correlation, Prediction repeat: ", Nrep, "\n ")


source("simu_predCV_mixed.R")

save(res_r2_all, file = "../../data/simu/res_r2_mixed_weak.rda")
save(res_rmse_all, file = "../../data/simu/res_rmse_mixed_weak.rda")



### Output result

cat("\n \n")

betaE = 1; sig2E = 1
betaX = 1.5; sig2X = 1; alpha = 0.8
betaMu = 1.5; sig2Mu = 1; mu = 1


#theta_true = c(betaE, mu, betaMu, alpha, betaX)

theta_true = list(c(2,1,1,0.8, 1), c(1,1,1.5, 0.8, 1.5))

egs = c("strong", "weak")
for (k in 1:2)
{
  load(paste0("../../data/simu/estmixed_", egs[k], ".rda"))
  load(paste0("../../data/simu/res_r2_mixed_", egs[k], ".rda"))
  for (i in 1:2)
  {
    for (t in 1:2)
    {
      if (t==1) 
        cat(Ns[i], " & ")
      else
        cat(" & ")
      est = specify_decimal(sqrt(rowMeans((estRes0[[i]][[t]][[1]][c(1,3,4,6,7),] - theta_true[[k]])^2))*100, 1)
      CI = specify_decimal(rowMeans(estRes0[[i]][[t]][[2]])[c(1,3,4,6,7)]*100, 1)
      cc = paste(est, paste("(", CI, ")", sep = ""), sep = " ", collapse = " & ")
      r2 = paste(specify_decimal(res_r2_all[[i]][t,]*100, 1), collapse = " & ")
      cat(Times[t], "&",  cc, " &", r2, " \\\\ \n")
    }
  }
  cat('\\multicolumn{9}{c}{Example 2: Weak Spatial Dependence}\\\\ \n')
}


