
### This file is for:
### fix sigma to estimate spatial dependence para & mean para

setwd("E:/OneDrive/1. Academic/spatial_autoreg/code/simu/")

### format specify
specify_decimal = function(x, k) gsub('\\s+','',format(round(x, k), nsmall=k))                                   ### format function for keeping 2 digits after


source("../../real_data/city/realFunc_city2.R")
source("../../real_data/city/realFunc_city_quad.R")

### source the functions
source("../simulator.R") # real data functions
source("../est_infer_func.R") # real data functions
source("../predictFunc2.R") # prediction functions (similar to kriging functions)


### check the covriance values
loc = simu.loc(300)
dist_loc2 = as.matrix(dist(loc))
aa = cov.exp(dist_loc2, 1,1)
hist(aa[upper.tri(aa)])
quantile(aa[upper.tri(aa)])
sum(aa[upper.tri(aa)]>0.8)

### show spatial correlation levels
png("E:/OneDrive/1. Academic/spatial_autoreg/tex/spatial_cov.png", width = 800, height = 550)                           
par(mfrow = c(1,2))
xx = seq(0, 5, 0.1)
plot(xx, cov.exp(xx, 2,1), type = "l", xlab = "Distance", ylab = "Spatial Dependence", lty = 1, lwd = 2, cex = 5, cex.lab = 1.5, cex.axis = 1.5)
lines(xx, cov.exp(xx, 1,1), type = "l", col = "red", lty = 2,lwd = 2)
legend(2.5, 1, c(expression(beta == 2), expression(beta == 1)), col = c("black","red"), lty = 1:2, lwd = 2,
       y.intersp = 0.8, x.intersp = 0.5,
       cex = 1.5, bty='n')
plot(xx, cov.quad(xx, 1,1), type = "l", xlab = "Distance", ylab = "Spatial Dependence", lty = 1, lwd = 2, cex = 5, cex.lab = 1.5, cex.axis = 1.5)
lines(xx, cov.quad(xx, 2,1), type = "l", col = "red", lty = 2,lwd = 2)
legend(2.5, 1, c(expression(phi == 1), expression(phi == 2)), col = c("black","red"), lty = 1:2, lwd =  2,
       y.intersp = 0.8, x.intersp = 0.5, 
       cex = 1.5, bty='n')
dev.off()


