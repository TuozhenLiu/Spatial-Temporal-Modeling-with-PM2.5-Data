
library(plyr)
library(ggplot2)
library(ggmap)

setwd("E:/1. Academic/data/pm2.5data/org_data/")
dat = read.csv("Beijing.csv", stringsAsFactors = F)
loc = read.csv("station0.csv")[,1:3]
#load("beijing.rda")


source("F:/OneDrive/1. Academic/spatial_autoreg/code/real_data/krigingMuFunc.R")
source("F:/OneDrive/1. Academic/spatial_autoreg/code/simu20161213/func.R")
source("F:/OneDrive/1. Academic/spatial_autoreg/code/real_data/realFunc.R")
source("F:/OneDrive/1. Academic/spatial_autoreg/code/real_data/predictFunc.R")


summary(dat)

dat_sub = dat[dat$year==2015,]
dat_sub = dat_sub[order(dat_sub$month, dat_sub$day),]
sum(dat_sub$PM2_5==-99) # 有196个缺失
dat_sub$PM2_5[dat_sub$PM2_5==-99] = NA

dat_sub$date = as.Date(paste(dat_sub$year, dat_sub$month, dat_sub$day, sep = "-"))

dat_wide = reshape(dat_sub[,c("station", "date", "PM2_5")], idvar = "date", timevar = "station", direction = "wide")
quantile(apply(dat_wide[,-1], 1, function(x) mean(is.na(x))))
station_na = apply(dat_wide[,-1], 2, function(x) sum(is.na(x)))

dat_wide1 = dat_wide[, c(1, which(station_na<=30)+1)] ### 生成station-date的
# sum(dat_wide1$date[-1] - dat_wide1$date[-nrow(dat_wide1)]<=0) # 检查日期是否按顺序排列


### Estimate Rho function
dat_wide2 = apply(dat_wide1[,-1], 2, function(x){ ### impute the missing values
  if (sum(is.na(x))==0)
    return(x)
  ind = which(is.na(x))
  Xmat = cbind(x[ind-3], x[ind-2], x[ind-1], x[ind+1], x[ind+2], x[ind+3])
  x[ind] = apply(Xmat, 1, mean, na.rm = T)
  return(x)
})
station_na = apply(dat_wide2, 2, function(x) sum(is.na(x)))
dat_wide2 = dat_wide2[, station_na==0]
dim(dat_wide2)



Ymat = log(t(dat_wide2))[,1:365] ### transform to Ymat
hist(Ymat, xlab = "Log PM2.5", col = "dodgerblue")
plot(colMeans(Ymat), type = "l", xlab = "Day", ylab = "Average Response", lwd = 1.5)
#Ymat = Ymat - mean(log(dat_wide2)) # 变成log之后center
#summary(t(Ymat))

mean(Ymat)
rhosMu = estRhoMuIter(Ymat)

rhos_df = data.frame(station = gsub("PM2_5.", "", colnames(rhosMu)), rhos = rhosMu[2,])
hist(rhosMu[2,], xlab = "Autoregressive Coefficients", col = "dodgerblue")

### obtain mu and rhos

(mu = mean(rhosMu[1,]))
rhos = rhosMu[2,]
station_rho_loc = merge(loc, rhos_df, by = "station")
dist_loc2 = as.matrix(dist(station_rho_loc[,2:3]*10))
#write.csv(station_rho_loc, file = "../output_data/station_rho_loc.csv", row.names = F)


### for estimation theta_x
#X1 = qnorm((rhos+1)/2)
#alpha = mean(X1)
#(thetaX = c(alpha, lse.X(X1 - alpha, dist_loc2)))
thetaX = estThetaRho(rhos, dist_loc2)

# quantile(apply(dat_wide2, 1, function(x) mean(is.na(x)))) # 检查是否还有缺失值

### filter to obtain eps
Time = ncol(Ymat)
eps = filter(Ymat, rhos = rhosMu[2,], mu = rhosMu[1,])
hist(eps, xlab = "Residuals", col = "dodgerblue")

station_eps_loc = data.frame(station_rho_loc[,-4], epsMean = eps[, Time-1])
#write.csv(station_eps_loc, 
#          file = "../output_data/station_eps_loc.csv", row.names = F)


### for estimation theta_e
(thetaE = lse.theta(Ymat, dist_loc2, rhos = rhosMu[2,], mu = rhosMu[1,]))




### Kriging surface

summary(station_rho_loc)

lattice_num = 20
loc_pred = cbind(lat = rep(seq(39.52, 40.53, length.out = lattice_num), each = lattice_num), 
                  lon = rep(seq(116, 117.1, length.out = lattice_num), lattice_num))

dist_loc2_all = as.matrix(dist(rbind(loc_pred, station_rho_loc[,2:3])*10))

loc_train = station_rho_loc[,2:3]*10
loc_pred = loc_pred*10

dist_pred = dist_loc2_all[(lattice_num^2+1):(lattice_num^2+nrow(Ymat)), 1:lattice_num^2]

(Yt_pred1 = krigingMu.typeI(Y_train = Ymat, rhos_train = rhos, 
                            loc_train = loc_train, loc_pred = loc_pred, 
                            thetaE, thetaX, mu))

(Yt_pred2 = krigingMu.typeII(Y_train = Ymat, rhos_train = rhos, 
                           dist_train = dist_loc2, dist_pred = dist_pred, 
                           thetaE, thetaX, mu = mu))
par(mfrow = c(1,3))
hist(Yt_pred1, xlab = "Type I Kriging", main = "", col = "dodgerblue", cex = 1.5)
hist(Yt_pred2, xlab = "Type II Kriging", main = "", col = "dodgerblue")
hist(Ymat[,ncol(Ymat)], xlab = "PM2.5 of 2015.12.31", main = "", col = "indianred1")


Yt_pred1_mat = trans2mat(Yt_pred1, lattice_num)
Yt_pred2_mat = trans2mat(Yt_pred2, lattice_num)

par(mfrow = c(1,1))
heatmap(exp(Yt_pred1_mat), Rowv = NA, Colv = NA, scale = "none")
heatmap(exp(Yt_pred2_mat), Rowv = NA, Colv = NA, scale = "none")

Yt_loc_pred = data.frame(loc_pred/10, Yt_pred1 = exp(Yt_pred1), Yt_pred2 = exp(Yt_pred2))

write.csv(Yt_loc_pred, "../output_data/Yt_loc_pred.csv", row.names = F)
