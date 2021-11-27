
setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/")
load("rda/Ymat.rda")

### 去除季节以及趋势特征

library(TTR)

pm_ts = colMeans(Ymat)
pm_ts = ts(colMeans(Ymat),start=c(1))
plot(pm_ts)

pm_sma = SMA(pm_ts,n=30)
#plot.ts(pm_sma)
pm_sma = as.numeric(pm_sma)

Ymat_detrend = t(apply(Ymat[,!is.na(pm_sma)], 1, function(x) x-pm_sma[!is.na(pm_sma)]))
#plot(colMeans(Ymat_detrend), type = "l")

save(Ymat_detrend, file = "rda/Ymat_detrend.rda")


### the mean PM2.5 over the year
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
begin_date = as.Date("2015-01-01") 
end_date = as.Date("2015-12-31") 
dd = as.Date(begin_date:end_date, origin = "1970-01-01")

png("E:/OneDrive/1. Academic/spatial_autoreg/tex/hist_avePM25.png", width = 800, height = 550)    
par(mfrow = c(1,2))
hist(Ymat, xlab = "Log PM2.5", main = "", col = "grey", cex.axis = 1.2, cex = 1.2, cex.lab = 1.2)
plot(dd, colMeans(Ymat), type = "l", xlab = "Day", ylab = "Average Response", 
     lwd = 1.2, cex.axis = 1.2, cex = 1.2, cex.lab = 1.2)
lines(dd, pm_sma, col = "red", lwd = 2, cex.axis = 1.2, cex = 1.2, cex.lab = 1.2)
dev.off()

