# 本文件用于intro部分描述分析（地理相关性+时间相关性）

library(plyr)
library(ggplot2)
library(ggmap)

setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/")

load("rda/data_wide2.rda")
dat = read.csv("data15.csv", stringsAsFactors = F, sep = ";")
load("rda/loc.rda")
load("rda/dat_sub.rda")



### 地理相关性

# 求城市PM2.5均值
pm_mean_df = ddply(dat_sub, .(citycode), function(x) pm_mean = mean(x$PM2_5, na.rm = T)); colnames(pm_mean_df)[2] = "pm_mean"
citycodes = intersect(pm_mean_df$citycode, loc$citycode)
pm_mean_df = pm_mean_df[is.element(pm_mean_df$citycode, citycodes),]
loc1 = loc[match(pm_mean_df$citycode, loc$citycode),]
summary(loc1)

# 合并最近城市的PM2.5
loc_dist = as.matrix(dist(loc1[,3:4]))
near_citycodes = apply(loc_dist, 1, function(x) order(x, decreasing = F)[2])
pm_mean_df$near_pm_mean = pm_mean_df$pm_mean[near_citycodes]

# 写入csv 方便excel powermap作图
pm_mean_df1 = merge(pm_mean_df, loc1, by = "citycode")
write.csv(pm_mean_df1[,c("lon", "lat", "pm_mean")] , file = "city_pm_mean.csv", row.names = F)

# 地理相关性展示
par(mfrow = c(1,1), mai = c(1,0.6,1,0.7)+0.5)
plot(pm_mean_df$pm_mean, pm_mean_df$near_pm_mean, pch = 20, 
     xlab = "Average PM2.5", ylab = "Average PM2.5 of Nearest citycodes")


### 时序相关性
library(showtext)
X11()
nn = c("Beijing", "Shanghai", "Xi'an", "Chengdu")
data_wide2_sub = dat_wide2[,grep("(北京)|(上海)|(西安)|(成都)", colnames(dat_wide2))]

pdf("E:/OneDrive/1. Academic/spatial_autoreg/tex/ts.pdf", onefile = T, width = 8, height = 11)
par(mfrow = c(4,2), mai = c(0.85,1,1,1)-0.12)
for (i in 1:4)
{
  plot(data_wide2_sub[,i], type = "l", xlab = "Day", 
       ylab = "PM 2.5", main = nn[i], ylim = c(0, 500), lwd = 1.5)
  pacf(data_wide2_sub[,i], xlab = "Lag", ylab = "PACF", main = nn[i])
}
dev.off()


