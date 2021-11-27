### 本文件进行数据预处理及描述分析

library(showtext)
library(plyr)

setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/")
dat = read.csv("data15.csv", stringsAsFactors = F, sep = ";") # 读入全国数据集
loc = read.csv("location15.csv") # 读入经纬度数据
save(loc, file = "rda/loc.rda")

### source the functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/simulator.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/est_infer_func.R") # real data functions
source("E:/OneDrive/1. Academic/spatial_autoreg/code/predictFunc2.R") # prediction functions (similar to kriging functions)

### 保留相关列
dat_sub = dat[, c("citycode", "cityname", "months", "days", "pm25", "DATES")]
colnames(dat_sub) = c("citycode", "cityname", "month", "day", "PM2_5", "date")
dat_sub$date = as.Date(dat_sub$date) # 转换日期格式


### 城市+月+日
# aa = dat_wide[, grep("杭州", colnames(dat_wide))] # 出现了重复
# table(apply(aa, 1, function(x) sum(is.na(x))))

dat_sub = dat_sub[order(dat_sub$citycode, dat_sub$month, dat_sub$day),] # 按照城市名称排序
dat_sub$cityname[dat_sub$cityname=="杭州"] = "杭州市"
dat_sub$cityname[dat_sub$cityname=="台州"] = "台州市"

sum(dat_sub$PM2_5==-1) # 有940个缺失
dat_sub$PM2_5[dat_sub$PM2_5==-1] = NA # 将取值为-1的标识为-1
save(dat_sub, file = "rda/dat_sub.rda")

### 长表变宽表
dat_wide = reshape(dat_sub[,c("cityname", "date", "PM2_5")], 
                   idvar = "date", timevar = "cityname", direction = "wide")

quantile(apply(dat_wide[,-1], 1, function(x) mean(is.na(x)))) # 日期缺失
station_na = apply(dat_wide[,-1], 2, function(x) sum(is.na(x))) # 站点缺失
sort(station_na, decreasing = T)[1:5]

### 查看站点缺失分布
X11()
hist(station_na, xlab = "缺失天数", col = "grey")
station_na[which.max(station_na)]

# # the citynames used in dat_wide 
# city_names1 = gsub("PM2_5.", "", colnames(dat_wide[,-1]))
# loc0 = loc[match(city_names1, loc$cityname),]
# save(dat_wide, file = "rda/data_wide.rda")
# save(loc0, file = "rda/loc0.rda")

### 取缺失<=5的站点
dat_wide1 = dat_wide[, c(1, which(station_na<=5)+1)] ### 保留缺失值小于5的站点，生成station-date的宽表数据
# sum(dat_wide1$date[-1] - dat_wide1$date[-nrow(dat_wide1)]<=0) # 检查日期是否按顺序排列
city_names1 = gsub("PM2_5.", "", colnames(dat_wide1[,-1]))
loc1 = loc[match(city_names1, loc$cityname),]
save(dat_wide1, file = "rda/data_wide1.rda")
save(loc1, file = "rda/loc1.rda")
sum(is.na(loc1))


### 此时还有一部分缺失值，需要进一步插补，此处取前3天&后3天的数据
dat_wide2 = apply(dat_wide1[,-1], 2, function(x){ ### impute the missing values
  if (sum(is.na(x))==0)
    return(x)
  ind = which(is.na(x))
  indm = cbind(ind-3, ind-2, ind-1, ind+1, ind+2, ind+3)
  indm[indm<1] = 1
  Xmat = cbind(x[indm[,1]], x[indm[,2]], x[indm[,3]], 
               x[indm[,4]], x[indm[,5]], x[indm[,6]])
  x[ind] = apply(Xmat, 1, mean, na.rm = T)
  return(x)
}) ### 291个站点
station_na = apply(dat_wide2, 2, function(x) sum(is.na(x))) # 填补完之后，如果再出现缺失值的站点代表该站点连续7天有缺失值，应该删去该站点
dat_wide2 = dat_wide2[, station_na==0]
dim(dat_wide2) # dat_wide2 为最终站点
save(dat_wide2, file = "rda/data_wide2.rda")

### response matrix: log(PM2.5)
Ymat = log(t(dat_wide2))[,1:365] ### transform to Ymat
hist(Ymat, xlab = "Log PM2.5", col = "dodgerblue")


### 保留response 矩阵
cityname = gsub("PM2_5.", "", rownames(Ymat)) # city names list
mean(Ymat)
Ymat = Ymat[is.element(cityname, loc$cityname),] # 匹配有经纬度信息的
dim(Ymat)
save(Ymat, file = "rda/Ymat.rda")


