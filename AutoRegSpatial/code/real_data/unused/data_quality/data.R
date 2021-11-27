setwd("E:/1. Academic/data/pm2.5data/org_data/China_2015/rda/")

### load数据
load("data_wide1.rda") # 291 个城市 365天的 PM2.5记录
write.csv(dat_wide1, "../data_wide1.csv", row.names = F)


load("loc1.rda") # 对应这291个城市的经纬度数据
write.csv(loc1, "../loc1.csv", row.names = F)

dist_loc2 = as.matrix(dist(loc1[,c("lon", "lat")])) # 计算291个城市之间的

### 写出回归的code