### 每一天单独看地理相关性
citycodes1 = intersect(loc$citycode, gsub("PM2_5.", "", colnames(dat_wide2)))
dat_wide3 = dat_wide2[is.element(gsub("PM2_5.", "", colnames(dat_wide2)), citycodes1),]
loc2 = loc[match(gsub("PM2_5.", "", colnames(dat_wide3)), loc$citycode),]
summary(loc2)
loc_dist = as.matrix(dist(loc2[,3:4]))
near_citycodes = apply(loc_dist, 1, function(x) order(x, decreasing = F)[2])
near_city_mat = dat_wide3[,near_citycodes]
far_citycodes = apply(loc_dist, 1, function(x) order(x, decreasing = T)[1])
far_city_mat = dat_wide3[,far_citycodes]

par(mfrow = c(1,1), mai = c(1,0.6,1,0.7)+0.5)
plot(log(dat_wide3[,3]), log(near_city_mat[,3]), pch = 20, 
     xlab = "Average PM2.5", ylab = "Average PM2.5 of Nearest citycodes") ### 分别location的地理相关性非常明显

### 画出每个月每天全国的boxplot
for (i in 1:12)
{
  show(i)
  png(paste0("F:/OneDrive/1. Academic/spatial_autoreg/report/boxplot_China/LogPM25_month",i,".png"), width = 800, height = 480)
  boxplot(log(dat_sub$PM2_5[dat_sub$month==i])~dat_sub$day[dat_sub$month==i])
  dev.off()
}