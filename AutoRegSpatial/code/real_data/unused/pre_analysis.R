
library(plyr)
library(ggplot2)
library(ggmap)

setwd("E:/1. Academic/data/pm2.5data/org_data/")
dat = read.csv("Beijing.csv", stringsAsFactors = F)
loc = read.csv("station0.csv")[,1:3]
load("beijing.rda")



summary(dat)

dat_sub = dat[dat$year==2015,]
dat_sub = dat_sub[order(dat_sub$month, dat_sub$day),]
sum(dat_sub$PM2_5==-99) # 有196个缺失
dat_sub$PM2_5[dat_sub$PM2_5==-99] = NA
dat_sub = dat_sub[!is.na(dat_sub$PM2_5),]


loc1 = loc[is.element(loc$station, dat_sub$station),]

plot(c(min(loc1$lon), max(loc1$lon)), c(min(loc1$lat), max(loc1$lat)), type = "n", axes = FALSE)
text(loc1$lon, loc1$lat, gsub("Beijing_", "", loc1$station), cex = 0.8)

locs = paste0("Beijing_", c("Changping", "Fengtai", "Daxing", "Nongzhanguan"))
par(mfrow = c(4,2), mai = c(0.85,1,1,1)-0.12)
for (ll in locs)
{
  dat_sub1 = dat_sub[dat_sub$station==ll,]
  plot(dat_sub1$PM2_5, type = "l", xlab = "Day", 
       ylab = "PM 2.5", main = gsub("Beijing_", "", ll), ylim = c(0, 600), lwd = 1.5)
  pacf(dat_sub1$PM2_5, xlab = "Lag", ylab = "PACF", main = gsub("Beijing_", "", ll))
}







par(mfrow = c(1,1))
station_mean = ddply(dat_sub, .(station), function(x) pm_mean = mean(x$PM2_5))
colnames(station_mean)[2] = "pm_mean"
pm_mean = station_mean$pm_mean; names(pm_mean) = gsub("Beijing_", "", station_mean$station)
par(mai = c(2.5,1.5,1,1))
barplot(sort(pm_mean, decreasing = T), las = 3, col = heat.colors(40), ylab = "mean of PM2.5")




pm_mean_df = ddply(dat_sub, .(station), function(x) pm_mean = mean(x$PM2_5)); colnames(pm_mean_df)[2] = "pm_mean"
loc1 = loc1[match(pm_mean_df$station, loc1$station),]
loc_dist = as.matrix(dist(loc1[,2:3]))
near_stations = apply(loc_dist, 1, function(x) order(x, decreasing = F)[2])
pm_mean_df$near_pm_mean = pm_mean_df$pm_mean[near_stations]
far_stations = apply(loc_dist, 1, function(x) order(x, decreasing = T)[1])
pm_mean_df$far_pm_mean = pm_mean_df$pm_mean[far_stations]

par(mfrow = c(1,2), mai = c(1,0.6,1,0.7)+0.5)
plot(pm_mean_df$pm_mean, pm_mean_df$near_pm_mean, pch = 20, 
     xlab = "Average PM2.5", ylab = "Average PM2.5 of Nearest Stations")
plot(pm_mean_df$pm_mean, pm_mean_df$far_pm_mean, pch = 20, 
     xlab = "Average PM2.5", ylab = "Average PM2.5 of Farthest Stations")






station_mean_loc = merge(loc, station_mean, by = "station")
write.csv(station_mean_loc, "../output_data/station_mean_loc15.csv", row.names = F)



colormap = c("Blue","Green","Yellow","Red","White")  #设置热力图各级颜色

#beijing = get_map(location="Beijing",  maptype='terrain', source = "osm") #交通图


ggplot(data = station_mean_loc,
       aes(x = lon, y = lat, size = (pm_mean), colour = pm_mean), alpha = 0.5) + 
  geom_point() + scale_size(range = c(1, 15)) + 
  scale_colour_gradient(low = "white", high = "red")+
  theme(plot.margin = unit(c(0,1,1,0), "cm"),
        axis.text=element_text(size=17, family = "wqy-microhei"),
        axis.title.x = element_text(vjust=-2),
        axis.title=element_text(size=20,face="bold"))

ggmap(roadmap, extent='device') + 
  geom_point(data = station_mean_loc,
             aes(x = lon, y = lat, size = (pm_mean), colour = pm_mean)) + 
  scale_size(range = c(1, 10)) + 
  scale_colour_gradient(low = "white", high = "black")









