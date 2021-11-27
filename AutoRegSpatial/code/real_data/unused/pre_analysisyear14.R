
library(xlsx)
dat14 = dat[dat$year==2014,]
dat14 = dat14[order(dat14$month, dat14$day),]
sum(dat14$PM2_5==-99) # 有196个缺失
dat14$PM2_5[dat14$PM2_5==-99] = NA
dat14$date = as.Date(paste(dat14$year, dat14$month, dat14$day, sep = "-"))

dat14_sub = dat14[dat14$date>=as.Date("2014-11-08")&dat14$date<=as.Date("2014-11-15"),]

station_1411_loc = merge(dat14_sub[,-(2:5)], loc, by = "station")
write.csv(station_1411_loc, paste0("../output_data/station_PM25_1411_loc.csv"), row.names = F)


dd = as.Date("2014-11-08")
while (dd<as.Date("2014-11-15"))
{
  df = station_1411_loc[station_1411_loc$date==dd,]
  write.csv(df, paste0("../output_data/station_PM25_1411_loc_", dd, ".csv"), row.names = F)
  dd = dd + 1
}

