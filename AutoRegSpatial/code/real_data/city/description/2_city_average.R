require(leafletCN)
setwd("E:/OneDrive/1. Academic/spatial_autoreg/data/")

tb = read.csv('data_wide1.csv') # Switch the data in notepad++ to ANSI encoding first
tb = tb[,-1]
nms = names(tb)
nms = gsub('PM2_5\\.','',nms)
pm2.5 = colMeans(tb, na.rm = TRUE)

dat = data.frame(city = nms, PM2.5 = pm2.5)
dat$city = as.character(dat$city)
rownames(dat) = NULL

#dat$city = as.character(dat$city)
geojsonMap(dat, "city", ~city, ~PM2.5,
           popup =  paste0(dat$city, ":", dat$PM2.5),
           palette = "Reds", legendTitle = "PM 2.5")


geojsonMap(dat, "city", ~city, ~PM2.5,na.color ="white",
           popup =  paste0(dat$city, ":", dat$PM2.5),
           palette = "Greys", legendTitle = "PM 2.5")

