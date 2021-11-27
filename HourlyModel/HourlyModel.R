#小时数据插补
#边做边ppt
#1.小时数据梳理，汇报城市、站点、缺失情况（两三页ppt）
#2.复现北京代码，预测与插补
#3.分头做其他城市数据，写报告

library(tidyverse)
setwd("D:/Desktop/ppt")

# 1. 小时数据梳理

# load data 
dat <- read.csv("hourlydat/city_pol_h.csv", stringsAsFactors = F, sep = ";") %>%
  select(station_code, time, pm2_5) %>%
  # delete error obs.
  filter(station_code %in% unique(station_code[1:500000])) %>%
  mutate(time = as.Date(time), pm2_5 = as.numeric(pm2_5)) %>%
  # abnormal values to NA
  mutate(pm2_5 = ifelse(pm2_5 <= 0 | pm2_5>=1000, NA, pm2_5))

b = dat[which(is.na(dat$pm2_5)), ]

# 分月份缺失率：201906-201908连续三个月缺失，201912-202006缺失严重
month_obs = dat %>%
  mutate(month = format(time, "%Y%m")) %>%
  group_by(month) %>%
  summarise(total_obs = n())

# 分站点缺失率：G110009、G130104缺失严重，G110008较严重
station_obs = dat %>%
  filter(time < '2019-06-01' | (time >= '2019-09-01' & time < '2019-12-01')) %>%
  group_by(station_code) %>%
  summarise(total_obs = n()) %>%
  mutate(missing_rate = 1-total_obs/((3*365-31*3-30)*24))

day = c(31,28,31,30,31,30,31,31,30,31,30,31)
get_day = function(month){
  y = as.numeric(substr(month, 1, 4))
  m = as.numeric(substr(month, 6, 7))
  if(m == 2){
    if(y == 2020) 
      return(29)
    else 
      return(28)
  }
  else
    return(day[m])
}

obs = dat %>%
  na.omit()%>%
  mutate(month = format(time, "%Y-%m")) %>%
  group_by(station_code, month)%>%
  summarise(obs = n()) %>%
  mutate(total_obs = sapply(month, get_day)*24, missing_rate = 1-obs/total_obs)

#obs2 = obs %>%
#  pivot_wider(names_from = station_code, values_from = total_obs)

#stations = unique(dat$station_code)
#months = format(as.Date('2017-01-01')%m+%months(0:41), "%Y-%m")

station_list = read.csv('station_code.csv')
station_name = station_list$station_name[c(6:18, 49:57)]

obs %>%
  select(station_code, month, missing_rate)%>%
  ungroup()%>%
  add_row(station_code = c(rep('G110001', 4), 'G110009', 'G130104'), 
          month = c(paste0('2019-0', 6:8), rep('2020-01', 3)), 
          missing_rate = rep(1, 6))%>%
  ggplot(aes(x = station_code, y = month, fill = missing_rate)) +
  geom_tile() +
  scale_x_discrete(labels = station_name)+
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_gradient(name = "missing_rate", low = "steelblue", high ="white") +
  labs(x='', y='', fill = 'Missing Rate') +
  theme(axis.text.y = element_text(face = "bold",size=10),
        axis.title.y = element_text(face = "bold",size=12),
        axis.text.x = element_text(angle = 45,face = "bold", 
                                   size=10, vjust = 0.5, 
                                   colour = c(rep('black', 13), rep('brown', 9))),
        panel.background = element_rect(fill = "transparent"),
        panel.border=element_rect(fill='transparent', 
                                  color='transparent'),
        axis.line = element_line(color = "black"))

a = obs%>%
  mutate(month = as.numeric(gsub('-', '', month)))%>%
  filter(month<201902)%>%
  filter(station_code != 'G110008', 
         station_code != 'G110009',
         station_code != 'G130104')

#乱码
dat2 <- read.csv("hourlydat/city_pol_h.csv", stringsAsFactors = F, sep = ";")[500000:522103,]

## Visualizations

all_stations <- read.csv("hourlydat/station_location.csv")
names(all_stations)[1] = 'station_code'
stations_loc <- all_stations %>%
  filter(station_code %in% paste0('G13010', 1:9))
leaflet(stations_loc) %>% addTiles() %>%
  addCircleMarkers(lng=~lon, lat=~lat)

#########################
# 北京连续
dat_bj = read.csv("hourlydat/city_pol_h.csv", stringsAsFactors = F, sep = ";") %>%
  filter(!station_code %in% paste0('G13010', 1:9))%>%
  # delete error obs.
  filter(station_code %in% unique(station_code[1:250000]))%>%
  filter(time<"2020-01-01")%>%
  select(-id)

dat_bj_long = gather(dat_bj, polluant, concentration, -station_code, -time)%>%
  mutate(concentration = as.numeric(concentration))%>%
  mutate(concentration = ifelse(concentration<=0 | concentration>2000, NA, concentration))

dat_bj_wide = pivot_wider(dat_bj_long, names_from = polluant, values_from = concentration)%>%
  arrange(station_code, time)

x = diff(dat_bj_wide$pm2_5)

xx = ifelse(x==0, 1, 0)

s = 0
ss = 0
for(i in 1:length(xx)){
  if(!is.na(xx[i]) & xx[i] == 1){
    s = s+1
  }
  else{
    ss = max(ss, s)
    s = 0
  }
}

for(i in 1:length(xx)){
  if(!is.na(xx[i]) & xx[i] == 1){
    s = s+1
  }
  else{
    ss = max(ss, s)
    s = 0
  }
  if(ss == 10)
    break
}

ss

which(duplicated(na.omit(dat_bj_wide[,-c(1,2)])))


so2 = diff(as.numeric(dat_bj$so2))
no2 = diff(as.numeric(dat_bj$no2))
no2 = diff(as.numeric(dat_bj$no2))

###################################################################
# seasonal model

# predict plot
load('pred-bj-ori.RData')
pre_data = rbind(pred[[1]], pred[[2]], pred[[3]], pred[[4]])

###
pre_data_win = pred[[4]]
date_win = dat_wide$time[which(dat_wide$season=='Winter')]
pre_data_win %>%
  select(pred_ori, true_ori)%>%
  `colnames<-`(c('pred', 'true'))%>%
  mutate(date=date_win[25:length(date_win)]) %>%
  gather(key = Type, value = pm2_5, -date) %>%
  ggplot(aes(x = date, y = pm2_5, color = Type)) +
  geom_line() +
  scale_colour_manual(values=c("goldenrod4", "cyan")) +
  labs(x = '', 
       y = "PM2.5", 
       title = 'Predicted Kriging(Beijing, 2018.11-2019.02)')+
  theme_classic()
###
spr = pred[[1]];sum = pred[[2]];fal = pred[[3]];win = pred[[4]]

mean(spr$true_ori)
mean(sum$true_ori)
mean(fal$true_ori)
mean(win$true_ori)

sum(abs(spr$pred_ori-spr$true_ori)>30)/length(spr$pred_ori)
sum(abs(sum$pred_ori-sum$true_ori)>30)/length(sum$pred_ori)
sum(abs(fal$pred_ori-fal$true_ori)>30)/length(fal$pred_ori)
sum(abs(win$pred_ori-win$true_ori)>30)/length(win$pred_ori)

sum(abs(spr$pred_ori-spr$true_ori)/spr$true_ori>0.8)/length(spr$pred_ori)
sum(abs(sum$pred_ori-sum$true_ori)/sum$true_ori>0.8)/length(sum$pred_ori)
sum(abs(fal$pred_ori-fal$true_ori)/fal$true_ori>0.8)/length(fal$pred_ori)
sum(abs(win$pred_ori-win$true_ori)/win$true_ori>0.8)/length(win$pred_ori)

###############################
load('pred-sjz-ori.RData')
pre_data = rbind(pred[[1]], pred[[2]], pred[[3]], pred[[4]])

pre_data_win = pred[[1]]
date_win = dat_wide$time[which(dat_wide$season=='Spring')]
pre_data_win %>%
  select(pred_ori, true_ori)%>%
  `colnames<-`(c('pred', 'true'))%>%
  mutate(date=date_win[25:length(date_win)]) %>%
  gather(key = Type, value = pm2_5, -date) %>%
  ggplot(aes(x = date, y = pm2_5, color = Type)) +
  geom_line() +
  scale_colour_manual(values=c("goldenrod4", "cyan")) +
  labs(x = '', 
       y = "PM2.5", 
       title = 'Predicted Kriging(Shijiazhuang, 2018.02-2018.04)')+
  theme_classic()
