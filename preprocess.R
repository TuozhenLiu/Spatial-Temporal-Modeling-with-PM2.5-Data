library(tidyverse)
setwd("~/Desktop/大气")

dat <- read.csv('pollution_day.csv') %>% select(-id, -city_en)
dat_match <- read.xlsx('省市对照.xlsx')



map_dbl(dat, ~sum(is.na(.))) # 原有数据集没有缺失值

dat %>%
  mutate(year = substr(date, 1, 4),
         month = substr(date, 6, 7)) %>%
  group_by(city, year, month) %>%
  summarise(观测数量 = n()) -> temp #

temp %>%
  group_by(year, month) %>%
  summarise(月度城市数量 = length(unique(city))) %>%
  arrange(year, month)-> temp2

# 取2019-08的73个城市
dat_1 <- temp %>% filter(year == '2019', month == '08') %>%
  ungroup() %>% select(city) %>% unique() %>%
  left_join(dat %>% mutate(year = substr(date, 1, 4), month = substr(date, 6, 7)))

dat_1 %>%
  mutate(year = substr(date, 1, 4),
         month = substr(date, 6, 7)) %>%
  group_by(city, year, month) %>%
  summarise(观测数量 = n()) -> temp_1 

temp_1 %>%
  group_by(year, month) %>%
  summarise(月度城市数量 = length(unique(city))) %>%
  arrange(year, month)-> temp2_1

# 取2015-02的72个城市
dat_2 <- temp_1 %>% filter(year == '2015', month == '02') %>%
  ungroup() %>% select(city) %>% unique() %>%
  left_join(dat_1 %>% mutate(year = substr(date, 1, 4), month = substr(date, 6, 7)))

unique(dat_2$city)

dat_2 %>%
  mutate(year = substr(date, 1, 4),
         month = substr(date, 6, 7)) %>%
  group_by(city, year, month) %>%
  summarise(观测数量 = n()) -> temp_2

temp_2 %>%
  group_by(year, month) %>%
  summarise(月度城市数量 = length(unique(city))) %>%
  arrange(year, month)-> temp2_2

dat_2 %>%
  mutate(year = as.numeric(year), month = as.numeric(month)) %>%
  filter(year>2014) -> dat_new

# 补齐date
dat_new %>% select(city) %>% unique() %>%
  expand_grid(dat_new %>% select(date, year, month) %>% unique()) %>%
  left_join(dat_new) -> dat_new

map_dbl(dat_new, ~sum(is.na(.))) # 总缺失数

#load('数据和报告/2020数据/数据/data_new.RData')

dat_new %>%
  left_join(dat_match %>% rename(province = 省, city = 经营地) %>% 
              mutate(city = ifelse(str_detect(city, '市$'), str_extract(city, '.*(?=市)'), city))) -> temp

map_dbl(temp, ~sum(is.na(.)))

temp %>%
  filter(is.na(province)) %>%
  select(province, city) %>%
  unique() %>% view()

temp$province[which(temp$city %in% c('鄂尔多斯', '锡林郭勒'))] = '内蒙古省'
temp$province[which(temp$city == '延边')] = '吉林省'

dat_new <- temp

# 补时间之前的原来观测数
dat_new %>%
  filter(!is.na(o3)) %>% group_by(year, city) %>%
  summarise(n = n()) %>%
  pivot_wider(city, names_from = year, values_from = n) %>%
  write.xlsx('六地区原始数据观察数_城市.xlsx')

# 缺失率
dat_new %>%
  group_by(province, city, year, month) %>%
  summarise(total_amount = n(), valid_amount = sum(!is.na(o3)), 
            missing_rate = 1 - valid_amount / total_amount) -> dat_1

# 变异系数
dat_new %>%
  select(-c(date, quality)) %>%
  group_by(province, city, year, month) %>%
  summarise_all(list(avg = ~mean(., na.rm = T), std = ~sd(., na.rm = T), 
                cv = ~sd(., na.rm = T) / mean(., na.rm = T))) %>%
  full_join(dat_1)-> dat_month

save(dat_month, file = '72city_2015-2020数据.RData')
