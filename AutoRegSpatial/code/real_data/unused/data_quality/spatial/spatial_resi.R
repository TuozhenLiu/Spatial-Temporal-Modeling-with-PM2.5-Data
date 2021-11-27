load("rda/Yt_krigs.rda")
begin_ind = 1; end_ind = 365

pred_inds = grep("(北京)|(天津)|(上海)|(南京)|(保定)|(石家庄)|(成都)|(广州)|(深圳)", loc1$cityname)
resi = Ymat[pred_inds, begin_ind:end_ind] - Yt_krigs[pred_inds,]

cities = as.character(loc1$cityname[pred_inds]) #c('北京', '天津', '上海', '南京', '保定', '石家庄', '成都', '广州', '深圳')
par(mfrow = c(1,2))

for (i in 1:9)
{
  plot(Yt_krigs[pred_inds,][i,], resi[i,], xlab = "Kriging 值", 
       ylab = "Kriging Residuals", main = cities[i])
  
  cat(cities[i], " : ",  mean(resi[i,]),  # residual的均值
      mean(resi[i,]<0), # 预报值低于kriging值的比例
      1-mean(resi[i,]^2)/var(Ymat[pred_inds, ][i,]), "\n ") 
  
}

city_fac = factor(rep(cities, each = ncol(resi)), 
                     levels = cities)
df = data.frame(city = city_fac,
                residual = as.vector(t(resi)))
op = par(bg = "#EFEFEF")
boxplot(residual~city_fac, data = df, col = rep(c("#F38630", "indianred","#4ECDC4", "#556270","#CDB380"), each = 2))
