library(ggmap)

load("beijing.rda")
### draw the rhos map
ggmap(roadmap, extent='device') + 
  geom_point(data = station_rho_loc,
             aes(x = lon, y = lat, size = rhos, colour = rhos)) + 
  scale_size(range = c(1, 10)) + 
  scale_colour_gradient(low = "white", high = "black")

write.csv(station_rho_loc, file = "../output_data/station_rho_loc.csv", row.names = F)


ggmap(roadmap, extent='device') + 
  geom_point(data = station_eps_loc,
             aes(x = lon, y = lat, size = epsMean, colour = epsMean)) + 
  scale_size(range = c(1, 10)) + 
  scale_colour_gradient(low = "white", high = "black")
