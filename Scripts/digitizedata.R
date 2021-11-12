# Get data from Paper: Mieczan et al. 2015

library(digitize)
spring = digitize("Mieczan2015a.PNG", x1 = -32, x2 = -20, y1 = -9, y2 = 1)

spring

summer = digitize("Mieczan2015b.PNG", x1 = -32, x2 = -20, y1 = -9, y2 = 1)

summer

fall = digitize("Mieczan2015c.PNG", x1 = -32, x2 = -20, y1 = -9, y2 = 1)

# Assign point IDS

spring$ID = c("Sph_ang", "Phycof","OMB", "Rotif", "Copep", "Hya_pap")
spring$season = "spring"

summer$ID = c("Sph_ang", "Phycof","OMB", "Copep", "Rotif", "Hya_pap")
summer$season = "summer"

fall$ID = c("Sph_ang", "Phycof","OMB", "Rotif", "Copep", "Hya_pap")
fall$season = "fall"

rbind(spring, summer, fall) %>% write_rds("Data/Mieczan2015.rds")
