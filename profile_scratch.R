# profile scratch sheet - develop blowout profiles here

# set fresh profile
setwd('C:/Users/tom/Dropbox/RADS/')
pro <- read.table ('profile.csv', header = T, sep = ',')

pro$bl2 <- NA
brink_dist <- 100
bottom_elev = pro$dtm[pro$dists == 0.0] - 4.0

#####################################################
# upwind profile
subs = (pro$dists <= 0.0)
upwind_power = 2.0
upwind_slope = 0.005
pro$bl2[subs] = bottom_elev + (upwind_slope * (-1.0 * pro$dists[subs])^upwind_power)

# stoss profile
subs = (pro$dists > 0.0) & (pro$dists < brink_dist)
stoss_power = 1.3
stoss_slope = 0.039
pro$bl2[subs] = bottom_elev + (stoss_slope * pro$dists[subs]^stoss_power)

# slipface profile
slipface_jump = -0.6745085168424266
brink_elev = bottom_elev + (stoss_slope * brink_dist^stoss_power)
subs = (pro$dists >= brink_dist)
pro$bl2[subs] = brink_elev + ((pro$dists[subs] - brink_dist) * slipface_jump)        

# print stoss slope
(pro$bl2[pro$dists == 100] - pro$bl2[pro$dists == 0.0]) / brink_dist



#####################################################
# PLOT!
# get the bounds and id
max_y = max(pro$dtm) + 15.0
min_y = min(pro$bsmt)
plot (pro$dists, pro$bl2, ylim = c(min_y, max_y),
      ylab = 'elevation', xlab = 'distance',
      type = 'l', col = 'red')
lines (pro$dists, pro$dtm, col = 'blue')
lines (pro$dists, pro$bsmt, col = 'black')
