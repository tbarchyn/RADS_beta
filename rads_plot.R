# Reactivation Advance Dimensioning System - RADS
# Developed at the University of Calgary
# This software has no warranty whatsoever.

# Copyright Thomas E. Barchyn, 2015

# This file is part of RADS.

# RADS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# RADS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with RADS.  If not, see <http://www.gnu.org/licenses/>.


# rads_plot.R - this plots rads profiles

#setwd('C:/Users/tom/Dropbox/RADS/')

plot_profile <- function (){
    # read the file
    pro <- read.table ('profile.csv', header = T, sep = ',')

    # get the bounds and id
    max_y = max(pro$dtm) + 15.0
    min_y = min(pro$bsmt[pro$bsmt > -1000]) # sometimes basement comes back with nas
    id <- pro$id[1]
    
    # plot the profile   
    plot (pro$dists, pro$bl, ylim = c(min_y, max_y),
          ylab = 'elevation', xlab = 'distance',
          type = 'l', col = 'red',
          main = paste ('RADS system profile probe for id:', id))
    lines (pro$dists, pro$dtm, col = 'blue')
    lines (pro$dists, pro$bsmt, col = 'black')
}

png('profile.png', width = 1000, height = 500)
plot_profile()
dev.off()
