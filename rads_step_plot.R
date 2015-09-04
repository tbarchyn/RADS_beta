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

# rads_step_plot.R - this plots rads profiles for step outputs

print ('plotting profiles')
setwd (commandArgs(trailingOnly = T)[1])

#setwd('C:/Users/tom/Dropbox/RADS/')

get_year_from_filename <- function (filename) {
    # convenience function to get the year from a filename
    year <- as.numeric(strsplit (as.character(strsplit(filename,'.csv')), '_')[[1]][2])
    return (year)
}

plot_profile <- function (filename, xlims, ylims) {
    # function to plot profiles
    
    # filename = filename of the csv file
    # xlims = the xlims of the plot
    # ylims = the ylims of the plot
    
    # read the file
    pro <- read.table (filename, header = T, sep = ',')

    # get the id
    id <- pro$id[1]
    
    png_name <-  paste(strsplit(filename,'.csv'), '.png', sep = '')
    year <- get_year_from_filename (filename)
    
    png(png_name, width = 1000, height = 500)
    
    # plot the profile   
    plot (pro$dists, pro$bl, ylim = ylims, xlim = xlims,
          ylab = 'elevation', xlab = 'distance',
          type = 'l', col = 'red',
          main = paste ('RADS system STEP profile, id:', id, 'year:', year))
    lines (pro$dists, pro$dtm, col = 'blue')
    lines (pro$dists, pro$bsmt, col = 'black')
    dev.off()
}

file_list = list.files (pattern = '.csv')

# read in the last profile to get the plot limits
years <- 0
for (filename in file_list) {
    year <- get_year_from_filename (filename)
    years <- c(years, year)
}

max_year <- max(years)
filename <- paste ('t_', max_year, '.csv', sep = '')
pro <- read.table (filename, header = T, sep = ',')
max_y = max(pro$dtm) + 15.0
min_y = min(pro$bsmt[pro$bsmt > -1000]) # sometimes basement comes back with nas
max_x = max(pro$dists)
min_x = min(pro$dists)
ylims = c(min_y, max_y)
xlims = c(min_x, max_x)

# now, create the plots
for (filename in file_list) {
    plot_profile (filename, xlims, ylims)
}

# run quantile plots
source ('../rads_step_quants.R')
run_quants()

print ('plotting complete')


