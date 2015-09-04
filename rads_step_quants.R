# RADS step output analysis script
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


# PURPOSE: the idea here is to run this to calculate the range of outcomes
# from a given point against the range of empirical dhdt results.

# this should wrap up all the respective variables (flux, vpeak, variability
# into a basic idea of probabilities.


run_quants <- function () {
    # function to get the current empirical quantile distribution and throw some
    # advance quantile plots into the directory

    # read in the empirical quantiles, and check the read
    quants <- read.table ('D:/consulting/mbss/model/dhdt/climate_8/dhdt_quantiles.csv', header = T, sep = ',')
    quants$probs <- seq (0,1,0.01)

    # read in the step output
    steps <-  read.table ('../rads_step_output.csv', header = T, sep = ',')

    # get the kick age and brink dists
    quants$vegtime <- NA
    quants$brink_dist <- NA
    
    for (i in 1:nrow(quants)) {
        # loop down through the steps and find the first request year that corresponds
        # with a deporate that is below vpeak
        
        for (j in 2:nrow(steps)) {
            target_vpeak <- quants$quants[i]
            # check for convergence to be explicit
            if (steps$time_converged[j] == 1 | steps$time_adv_converged[j] == 1) {
                if (steps$deporate[j] < quants$quants[i]) {
                    quants$vegtime[i] <- steps$time_request[j]
                    quants$brink_dist[i] <- steps$brink_dist[j]
                    break
                }
            }
        }
    }
    plot_quants (quants, steps)
}

plot_quants <- function (quants, steps) {
    # function to plot quantiles
    
    # quants = quantiles of vpeak
    # steps = stepped outputs
    
    # plot advance
    png ('advance.png', width = 1000, height = 500)
    plot (steps$time, steps$brink_dist, col = 'blue', type = 'l',
        ylab = 'brink distance (m)', xlab = 'time (a)')
    dev.off()

    # plot deporate
    png ('deporate.png', width = 1000, height = 500)
    plot (steps$time, steps$deporate, col = 'orange', type = 'l',
        ylab = 'deporate (m/a)', xlab = 'time (a)')    
    dev.off()
    
    # plot the quantiles of time till vegetation
    png ('q_advance.png', width = 1000, height = 500)
    plot (quants$probs, quants$vegtime, xlab = 'quantiles of vpeak (m/a)',
        pch = 19, col = 'black', ylab = 'time (a)')
    abline (v = 0.5, col = 'black')
    abline (v = c(0.25, 0.75), col = 'red')
    abline (v = c(0.1, 0.9), col = 'grey')
    dev.off()
    
    # plot the quantiles of brink dist
    png ('q_brink_dist.png', width = 1000, height = 500)
    plot (quants$probs, quants$brink_dist, xlab = 'quantiles of vpeak (m/a)',
        pch = 19, col = 'grey', ylab = 'brink_dist (m)')
    abline (v = 0.5, col = 'black')
    abline (v = c(0.25, 0.75), col = 'red')
    abline (v = c(0.1, 0.9), col = 'grey')
    dev.off()
    
    postscript ('q_brink_dist.eps')
    plot (quants$probs, quants$brink_dist, xlab = 'quantiles of vpeak (m/a)',
        pch = 19, col = 'grey', ylab = 'brink_dist (m)')
    abline (v = 0.5, col = 'black')
    abline (v = c(0.25, 0.75), col = 'red')
    abline (v = c(0.1, 0.9), col = 'grey')
    dev.off()
    
    png ('q_brink_dist_new.png', width = 1200, height = 1200, res = 216)
    plot (quants$probs, quants$brink_dist, xlab = 'Quantiles, Cumulative probability',
        type = 'l', col = 'red', ylab = 'Simulated max. advance (m)')
    abline (v = 0.5, col = 'black')
    abline (v = c(0.25, 0.75), col = 'blue')
    abline (v = c(0.1, 0.9), col = 'grey')
    dev.off()
}




