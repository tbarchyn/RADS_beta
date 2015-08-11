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


# rads_run.py - test parameters for the RADS system

import os
execfile ('rads.py')

##########################################################
# PARAMETERS
input_topo = 'dtm.tif'              # topography raster
input_bsmt = 'bsmt.tif'             # basement raster
point_file = 'eval_points.csv'      # points to evaluate (id, x_loc, y_loc)

# set environmental parameters
wind_dir = 301.0                    # formative wind direction
flux = 22.0                         # volumetric blowout cross-brink flux
vpeak = 0.9                         # vpeak in m / a
gamma = 15.0                        # characteristic gamma in degrees (must be positive)
n_processes = -1                    # number of processes to run in parallel

# step outputs
max_years = 50                      # set max number of years for stepping

# set flux montecarlo parameters
flux_meanlog = 3.3120468            # mean of lognormal distribution
flux_sdlog = 0.2689315              # sd of lognormal distribution
flux_minbound = 15.091557           # minimum bound
flux_maxbound = 48.26884            # maximum bound
montecarlo_iters = 50               # number of iterations to run

# debug parameters
debug_time = False                  # debug time convergence
debug_vpe = False                   # debug vpe convergence
debug_mbe = False                   # extra mbe convergence

##########################################################
# INSTRUCTIONS:
# probe 12: probe point with id 12 - probing fits to vpeak and creates plot and prints status
# run: run the full file in parallel
# step 12: step forward in time at point with id 12
# step_montecarlo 12: step forward in time with the montecarlo parameter pulls
# run_montecarlo 12: run the specified montecarlo iters for a given site
##########################################################
# EXECUTE!
if __name__ == '__main__':
    rads_execute ()




