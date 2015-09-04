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


# rads_profiler.py: interpolation engine to produce profiles

class rads_profiler:
    """
    Class to contain and extent rads profiles for use in depth modeling
    """
    def __init__ (self, dtm, wind_dir):
        """
        Constructor initializes an interpolator that manages profile interpolation
        
        dtm = the dtm ref_raster object
        wind_dir = the compass azimuth associated with resultant wind direction
        """
        
        # flip the dtm
        dtm_flip = np.flipud (dtm.ras)          # flip the raster
        y_index_flip = dtm.y_index[::-1]        # flip the y index
        
        # create interpolator
        self.interpolator = RegularGridInterpolator (points = (y_index_flip, dtm.x_index), values = dtm_flip,
                                                    method = 'linear', bounds_error = True, fill_value = None)
        
        # set internal jumps
        self.wind_dir_rad = wind_dir * (math.pi / 180.0)
        
        self.inc = 0.5                                  # inc point spacing in meters
        self.step = 20.0                                # step spacing
        return
        
    def profile (self, x_loc, y_loc, steps, direction):
        """
        Method to output a profile from the location and return the distances
        along the profile and corresponding elevations. This goes from 0.0 distance
        to steps*self.step - self.inc as it uses
        np.arange to make arrays for interpolation. Note that the dists are negative
        if upwind. The 'zero point' is not included in upwind profiles either.
        
        x_loc = x location of profile start
        y_loc = y location of profile start
        steps = the number of step_spacings to go away from point
        direction = 'upwind' or 'downwind'
        """
        
        # compute distance arrays
        if direction == 'upwind':
            dists = np.arange (self.inc, (steps * self.step), self.inc)
            dists = -1.0 * dists        # make negative (upwind)
            dists = dists[::-1]         # sort from furthest upwind to closest
        else:
            dists = np.arange (0.0, (steps * self.step), self.inc)
        
        # compute point locations
        x_locs = x_loc - (dists * math.sin (self.wind_dir_rad))
        y_locs = y_loc - (dists * math.cos (self.wind_dir_rad))

        # interpolate the points, the interpolator should push error if out of raster bounds
        points = np.array (zip (y_locs, x_locs))         # note this is reversed as referencing is row, col
        
        try:
            elevs = self.interpolator (points)
            if len (elevs[np.isnan(elevs)]) > 0:
                elevs = x_locs * np.nan                  # if any values are missing, return all missing values
            
        except:
            elevs = x_locs * np.nan  
        
        return (dists, elevs)   
        
    def get_elev (self, x_loc, y_loc):
        """
        Function to interpolate the elevation at the point given
        
        x_loc = x location of profile start
        y_loc = y location of profile start
        """
        
        point = np.array([y_loc, x_loc])
        elev = self.interpolator (point)
        return (elev[0])




