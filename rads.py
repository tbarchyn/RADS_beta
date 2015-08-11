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


# rads.py: this file is the lead import for operating rads

import sys
import os
import shutil
import math
import numpy as np
import osr
import gdal
from joblib import Parallel, delayed
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize, brent

try:
    execfile ('../python/gdal_raster_io.py')
except:
    execfile ('./gdal_raster_io.py')
    
execfile ('./rads_profiler.py')
execfile ('./rads_execute.py')
execfile ('./rads_blowout.py')
execfile ('./rads_environment.py')



