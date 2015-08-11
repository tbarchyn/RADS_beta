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


# rads_blowout.py: this contains the blowout comparison engine

class blowout:
    """
    Generic blowout class that contains dimension minimization engines and related dimensional functions
    """
    def __init__ (self, id, x_loc, y_loc, dtm_pro, bsmt_pro, flux, vpeak, gamma):
        """
        Constructor populates some basic internals and gets initial profiles up and downwind
        
        id = a local id to store and report back
        x_loc = the x_location
        y_loc = the y location
        dtm_pro = the dtm rads_profiler object
        bsmt_pro = the bsmt rads_profiler object
        flux = the sediment flux
        vpeak = maximum vegetation deposition tolerance
        gamma = characteristic gamma associated with blowout slipfaces
        """

        # populate internal variables
        self.id = id
        self.x_loc = x_loc
        self.y_loc = y_loc
        self.dtm_pro = dtm_pro
        self.bsmt_pro = bsmt_pro
        self.flux = flux
        self.vpeak = vpeak
        self.gamma = gamma
        
        # internal geometry parameters
        self.dig = np.nan
        self.slipface_start_coord = np.nan
        self.uw_edge_coord = np.nan
        self.dw_edge_coord = np.nan
        self.target_coord = np.nan
        self.brink_dist = np.nan
        self.erosion = np.nan
        self.deposition = np.nan
        self.slip_bottom_elev = np.nan
        self.time = np.nan
        self.time_request = np.nan
        self.deporate = np.nan
        self.slip_area = np.nan
        self.advance = 0.0                  # NOTE: upon initialization, advance = 0
        
        # convergence error flags and results
        self.vpe_iter = 0                   # set vpe convergence iterations
        self.mbe_iter = 0                   # set mpe convergence iterations
        self.adv_iter = 0                   # set adv convergence iterations
        self.time_iter = 0                  # set time convergence iterations
        self.time_adv_iter = 0              # set time_adv convergence iterations
        self.vpe_converged = False          # set success flag after successful convergence
        self.mbe_converged = False          # set success flag after successful convergence
        self.adv_converged = False          # set success flag after successful convergence
        self.time_converged = False         # set success flag after successful convergence
        self.time_adv_converged = False     # set success flag after successful convergence
        
        # error codes
        self.no_surface_breakout = False    # vol_calc error: lack of surface breakout flag
        self.cant_extend_dw = False         # vol_calc_error: can't get enough dtm to calc properly
        self.cant_extend_uw = False         # vol_calc_error: can't get enough dtm to calc properly
        self.basement_hit = False           # basement hit
        self.status = True                  # fit status    

        # get surface and basement elevations at target cell
        self.surface_elev = self.dtm_pro.get_elev (x_loc = x_loc, y_loc = y_loc)
        self.bsmt_elev = self.bsmt_pro.get_elev (x_loc = x_loc, y_loc = y_loc)
    
        # get initial profiles
        self.dists = np.nan          # the distance array (negative is upwind)
        self.dtm = np.nan            # the dtm surface at corresponding distances
        self.bsmt = np.nan           # the bsmt surface at corresponding distances
        
        self.dists_uw = np.array([np.nan])
        self.dists_dw = np.array([np.nan])
        self.dtm_uw = np.array([np.nan])
        self.dtm_dw = np.array([np.nan])
        self.bsmt_uw = np.array([np.nan])
        self.bsmt_dw = np.array([np.nan])

        self.steps_uw = 0               # set initial step counter
        self.steps_dw = 0               # set initial step counter
        uw_status = self.extend_profile ('upwind')
        dw_status = self.extend_profile ('downwind')


        if not uw_status or not dw_status:
            self.status = False         # this site is not good, we can't get a prelim
                                        # profile
        
        return
        
    def extend_profile (self, direction):
        """
        Internal method to revise the internal profiles. this tries to extend, but
        if it goes off the raster, it returns False, and does not update internal
        profile variables.
        
        direction = 'upwind' or 'downwind'
        """
        
        if direction == 'upwind':
            _steps_uw = self.steps_uw + 1
            _dists_uw, _dtm_uw = self.dtm_pro.profile (self.x_loc, self.y_loc, _steps_uw, 'upwind')
            _dists_uw, _bsmt_uw = self.bsmt_pro.profile (self.x_loc, self.y_loc, _steps_uw, 'upwind')
            
            # check for any na's in the profile
            if len(_dtm_uw[np.isnan(_dtm_uw)]) > 0:
                success = False
            else:
                success = True
                self.steps_uw = _steps_uw
                self.dists_uw = _dists_uw
                self.dtm_uw = _dtm_uw
                self.bsmt_uw = _bsmt_uw
        
        else:
            _steps_dw = self.steps_dw + 1
            _dists_dw, _dtm_dw = self.dtm_pro.profile (self.x_loc, self.y_loc, _steps_dw, 'downwind')
            _dists_dw, _bsmt_dw = self.bsmt_pro.profile (self.x_loc, self.y_loc, _steps_dw, 'downwind')
        
            # check for any na's in the profile
            if len(_dtm_dw[np.isnan(_dtm_dw)]) > 0:
                success = False
            else:
                success = True
                self.steps_dw = _steps_dw
                self.dists_dw = _dists_dw
                self.dtm_dw = _dtm_dw
                self.bsmt_dw = _bsmt_dw        
        
        if success:
            # reconstruct profile with results and internals
            self.dists = np.append (self.dists_uw, self.dists_dw)
            self.dtm = np.append (self.dtm_uw, self.dtm_dw)
            self.bsmt = np.append (self.bsmt_uw, self.bsmt_dw)
        
        return (success)
    
    def bl_profile (self, dig, brink_dist, advance):
        """
        Internal method to create blowout profile along the present array of dists. This
        uses parabolic esque curves to create a loose approximation of the blowout.
        
        dig = the assigned test dig
        brink_dist = the distance from the focal cell to the brink
        advance = the downwind advance of the bottom of the slipface
        """
        
        self.bl = self.dists * np.nan
        
        # in the case of no advance, set stoss and upwind bottoms
        if advance == 0.0:
            stoss_bottom_elev = self.surface_elev - dig
            upwind_bottom_elev = self.surface_elev - dig
        # else, in the case of advance, set stoss and upwind bottoms with interpolation
        else:
            stoss_bottom_elev = np.interp (advance, xp = self.dists, fp = self.bsmt)
            upwind_bottom_elev = self.bsmt_elev

        # upwind profile: from dist = 0 upwind
        subs = self.dists <= 0.0
        upwind_power = 2.0
        upwind_slope = 0.005
        self.bl[subs] = upwind_bottom_elev + (upwind_slope * (-1.0 * self.dists[subs])**upwind_power)

        # middle profile: from dist = 0 to dist < advance
        subs = (self.dists > 0.0) & (self.dists < advance)
        self.bl[subs] = self.bsmt[subs]
        
        # stoss profile: from dists >= advance to dists < brink_dist
        subs = (self.dists >= advance) & (self.dists < brink_dist)
        stoss_power = 1.3
        stoss_slope = 0.039
        self.bl[subs] = stoss_bottom_elev + (stoss_slope * ((self.dists[subs] - advance)**stoss_power))

        # slipface profile
        slipface_jump = -0.6745085168424266
        self.brink_elev = stoss_bottom_elev + (stoss_slope * ((brink_dist - advance)**stoss_power))
        subs = self.dists >= brink_dist
        self.bl[subs] = self.brink_elev + ((self.dists[subs] - brink_dist) * slipface_jump)        

        return
    
    def calc_vols (self, dig, brink_dist, advance):
        """
        Method to create a blowout profile, ensure it surfaces, and return erosion and
        deposition volumes. This method runs a number of exception checks to evaluate
        
        brink_dist = the distance from the focal cell to the brink
        dig = the assigned test dig
        advance = the assigned test advance
        """
        
        # reset error flags and assigned internal variables
        self.no_surface_breakout = False
        self.cant_extend_dw = False
        self.cant_extend_uw = False
        self.erosion = numpy.nan
        self.deposition = numpy.nan
        self.time = numpy.nan
        
        self.target_coord = numpy.nan
        self.uw_edge_coord = numpy.nan
        self.dw_edge_coord = numpy.nan
        self.slipface_start_coord = numpy.nan
        
        # 1) check for an upwind surface breakout
        while True:
            self.bl_profile (dig, brink_dist, advance)
            diff = self.bl - self.dtm
            if len(self.bl[(self.dists < 0.0) & (diff > 0.0)]) > 0:
                break                           # we surfaced!
            else:
                uw_extend_success = self.extend_profile ('upwind')
                if not uw_extend_success:
                    self.cant_extend_uw = True
                    return                      # we can't get a profile upwind
        
        # 2) check that profile extends past our brink_dist
        while brink_dist > np.max (self.dists):
            # if the brink distance is greater than our dists, we for sure need to extend
            dw_extend_success = self.extend_profile ('downwind')
            if not dw_extend_success:
                self.cant_extend_dw = True
                return                          # we can't get a profile downwind
            
        self.bl_profile (dig, brink_dist, advance)      # recalculate the profile
        diff = self.bl - self.dtm                       # recalculate the diff

        # 3) check for a downwind surface breakout
        if len(diff[(self.dists > 0.0) & (diff > 0.0)]) == 0:
            self.no_surface_breakout = True
            return                              # no surface breakout!
        
        # 4) check to make sure the profile ends up below the topo (end of slipface)
        while True:
            if (diff[len(diff) - 1] < 0.0):
                break       # the last point in the diff is below, passed check!
            else:
                success = self.extend_profile ('downwind')
                if not success:
                    self.cant_extend_dw = True
                    return
                
                self.bl_profile (dig, brink_dist, advance)  # recalculate the profile
                diff = self.bl - self.dtm                   # recalculate the diff
                    
        # 5) get the target coordinate
        target_coord_return = np.where ((self.dists < 0.00001) & (self.dists > -0.00001))
        self.target_coord = int(target_coord_return[0])     # this obtuse coercion is unfortunately pythonic
        
        # 6) determine the upwind edge with an upwind walk
        for i in range(self.target_coord, -1, -1):
            if diff[i] < 0.0 and diff[i-1] > 0.0:
                self.uw_edge_coord = i         # this is the last point in the blowout
                break
        
        # catch the exception where we can't find the upwind edge
        if np.isnan(self.uw_edge_coord):
            if debug_vpe:
                self.print_status()
            return                          # we can't do much with this site
        
        # 7) determine the start of the slipface
        try:
            self.slipface_start_coord = np.min (np.where (self.dists > brink_dist))
        except:
            # ERROR occurs here when there is no slipface present!
            # just return out of here as we can't do much with this brink_dist value
            return
        
        # 8) now, get the end of the slipface
        for i in range(self.slipface_start_coord, len(self.dists)):
            if diff[i-1] > 0.0 and diff[i] < 0.0:
                self.dw_edge_coord = i
                break
         
        # catch the exception where we can't find the slipface edge
        if np.isnan(self.dw_edge_coord):
            if debug_vpe:
                self.print_status()
            return                          # we can't do much with this site
         
        # 9) calculate erosion and deposition, and time it would take to perform excavation
        blowout_diffs = diff[self.uw_edge_coord:self.dw_edge_coord]             # subset diffs
        self.erosion = -1.0 * np.sum (blowout_diffs[blowout_diffs < 0.0])       # calc erosion sum
        self.deposition = np.sum (blowout_diffs[blowout_diffs > 0.0])           # calc depo sum
        self.time = ((self.erosion + self.deposition) / 2.0) / self.flux        # calc time
        
        # debug outputs
        if debug_mbe:
            print ('brink_dist: ' + str(self.brink_dist))
            print ('erosion: ' + str(self.erosion))
            print ('deposition: ' + str(self.deposition))
            print ('slipface_start_coord: ' + str(self.slipface_start_coord))
            print ('uw_edge_coord: ' + str(self.uw_edge_coord))
            print ('dw_edge_coord: ' + str(self.dw_edge_coord))
            print ('target_coord: ' + str(self.target_coord))
            print ('-------------------------------------------------------')
                        
        return
        
    def fit_mbe (self, dig, advance):
        """
        Method to minimize the brink_dist such that the mass balance error
        is close to 0.0. This numerically tries different brink distances
        with a given dig until the mbe is as low as it can go. We cannot fit a brink
        dist that is less than the advance + 1.0, we need at least a 1 m long stoss
        slope (in reality this is unlikely to be a problem).
        
        dig = the assigned test dig
        advance = the assigned test advance
        """
        
        mbe_tolerance = 0.1
        min_bracket_start = advance + 1.0       # the brink_dist cannot be less than advance + 1.0!
        max_bracket_start = 500.0
        min_bracket = min_bracket_start
        max_bracket = max_bracket_start
        max_iterations = 50
        iter = 0
        guess = max_bracket
        
        # reset the assigned coordinates
        self.mbe_converged = False
        self.mbe_iter = numpy.nan
        self.brink_dist = numpy.nan
        self.mbe_iter = numpy.nan
        
        # check the dig, and ensure we have not dug into the basement
        if dig >= (self.surface_elev - self.bsmt_elev):
            # if we have dug into the basement, still converge, but
            # set the basement hit flag high
            dig = (self.surface_elev - self.bsmt_elev)
            self.basement_hit = True
        else:
            self.basement_hit = False
        
        if self.status:
            while True:
                iter = iter + 1
                if iter > max_iterations:
                    brink_dist = np.nan
                    break
                error = mbe_tolerance + 1.0                 # set high, then reassign
                old_guess = guess                           # record previous guess
                guess = (min_bracket + max_bracket) / 2.0   # guess halfway
                if abs(old_guess - guess) < 0.0001:
                    brink_dist = numpy.nan
                    break                                   # we are trying the same value, break
                
                self.calc_vols (dig, guess, advance)        # calc volumes
                
                if self.no_surface_breakout:
                    min_bracket = guess             # we are too low, up the min_bracket
                else:
                    if self.cant_extend_uw or self.cant_extend_dw:
                        max_bracket = guess         # we are too high, get lower
                    else:
                        mbe = self.erosion - self.deposition
                        error = abs(mbe)            # get difference
                        if error < mbe_tolerance:
                            brink_dist = guess
                            break                   # we have succeeded!
                        else:
                            if mbe > 0.0:
                                min_bracket = guess  # we are too low, get higher
                            else:
                                max_bracket = guess  # we are too high, get lower
        
            # check the values aren't too close to the bounds
            if not np.isnan (brink_dist) and (brink_dist < (min_bracket_start + 0.1)):
                brink_dist = np.nan            # too close to min bracket
            if not np.isnan (brink_dist) and (brink_dist > (max_bracket_start - 0.1)):
                brink_dist = np.nan            # too close to max bracket
        else:
            brink_dist = np.nan                # status = false
        if not np.isnan(brink_dist):
            self.mbe_converged = True               # we succeeded
        else:
            self.mbe_converged = False              # we did not succeed

        self.mbe_iter = iter            # set convergence iterations for optimal fit
        self.brink_dist = brink_dist    # record the value
        return
    
    def calc_deporate (self):
        """
        Internal method to calculate the deposition rate and populate internal variables
        """
        repose = 34.0                       # repose in degrees

        # use more accurate slipface sizing to ensure the minimization converges properly
        # this calculates the slipface bottom edge with linear approximation between the
        # relevant distance samples.
        dtm1 = self.dtm[self.dw_edge_coord - 1]                 # get the dtm elev before cross over
        dtm2 = self.dtm[self.dw_edge_coord]                     # get the dtm elev after cross over
        bl1 = self.bl[self.dw_edge_coord - 1]                   # get blowout profile before cross over
        bl2 = self.bl[self.dw_edge_coord]                       # get blowout profile after cross over
        x = 1.0 - ((dtm2 - bl2) / ((bl1-bl2)-(dtm1-dtm2)))      # get the cross over fraction along the sample
        self.slip_bottom_elev = (bl1 * (1.0-x)) + (bl2 * x)     # get weighted average
        slip_ht = self.brink_elev - self.slip_bottom_elev
        slip_wd = slip_ht / math.tan (repose * (math.pi / 180.0))
        
        # calculate the area of the slipface and characteristic deposition rate
        b_up = 1.0
        b_down = b_up + (2.0 * slip_wd * math.tan (self.gamma * math.pi / 180.0))
        if slip_wd == 0.0:
            self.deporate = np.nan
            return
        
        self.slip_area = ((b_up + b_down) / 2.0) * slip_wd
        
        if self.slip_area == 0.0:
            self.deporate = np.nan
            return
                
        self.deporate = self.flux / self.slip_area
        return
        
    def calc_vpe (self, dig, advance):
        """
        Internal method to calculate the deviation between vpeak and the predicted
        vpeak from a given blowout profile. Retuns nan if the fit mbe routine returns
        nan
        
        dig = the assigned test dig
        advance = the assigned test advance
        """

        self.fit_mbe (dig, advance)         # fit the profile
        
        if not np.isnan (self.brink_dist):
            self.calc_deporate()
            vpe = self.deporate - self.vpeak
            
            if debug_vpe:
                print ('--------------------------------------------------------------------')
                print ('VPE CONVERGENCE REPORT: brink_dist: ' + str(self.brink_dist) + ', advance: ' + str(advance))
                print ('--------------vpe: ' + str(vpe))
                print ('--------------slip_wd: ' + str(slip_wd))
                print ('--------------slip_ht: ' + str(slip_ht))
                print ('--------------slip_area: ' + str(self.slip_area))
                print ('--------------deporate: ' + str(self.deporate))
                print ('--------------slip_bottom_elev: ' + str(self.slip_bottom_elev))
            
            return (vpe)
        else:
            return (np.nan)     # code for not possible to fit profile
        
    def calc_time_error (self, dig, advance):
        """
        Internal method to calculate the time error from a given fit
                
        dig = the assigned test dig
        advance = the assigned test advance
        """
        self.fit_mbe (dig, advance)         # fit the profile
         
        if np.isnan(self.brink_dist):
            time_error =  np.nan
        else:
            time_error = self.time_request - self.time 
            
        if debug_time:
            print ('--------------------------------------------------------------------')
            print ('TIME CONVERGENCE REPORT: brink_dist: ' + str(self.brink_dist) + ', advance: ' + str(advance))
            print ('--------------time_error: ' + str(time_error))
            print ('--------------request dig: ' + str(dig))
            print ('--------------request advance: ' + str(advance))

        return (time_error)
        
    def fit_evolution (self, variable):
        """
        Method to fit various digs until we get one with values corresponding to
        a vpeak, this will correspond with the maximum advance at a site.
        
        If the converged dig looks like it is deeper than the basement, we call
        the fit_advance function to continue with the fit under various advances
        
        variable = the variable to fit
        """
        
        vpe_tolerance = 0.05                    # in m/a
        time_tolerance = 0.05                   # in years
        min_bracket_start = 0.5
        max_bracket_start = 25.0
        min_bracket = min_bracket_start
        max_bracket = max_bracket_start
        max_iterations = 50
        iter = 0
        guess = max_bracket
        
        # assign tolerances
        if variable == 'vpe':
            tolerance = vpe_tolerance
            self.vpe_converged = False
            self.vpe_iter = numpy.nan
        elif variable == 'time':
            tolerance = time_tolerance
            self.time_converged = False
            self.time_iter = numpy.nan
        else:
            print ('ERROR: unrecognized fit_dig variable')
            sys.exit()

        # reset the assigned coordinates
        self.dig = numpy.nan
        self.advance = 0.0      # ensure start of vpe fit involves no advance
        advance = 0.0           # we start the vpe fit with digs, so advance is 0.0
        
        if self.status:
            while True:
                iter = iter + 1
                if iter > max_iterations:
                    opt_dig = np.nan
                    break
                error = tolerance + 1.0                     # set high, then reassign
                old_guess = guess
                guess = (min_bracket + max_bracket) / 2.0   # guess halfway
                if abs(old_guess - guess) < 0.0001:
                    opt_dig = numpy.nan
                    break                                   # we are trying the same value, break
                
                if variable == 'vpe':
                    r = self.calc_vpe (guess, advance)      # calc vpeak error
                if variable == 'time':
                    r = self.calc_time_error (guess, advance)# calc time error
                    
                if np.isnan(r):
                    max_bracket = guess                     # we are too deep
                else:
                    error = abs(r)
                    if error < tolerance:
                        opt_dig = guess                     # we have succeeded
                        break
                    else:
                        if r > 0.0:
                            if self.basement_hit:
                                opt_dig = guess             # we hit basement!
                                self.dig = opt_dig          # set the dig for the advance fit routine                
                                self.fit_advance (variable) # run the advance fit routine!
                                self.basement_hit = True    # ensure we exit with positive basement hit
                                break                       
                            else:
                                min_bracket = guess         # we are too shallow
                        else:
                            max_bracket = guess             # we are too deep

            # check those values that they are not too close to the brackets
            if not np.isnan(opt_dig) and (opt_dig < min_bracket_start + 0.1):
                opt_dig = np.nan            # too close to min bracket
            if not np.isnan(opt_dig) and (opt_dig > max_bracket_start - 0.1):
                opt_dig = np.nan            # too close to max bracket
        else:
            opt_dig = np.nan                # status == false
        
        # set convergence flags and record convergence iterations
        if variable == 'vpe':
            if not np.isnan(opt_dig):
                self.vpe_converged = True       # we succeeded!
            else:   
                self.vpe_converged = False      # we did not succeed
            self.vpe_iter = iter                # set convergence iterations
        
        if variable == 'time':
            if not np.isnan(opt_dig):
                self.time_converged = True      # we succeeded!
            else:   
                self.time_converged = False     # we did not succeed
            self.time_iter = iter               # set convergence iterations
        
        self.dig = opt_dig                      # record the optimum dig
        return
        
    def fit_advance (self, variable):
        """
        Internal method to minimize the extension to match a given vpeak, this continues from
        the dig minimization routine to move the dune out, and away from the target point. We
        hold the dig constant and try progressive advances until we find the best advance.
                
        variable = the variable to fit
        """
        
        adv_tolerance = 0.05
        time_adv_tolerance = 0.05
        min_bracket_start = 0.0
        max_bracket_start = 500.0
        min_bracket = min_bracket_start
        max_bracket = max_bracket_start
        max_iterations = 50
        iter = 0
        guess = max_bracket
        
        if variable == 'vpe':
            tolerance = adv_tolerance
            self.adv_converged = False
            self.adv_iter = numpy.nan
        if variable == 'time':
            tolerance = time_adv_tolerance
            self.time_adv_converged = False
            self.time_adv_iter = numpy.nan
        
        # reset the assigned coordinates
        self.advance = numpy.nan        # we set internal advance to be unassigned
        
        if self.status:
            while True:
                iter = iter + 1
                if iter > max_iterations:
                    advance = np.nan
                    break
                error = tolerance + 1.0                 # set high, then reassign
                old_guess = guess
                guess = (min_bracket + max_bracket) / 2.0   # guess halfway
                if abs(old_guess - guess) < 0.0001:
                    advance = numpy.nan
                    break                                   # we are trying the same value, break
                
                if variable == 'vpe':
                    r = self.calc_vpe (self.dig, guess)       # calc vpeak error with previously assigned dig
                if variable == 'time':
                    r = self.calc_time_error (self.dig, guess)# calc time error
                    
                if np.isnan(r):
                    max_bracket = guess                     # we are too far downwind probably
                else:
                    error = abs(r)
                    if error < adv_tolerance:
                        advance = guess                     # we have succeeded
                        break
                    else:
                        if r > 0.0:
                            min_bracket = guess             # advance is too far upwind
                        else:
                            max_bracket = guess             # advance is too far downwind

            # check those values that they are not too close to the brackets
            if not np.isnan(advance) and (advance < min_bracket_start + 0.1):
                advance = np.nan            # too close to min bracket
            if not np.isnan(advance) and (advance > max_bracket_start - 0.1):
                advance = np.nan            # too close to max bracket
        else:
            advance = np.nan                # status == false
        
        # set the convergence and iteration flags
        if variable == 'vpe':
            if not np.isnan(advance):
                self.adv_converged = True       # we succeeded!
            else:   
                self.adv_converged = False      # we did not succeed
            self.adv_iter = iter                # set convergence iterations
        if variable == 'time':
            if not np.isnan(advance):
                self.time_adv_converged = True  # we succeeded!    
            else:
                self.time_adv_converged = False # we did not succeed
            self.time_adv_iter = iter           # set convergence iterations
        
        self.advance = advance              # record the opt_adv
        return
    
    def write_profile (self, output_filename):
        """
        Method to write out the profile to an external file for analysis
        
        output_filename = the output filename
        """
        
        num_parameters = 5
        headers = 'dists,dtm,bsmt,bl,id'
        oput = np.zeros ((len(self.dists), num_parameters))
        oput = oput * np.nan
        
        # populate the rows
        oput[:,0] = self.dists
        oput[:,1] = self.dtm
        oput[:,2] = self.bsmt
        oput[:,3] = self.bl
        oput[:,4] = self.id
        
        np.savetxt (output_filename, oput, delimiter = ',', header = headers, comments = '')
        
        return
        
    def return_status (self):
        """
        Method to return a tuple of results to store in an output table
        """
        return (self.id,self.x_loc, self.y_loc, self.surface_elev, self.bsmt_elev, self.dig,
                self.brink_dist, self.erosion, self.deposition, self.flux, self.vpeak,
                self.time, self.time_request,
                self.deporate, self.slip_area, self.advance, self.status, self.vpe_converged,
                self.mbe_converged, self.adv_converged, self.time_converged, self.time_adv_converged,
                self.vpe_iter, self.mbe_iter, self.adv_iter, self.time_iter, self.time_adv_iter,
                self.no_surface_breakout, self.cant_extend_uw, self.cant_extend_dw, self.basement_hit)
    
    def print_status (self):
        """
        Method to print the internal status to review the present status of the blowout
        """
        print ('STATUS REPORT')
        print ('========================================')
        print ('id: ' + str(self.id))
        print ('x_loc: ' + str(self.x_loc))
        print ('y_loc: ' + str(self.y_loc))
        print ('surface_elev: ' + str(self.surface_elev))
        print ('bsmt_elev: ' + str(self.bsmt_elev))
        print ('dig: ' + str(self.dig))
        print ('brink_dist: ' + str(self.brink_dist))
        print ('erosion: ' + str(self.erosion))
        print ('deposition: ' + str(self.deposition))
        print ('slip_bottom_elev: ' + str(self.slip_bottom_elev))
        print ('time: ' + str(self.time))
        print ('time_request: ' + str(self.time_request))
        print ('deporate: ' + str(self.deporate))
        print ('slip_area: ' + str(self.slip_area))
        print ('advance: ' + str(self.advance))
        print ('----------------------------------------')
        print ('slipface_start_coord: ' + str(self.slipface_start_coord))
        print ('uw_edge_coord: ' + str(self.uw_edge_coord))
        print ('dw_edge_coord: ' + str(self.dw_edge_coord))
        print ('target_coord: ' + str(self.target_coord))
        print ('----------------------------------------')
        print ('status: ' + str(self.status))
        print ('vpe_converged: ' + str(self.vpe_converged))
        print ('mbe_converged: ' + str(self.mbe_converged))
        print ('adv_converged: ' + str(self.adv_converged))
        print ('time_converged: ' + str(self.time_converged))
        print ('time_adv_converged: ' + str(self.time_adv_converged))
        print ('vpe_iterations: ' + str(self.vpe_iter))
        print ('mbe_iterations: ' + str(self.mbe_iter))
        print ('adv_iterations: ' + str(self.adv_iter))
        print ('time_iterations: ' + str(self.time_iter))
        print ('time_adv_iterations: ' + str(self.time_adv_iter))
        print ('no_surface_breakout: ' + str(self.no_surface_breakout))
        print ('cant_extend_uw: ' + str(self.cant_extend_uw))
        print ('cant_extend_dw: ' + str(self.cant_extend_dw))
        print ('basement_hit: ' + str(self.basement_hit))
        print ('steps_uw: ' + str(self.steps_uw))
        print ('steps_dw: ' + str(self.steps_dw))
        print ('========================================')
        return
    