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


# rads_execute.py - main run function to run the dimensioning system

def rads_welcome ():
    """
    Function to print welcome message
    """
    print ('-----------------------------------------')
    print ('Welcome to RADS - the Reactivation Advance Dimensioning System')
    print ('Developed at the University of Calgary')
    print ('This program has no warranty whatsoever')
    print ('-----------------------------------------')
    return

def get_instructions ():
    """
    Function to get instructions from the command line
    """
    try:
        instruction = sys.argv[1]
        try:
            target_id = int (sys.argv[2])
            return (instruction, target_id)
        except:
            return (instruction, None)
    except:
        print ('ERROR: I need instructions!')
        sys.exit()
        
def rads_execute ():
    """
    Executes the calculation engine for the supplied points
    
    Check prototype rads_test.py parameter file for details. Note that this subtly depends on
    the scoping rules in python which allow external variable read access. To update in the
    future, perhaps pass variables as a dict, or list, or the path to the file to help make
    the program structure easier to refactor to another language.
    """
    
    rads_welcome()
    instruction, target_id = get_instructions()
    
    # read in the points from the point file
    points = np.loadtxt (point_file, skiprows = 1, delimiter = ',')
    if len(points.shape) == 1:
        points = np.array ([points])

    # read in the raster files
    dtm = ref_raster (input_topo)         # raw dtm
    bsmt = ref_raster (input_bsmt)        # raw basement
    
    oput_header = ('id,x,y,surface_elev,bsmt_elev,dig,brink_dist,erosion,deposition,flux,vpeak,time,time_request,deporate,slip_area,advance,status,vpe_converged,mbe_converged,adv_converged,time_converged,time_adv_converged,vpe_iter,mbe_iter,adv_iter,time_iter,time_adv_iter,no_surface_breakout,cant_extend_uw,cant_extend_dw,basement_hit')
    print ('Setup complete . . ')
    print ('-----------------------------------------')

    ###############################################################
    # PROBE A POINT
    if instruction == 'probe':
        point = points[points[:,0] == target_id,:]
        rads_probe_point (point[0,:], dtm, bsmt, wind_dir, vpeak, gamma)
    
    ###############################################################
    # RUN THE FULL POINT FILE
    elif instruction == 'run':
        # create a global dtm profiler object
        dtm_pro = rads_profiler (dtm, wind_dir)
        bsmt_pro = rads_profiler (bsmt, wind_dir)

        res = Parallel (n_jobs = n_processes, batch_size = 1, verbose = 10) (delayed (rads_process_point) 
                        (points[i,:], dtm_pro, bsmt_pro, wind_dir, vpeak, gamma) for i in range(0, points.shape[0]))
        oput = np.array (res)

        # uncomment for serial debugging
        #for i in range(0, points.shape[0]):
        #    res = rads_process_point (points[i,:], dtm_pro, bsmt_pro, wind_dir, vpeak, gamma)
        
        # output csv file
        np.savetxt ('rads_output.csv', oput, delimiter = ',', header = oput_header, comments = '')
        print ('Output file written . .')
    
    ###############################################################
    # STEP A POINT
    elif instruction == 'step':
        point = points[points[:,0] == target_id,:]
        oput = rads_step_point (point[0,:], dtm, bsmt, wind_dir, vpeak, gamma)

        # output csv file
        np.savetxt ('rads_step_output.csv', oput, delimiter = ',', header = oput_header, comments = '')
        print ('Output file written . .')
    
    ###############################################################
    # STEP A POINT MONTECARLO
    elif instruction == 'step_montecarlo':
        point = points[points[:,0] == target_id,:]
        oput = rads_step_point_montecarlo (point[0,:], dtm, bsmt, wind_dir, vpeak, gamma)
    
        # output csv file
        np.savetxt ('rads_step_montecarlo_output.csv', oput, delimiter = ',', header = oput_header, comments = '')
        print ('Output file written . .')
    
    ###############################################################
    # RUN A POINT MONTECARLO
    elif instruction == 'run_montecarlo':
        point = points[points[:,0] == target_id,:]
        
        # create a global dtm profiler object
        dtm_pro = rads_profiler (dtm, wind_dir)
        bsmt_pro = rads_profiler (bsmt, wind_dir)
        
        # process iterations in parallel
        res = Parallel (n_jobs = n_processes, batch_size = 1, verbose = 10) (delayed (rads_run_point_montecarlo) 
                        (points[0,:], dtm_pro, bsmt_pro, wind_dir, vpeak, gamma) for i in range(0, montecarlo_iters))
        
        oput = np.array (res)
        print (str(oput.shape))
        
        # create directory of results
        if os.path.isdir ('run_montecarlo_profiles'): 
            shutil.rmtree ('run_montecarlo_profiles')
            os.mkdir ('run_montecarlo_profiles')
        else:
            os.mkdir ('run_montecarlo_profiles')

        # output csv files
        for i in range(0, oput.shape[0]):
            filename = './/run_montecarlo_profiles//r_' + str(i) + '.csv'
            np.savetxt (filename, oput[i,:,:], delimiter = ',', header = oput_header, comments = '')
        
        print ('Output files written . .')
    
    else:
        print ('ERROR: undefined instruction')
        
    print ('-----------------------------------------')
    return

def rads_probe_point (point, dtm, bsmt, wind_dir, vpeak, gamma):
    """
    Function called to probe the fit of a given point
    
    point = a np array of length 3, with id, x, y positions
    dtm = the dtm object
    bsmt = the basement object
    wind_dir = the wind direction
    flux = the sediment flux
    vpeak = the vegetation deposition tolerance
    gamma = representative gamma value
    """
    
    dtm_pro = rads_profiler (dtm, wind_dir)
    bsmt_pro = rads_profiler (bsmt, wind_dir)

    # get id, x_loc, and y_loc
    id = point[0]
    x_loc = point[1]
    y_loc = point[2]
    
    # create the blowout, fit the vpe, print status and report profile
    bl = blowout (id, x_loc, y_loc, dtm_pro, bsmt_pro, flux, vpeak, gamma)
    bl.fit_evolution ('vpe')                    # fit the dig
    bl.print_status()                           # print results
    bl.write_profile ('profile.csv')            # output the results
    os.system ('Rscript.exe rads_plot.R')       # plot the results
    return
        
def rads_process_point (point, dtm_pro, bsmt_pro, wind_dir, vpeak, gamma):
    """
    Function called in parallel to process a given point
    
    point = a np array of length 3, with id, x, y positions
    dtm_pro = the dtm_pro object
    bsmt_pro = the basement_pro object
    wind_dir = the wind direction
    flux = the sediment flux
    vpeak = the vegetation deposition tolerance
    gamma = representative gamma value
    """

    # get id, x_loc, and y_loc
    id = point[0]
    x_loc = point[1]
    y_loc = point[2]
    
    bl = blowout (id, x_loc, y_loc, dtm_pro, bsmt_pro, flux, vpeak, gamma)
    bl.fit_evolution ('vpe')                    # fit the dig
    res = bl.return_status()                    # get the results
    return (res)

def rads_step_point (point, dtm, bsmt, wind_dir, vpeak, gamma):
    """
    Function to step a point forward in time with static parameters
    
    point = a np array of length 3, with id, x, y positions
    dtm = the dtm object
    bsmt = the basement object
    wind_dir = the wind direction
    flux = the sediment flux
    vpeak = the vegetation deposition tolerance
    gamma = representative gamma value
    """
    
    dtm_pro = rads_profiler (dtm, wind_dir)
    bsmt_pro = rads_profiler (bsmt, wind_dir)

    # get id, x_loc, and y_loc
    id = point[0]
    x_loc = point[1]
    y_loc = point[2]

    # create the blowout
    bl = blowout (id, x_loc, y_loc, dtm_pro, bsmt_pro, flux, vpeak, gamma)
    
    # set up profile directory
    if os.path.isdir ('step_profiles'): 
        shutil.rmtree ('step_profiles')
        os.mkdir ('step_profiles')
    else:
        os.mkdir ('step_profiles')
    
    # loop through years, forward in time
    for years in range(1, max_years + 1):
        bl.flux = flux                          # use constant flux
        bl.time_request = float(years)          # set the requested timefit
        bl.fit_evolution ('time')               # fit the blowout to time
        if not np.isnan (bl.brink_dist):
            bl.calc_deporate ()                 # populate internal parameters
        res = bl.return_status ()               # get a status
        if years == 1:
            oput = np.array([res])
        else:
            res = np.array([res])
            oput = np.vstack([oput, res])
        
        # output a profile if we had convergence
        if bl.time_converged or bl.time_adv_converged:
            profile_filename = os.path.join('./step_profiles', 't_' + str(years) + '.csv')
            bl.write_profile (profile_filename)
        
        if debug_time:
            bl.print_status()
            
        print ('completed iteration: ' + str(years))    
    
    # print status
    bl.print_status ()
    
    # run the R script to make plots 
    os.system ('Rscript.exe rads_step_plot.R step_profiles ' + str(id))
    
    return (oput)

def rads_step_point_montecarlo (point, dtm, bsmt, wind_dir, vpeak, gamma):
    """
    Function to step a point forward in time with montecarlo parameter
    pulls (where defined).
    
    point = a np array of length 3, with id, x, y positions
    dtm = the dtm object
    bsmt = the basement object
    wind_dir = the wind direction
    flux = the sediment flux
    vpeak = the vegetation deposition tolerance
    gamma = representative gamma value
    """
    
    dtm_pro = rads_profiler (dtm, wind_dir)
    bsmt_pro = rads_profiler (bsmt, wind_dir)

    # get id, x_loc, and y_loc
    id = point[0]
    x_loc = point[1]
    y_loc = point[2]

    # create the blowout
    bl = blowout (id, x_loc, y_loc, dtm_pro, bsmt_pro, flux, vpeak, gamma)
    
    # set up profile directory
    if os.path.isdir ('step_montecarlo_profiles'): 
        shutil.rmtree ('step_montecarlo_profiles')
        os.mkdir ('step_montecarlo_profiles')
    else:
        os.mkdir ('step_montecarlo_profiles')
    
    # loop through years, forward in time
    cum_ero = 0.0
    for years in range(1, max_years + 1):
        bl.flux = draw_flux()                   # get montecarlo flux
        cum_ero = cum_ero + bl.flux             # set the cumulative erosion
        bl.time_request = cum_ero / bl.flux
        bl.fit_evolution ('time')               # fit the blowout to time
        if not np.isnan (bl.brink_dist):
            bl.calc_deporate ()                 # populate internal parameters
        res = bl.return_status ()               # get a status
        if years == 1:
            oput = np.array([res])
        else:
            res = np.array([res])
            oput = np.vstack([oput, res])
        
        # output a profile if we had convergence
        if bl.time_converged or bl.time_adv_converged:
            profile_filename = os.path.join('./step_montecarlo_profiles', 't_' + str(years) + '.csv')
            bl.write_profile (profile_filename)
        
        if debug_time:
            bl.print_status()
            
        print ('completed iteration: ' + str(years))    
    
    # print status
    bl.print_status ()

    # run the R script to make plots 
    os.system ('Rscript.exe rads_step_plot.R step_montecarlo_profiles ' + str(id))
    
    return (oput)
    
def rads_run_point_montecarlo (point, dtm_pro, bsmt_pro, wind_dir, vpeak, gamma):
    """
    Function to step a point forward in time with montecarlo parameter
    pulls (where defined). This is executed in parallel repeatedly.
    
    point = a np array of length 3, with id, x, y positions
    dtm_pro = the dtm profiler object
    bsmt_pro = the basement profiler object
    wind_dir = the wind direction
    flux = the sediment flux
    vpeak = the vegetation deposition tolerance
    gamma = representative gamma value
    """

    # get id, x_loc, and y_loc
    id = point[0]
    x_loc = point[1]
    y_loc = point[2]

    # create the blowout
    bl = blowout (id, x_loc, y_loc, dtm_pro, bsmt_pro, flux, vpeak, gamma)
    
    # loop through years, forward in time
    cum_ero = 0.0                               # set cumulative erosion
    for years in range(1, max_years + 1):
        bl.flux = draw_flux()                   # get montecarlo flux
        cum_ero = cum_ero + bl.flux             # set the cumulative erosion
        bl.time_request = cum_ero / bl.flux
        bl.fit_evolution ('time')               # fit the blowout to time
        if not np.isnan (bl.brink_dist):
            bl.calc_deporate ()                 # populate internal parameters
        res = bl.return_status ()               # get a status
        if years == 1:
            oput = np.array([res])
        else:
            res = np.array([res])
            oput = np.vstack([oput, res])
        
        if debug_time:
            bl.print_status()
    
    return (oput)

    
    
    
    
    
    