# Reactivation Advance Dimensioning System - RADS
# Tom Barchyn

# rads_environment.py - functions for creating fake environment parameters from distributions, etc.

def draw_flux ():
    """
    Function to draw a flux from a predefined distribution, parameters defined in rads_run.py
    or calling program. The parameters needed are:
    
    flux_meanlog = 3.3120468            # mean of lognormal distribution
    flux_sdlog = 0.2689315              # sd of lognormal distribution
    flux_minbound = 15.091557           # minimum bound
    flux_maxbound = 48.26884            # maximum bound
    """
    
    # make a random draw from the lognormal distribution
    flux = np.random.lognormal (flux_meanlog, flux_sdlog, 1)
    if flux > flux_maxbound:
        flux = flux_maxbound
    if flux < flux_minbound:
        flux = flux_minbound
    return (flux)
    
    