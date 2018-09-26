# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 17:11:20 2018

@author: rstreet
"""

import numpy as np
import calculate_lsst_field_visibility_astropy

def calc_expected_visits(ra,dec,cadence,start_date,end_date):
    """Function to calculate the maximum possible number of visits to a 
    given pointing, given the expected cadence of observation and within
    the date ranges given, taking target visibility into account.
    
    Input:
    :param string ra:           RA, J2000.0, sexigesimal format
    :param string dec:          Dec, J2000.0, sexigesimal format
    :param float cadence:       Interval between successive visits in the 
                                same single filter in hours
    :param string start_date:   Start of observing window YYYY-MM-DD
    :param string start_date:   End of observation window YYYY-MM-DD
    
    Output:
    :param array n_visits:      Number of visits possible per night
    """
    
    (total_time_visible, hrs_visible_per_night) = calculate_lsst_field_visibility_astropy.calculate_lsst_field_visibility(ra,dec,start_date,end_date)
    
    n_visits = int(np.array(hrs_visible_per_night) / cadence)
    
    return n_visits