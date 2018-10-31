# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 15:38:39 2018

@author: rstreet
"""

from astropy.time import Time
import matplotlib.pyplot as plt

def simulate_bd_lightcurve():
    """Function to simulate LSST DDF observations of a variable brown dwarf"""
    
    params = get_args()
    
    lci = read_luhman_lc(params['ip_lc_file'])
    lcz = read_luhman_lc(params['ip_lc_file'])
    
    vis_windows = parse_visibility_windows(params)
    
    down_lci = downsample_lightcurve(params,lci,vis_windows)
    down_lcz = downsample_lightcurve(params,lcz,vis_windows)
    
    plot_lc(lc1,lc2,params)
    
def get_args():
    
    params = { }

    if len(argv) > 1:
        params['ip_lc_file'] = argv[1]
        params['zp_lc_file'] = argv[2]
        params['visibility_file'] = argv[3]
        params['year'] = argv[4]
        params['cadence'] = float(argv[5])
        params['plot_file'] = argv[6]
            
    else:
        params['ip_lc_file'] = input('Please enter the path to the SDSS-i lightcurve: ')
        params['zp_lc_file'] = input('Please enter the path to the PS-z lightcurve: ')
        params['visibility_file'] = input('Please enter the path to the file of visibility windows: ')
        params['year'] = input('Please enter the year of the observations: ')
        params['cadence'] = float(input('Please enter the cadence of the observations in mins: '))
        params['plot_file'] = input('Please enter the output path for the plot: ')
        
    params['cadence'] = params['cadence']/(60.0*24.0) # mins -> days
    
    return params

def read_luhman_lc(file_path):
    """Function to read in the data for a BD lightcurve"""
    
    if path.isfile(file_path) == False:
        raise IOError('Cannot find '+file_path)
        exit()
    
    lines = open(file_path,'r').readlines()
    
    data = []
    
    for l in lines:
        
        if l[0:1] != '#':
            
            entries = l.replace('\n','').split()
            
            data.append( [float(entries[1]), float(entries[2], float(entries[3])] )
            
    data = np.array(data)
    
    return data
    
def parse_visibility_windows(params):
    """Function to read in the visibility windows, transposing the dates so 
    that they overlap the year of the observations, for the purposes of
    comparison."""
    
    if path.isfile(params['visibility_file']) == False:
        raise IOError('Cannot find '+params['visibility_file'])
        exit()
    
    lines = open(file_path,'r').readlines()
    
    data = []

    for l in lines:
        
        if l[0:1] != '#':
            
            (date, start_time, date, end_time) = l.replace('\n','').split()
            
            start_date = Time(params['year']+date[4:]+'T'+start_time, 
                                  format='isot', scale='utc')
            
            end_date = Time(params['year']+date[4:]+'T'+end_time, 
                                  format='isot', scale='utc')
                                  
            data.append( [start_date.jd, end_date.jd] )
            
    return data

def downsample_lightcurve(params,lc,vis_windows):
    """Function to downsample a lightcurve, by selecting observations only
    within the given visibility windows, and at the cadence given"""
    
    obs_start = lc[:,0].min()
    obs_end = lc[:,0].max()
    
    tol = 15.0/(60.0*24.0)      # Tolerance on time matching, days
    
    down_lc = []
    
    for (start_window, end_window) in vis_windows:
        
        if start_window >= obs_start and end_window <= obs_end:
            
            intervals = np.arange(start_window,end_window,params['cadence'])
            
            for ts in intervals:
                
                idx = np.where(lc[:,0] >= (ts - tol))[0]
                jdx = np.where(lc[:,0] <= (ts + tol))[0]
                kdx = list(set(idx).intersection(set(jdx)))
                
                down_lc.append( lc[idx[0],:] )
    
    down_lc = np.array(down_lc)
    
    return down_lc

def plot_lc(lc1,lc2,params):
    """Function to output the downsampled lightcurves"""
    
    fig = plt.figure(1,(10,10))
    
    plt.errorbar(lc1[:,0], lc1[:,1], yerr=lc1[:,2])
    
    plt.errorbar(lc2[:,0], lc2[:,1], yerr=lc2[:,2])
    
    plt.xlabel('HJD')
    plt.ylabel('Instrumental mag')
    
    plt.savefig(params['plot_file'])
    
if __name__ == '__main__':
    
    simulate_bd_lightcurve()
    