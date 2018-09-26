# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 13:35:41 2018

@author: rstreet
"""
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon
import copy

def calculate_lsst_field_visibility(fieldRA,fieldDec,start_date,end_date,
                                    min_alt=30.0,dt_days=1.0,
                                    diagnostics=False):
        """Method to calculate the visibility of a given RA and Dec from LSST 
        over the course of a year
        
        Adapted from an example in the Astropy docs.
        
        Inputs:
            :param float fieldRA: Field RA in decimal degrees
            :param float fieldDec: Field Dec in decimal degrees
        """
        
        field = SkyCoord(fieldRA, fieldDec, unit=(u.hourangle, u.deg))
        
        lsst = EarthLocation(lat=-30.239933333333333*u.deg, 
                             lon=-70.7429638888889*u.deg, 
                             height=2663.0*u.m)
        
        total_time_visible = 0.0
        t_start = Time(start_date+' 00:00:00')
        t_end = Time(end_date+' 00:00:00')
        cadence = 0.0007   # In days
        
        n_days = int((t_end - t_start).value)
        
        dates = np.array([t_start + \
                TimeDelta(i,format='jd',scale=None) for i in range(0,n_days,1)])
        
        target_alts = []
        hrs_visible_per_night = []
        hrs_per_night = []
        
        for d in dates:
            
            t  = copy.copy(d)
            t.out_subfmt = 'date'
            tstr = t.value
            
            intervals = np.arange(0.0,1.0,cadence)
            
            dt = TimeDelta(intervals,format='jd',scale=None)
            
            ts = d + dt
    
            frame = AltAz(obstime=ts, location=lsst)
            
            altaz = field.transform_to(frame)
            
            alts = np.array((altaz.alt*u.deg).value)
            
            idx = np.where(alts > min_alt)[0]
            
            sun_altaz = get_sun(ts).transform_to(frame)
            
            sun_alts = np.array((sun_altaz.alt*u.deg).value)
            
            jdx = np.where(sun_alts < 12.0)[0]
            
            hrs_per_night.append( cadence*len(sun_alts[jdx])*24.0 )
            
            idx = list(set(idx).intersection(set(jdx)))
            
            if len(idx) > 0:
                
                ts_vis = ts[idx]
                
                tvis = cadence * len(ts_vis)
                
                total_time_visible += tvis
                
                target_alts.append(alts[idx].max())
                
                print('Target visible from LSST for '+str(round(tvis*24.0,2))+\
                        'hrs on '+tstr)
                        
                hrs_visible_per_night.append((tvis*24.0))
                
            else:
                
                target_alts.append(-1e5)
                
                hrs_visible_per_night.append(0.0)
                
                print('Target not visible from LSST on '+tstr)
            
        if diagnostics:
            
            fig = plt.figure(1)

            plt.plot(ts.jd, alts, 'k-')

            plt.plot(ts.jd, sun_alts, 'y-')
            
            plt.plot(ts.jd, [30.0]*len(ts), 'r-.')
            
            plt.fill_between(ts.jd, 0, 90,
                 sun_altaz.alt > -18*u.deg, color='y', zorder=0)

            plt.gcf().autofmt_xdate()
            
            (xmin,xmax,ymin,ymax) = plt.axis()
            #plt.axis([ts.jd[1500],ts.jd[2000],0.0,90.0])
            
            plt.xlabel('JD')

            plt.ylabel('Field altitude above horizon [deg]')

            plt.savefig('target_visibility.png')
        
        return total_time_visible, hrs_visible_per_night

if __name__ == '__main__':
    
    if len(argv) > 1:
        
        fieldRA = argv[1]
        fieldDec = argv[2]
        start_date = argv[3]
        end_date = argv[4]
        
    else:
        
        fieldRA = input('Please enter the RA in sexigesimal format, J2000.0: ')
        fieldDec = input('Please enter the Dec in sexigesimal format, J2000.0: ')
        start_date = input('Please enter the start date of the observing window, YYYY-MM-DD: ')
        end_date = input('Please enter the end date of the observing window, YYYY-MM-DD: ')
        
    (total_time_visible,hrs_per_night) = calculate_lsst_field_visibility(fieldRA,fieldDec,
                                                         start_date,end_date,
                                                         diagnostics=True)
    