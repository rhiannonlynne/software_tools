# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 21:38:32 2018

@author: rstreet
"""
import vizier_tools
import lsst_class
import wfirst_fov
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
from matplotlib.patches import Circle
from os import path
from sys import argv
import numpy as np

def optimize_survey_footprint():
    """Function to optimize the DDF survey footprint"""
    
    params = get_args()

    lsst = lsst_class.LSSTFootprint(ra_centre=params['ra'],
                                    dec_centre=params['dec'])
    
    wfirst = wfirst_fov.WFIRSTFootprint()
    
    params['fov'] = lsst.radius*2.0
    params['fov'] = 20.0
    params['lsst_pixscale'] = lsst.pixel_scale
    
    if 'catalog_file' not in params.keys():
        print(params)
        catalog = fetch_catalog_sources_within_image(params)
    
        catalog.write(path.join(params['red_dir'],'catalog.data'), 
                                format='ascii.basic', overwrite=True)
        
    else:

        catalog = Table.read(params['catalog_file'], format='ascii.basic')
    
    print(catalog)
    
    plot_catalog(params,catalog,lsst,wfirst)
    
def get_args():
    
    params = { }

    if len(argv) > 1:
        params['red_dir'] = argv[1]
        params['ra'] = argv[2]
        params['dec'] = argv[3]
        if len(argv) == 5:
            params['catalog_file'] = str(argv[4]).replace('-catalog=','')
            
    else:
        print('Parameters:')
        print('> python optimize_bulge_ddf_footprint.py red_dir  ra dec[sexigesimal] -catalog=file_path')
        exit()
    
    return params
    
def fetch_catalog_sources_within_image(params):
    """Function to extract the objects from the VPHAS+ catalogue within the
    field of view of the reference image, based on the metadata information."""
    
    # Radius should be in arcmin
    params['radius'] = (np.sqrt(params['fov'])/2.0)*60.0
    
    catalog = vizier_tools.search_vizier_for_sources(params['ra'], 
                                                       params['dec'], 
                                                        params['radius'], 
                                                        'VPHAS+')
    
    return catalog    

def plot_catalog(params,catalog,lsst,wfirst):
    """Function to plot the field of view"""
    
    c = SkyCoord(params['ra']+' '+params['dec'], unit=(u.hourangle, u.deg))
    
    w = WCS(naxis=2)
    w.wcs.crpix = [0.0, 0.0]
    w.wcs.cdelt = np.array([0.26/3600.0, 0.26/3600.0])
    w.wcs.crval = [c.ra.degree, c.dec.degree]
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    
    footprint = Circle((c.ra.degree/15.0, c.dec.degree), lsst.radius/2, 
                       edgecolor='yellow', facecolor='none', linewidth=5.0)

    #plt.subplot(projection=w)    
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    
    plt.subplots_adjust(left=None, bottom=0.2, 
                        right=None, top=None, 
                        wspace=None, hspace=None)
    try:
        plt.scatter(catalog['RAJ2000']/15.0, catalog['DEJ2000'], 
                s=(params['lsst_pixscale']/36000.0),
                edgecolor='black', facecolor=(0, 0, 0, 1.0))
    except KeyError:
        plt.scatter(catalog['_RAJ2000']/15.0, catalog['_DEJ2000'], 
                s=(params['lsst_pixscale']/36000.0),
                edgecolor='black', facecolor=(0, 0, 0, 1.0))
                
    ax.add_patch(footprint)
    
    wfirst.draw_footprint(ax, x=c.ra.degree/15.0, y=c.dec.degree)
    
    #plt.grid(color='white', ls='solid')
    plot_width = 2.0
    xmin = (c.ra.degree - 15.0*plot_width/2.0)/15.0
    xmax = (c.ra.degree + 15.0*plot_width/2.0)/15.0
    ymin = c.dec.degree - plot_width/2.0
    ymax = c.dec.degree + plot_width/2.0
    plt.axis([xmin,xmax,ymin,ymax])
    
    plt.ylabel('Dec')
    plt.xlabel('RA')

    plt.xticks(rotation=45.0)
    plt.yticks(rotation=45.0)
    
    plt.savefig(path.join(params['red_dir'],'sky_view.png'))
    
if __name__ == '__main__':
    
    optimize_survey_footprint()
    