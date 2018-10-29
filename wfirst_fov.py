# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 10:40:04 2018

@author: rstreet
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from sys import argv

class WFIRSTFootprint:
    
    def __init__(self):
        polygons = [ np.array([[-0.01981411, -0.00567634], 
                            [-0.08244734, 0.10480327], 
                            [0.02616297, 0.16866906], 
                            [0.08929735, 0.05732157], 
                            [-0.01981411, -0.00567634]]), 
                           np.array([[-0.13671453, -0.07317132], 
                            [-0.19971191, 0.03593987], 
                            [-0.09333858, 0.09966936],
                            [-0.03020550, -0.01167605],
                            [-0.13671453, -0.07317132]]),
                           np.array([[-0.26226768, -0.14566156],
                            [-0.32662835, -0.03618749], 
                            [-0.22185914, 0.02430662], 
                            [-0.15836219, -0.08567004], 
                            [-0.26226768, -0.14566156]]), 
                           np.array([[-0.05973719, 0.12946287], 
                            [-0.12150469, 0.24044163], 
                            [-0.01289466, 0.30430971], 
                            [0.05024008, 0.19296323], 
                            [-0.05973719, 0.12946287]]), 
                           np.array([[-0.17663707, 0.06196459], 
                            [-0.23840273, 0.17294098], 
                            [-0.13152991, 0.23580750],
                            [-0.06976251, 0.12482890], 
                            [-0.17663707, 0.06196459]]), 
                           np.array([[-0.30355509, -0.01016407], 
                            [-0.36431821, 0.09907748], 
                            [-0.26054970, 0.16130714], 
                            [-0.19928448, 0.05119708], 
                            [-0.30355509, -0.01016407]]), 
                           np.array([[-0.07281623,0.28010224],
                            [-0.13458423, 0.39107907], 
                            [-0.02510754, 0.45544909],
                            [0.03666190, 0.34447083], 
                            [-0.07281623, 0.28010224]]), 
                           np.array([[-0.18935090, 0.21396655], 
                            [-0.24975119, 0.32457547], 
                            [-0.14337765, 0.38831066], 
                            [-0.08160972, 0.27733403],
                            [-0.18935090, 0.21396655]]), 
                           np.array([[-0.31713549, 0.14133313],
                            [-0.37703327, 0.25107343], 
                            [-0.27326421, 0.31330682],
                            [-0.21249852, 0.20406430], 
                            [-0.31713549, 0.14133313]]),
                           np.array([[0.04955101, -0.12381585], 
                            [-0.01481417, -0.01433617],
                            [0.09429738, 0.04866159],
                            [0.15916344, -0.06168604],
                            [0.04955101, -0.12381585]]), 
                           np.array([[-0.06821700, -0.19180802], 
                            [-0.13171472, -0.08183094],
                            [-0.02520557, -0.02033586],
                            [0.03965952, -0.13068131],
                            [-0.06821700, -0.19180802]]), 
                           np.array([[-0.19427241, -0.26342878], 
                            [-0.25813391, -0.15482087],
                            [-0.15336241, -0.09432960], 
                            [-0.08986497, -0.20430612],
                            [-0.19427241, -0.26342878]]),
                           np.array([[0.14749294, -0.22545587], 
                            [0.08226173, -0.11647697],
                            [0.19224104, -0.05298073], 
                            [0.25760683, -0.16419344],
                            [0.14749294, -0.22545587]]),
                           np.array([[0.03058862, -0.29294629], 
                            [-0.03464149, -0.18396990], 
                            [0.07323614, -0.12284251], 
                            [0.13896730, -0.23268714], 
                            [0.03058862, -0.29294629]]),
                           np.array([[-0.09633612, -0.36506486], 
                            [-0.16056453, -0.25782348], 
                            [-0.05528978, -0.19820006], 
                            [0.00944017, -0.30630999], 
                            [-0.09633612, -0.36506486]]),
                           np.array([[0.27141457, -0.31209584], 
                            [0.20618352, -0.20311872],
                            [0.31666463, -0.14048928], 
                            [0.38189563, -0.24946824], 
                            [0.27141457, -0.31209584]]),
                           np.array([[0.15364197, -0.38008541],
                            [0.09324129, -0.26947695], 
                            [0.20112166, -0.20835031], 
                            [0.26288872, -0.31932701], 
                            [0.15364197, -0.38008541]]),
                           np.array([[0.02894538, -0.45206931], 
                            [-0.03614990, -0.34532926], 
                            [0.07049446, -0.28607275], 
                            [0.13472484, -0.39331510], 
                            [0.02894538, -0.45206931]]) ]
        
        self.footprint = []
        
        for ccd in polygons:

            ccd[:,0] = ccd[:,0] / 1.324
            ccd[:,1] = ccd[:,1] / 1.088
            
            self.footprint.append(ccd)
            
    def draw_footprint(self,ax,x,y):
        """Method to draw a WFIRST-WFI footprint on a pre-existing figure Axis,
        centered at the x, y coordinates given."""
        
        for ccd in self.footprint:
            
            ccd_on_sky = ccd.copy()
            
            ccd_on_sky[:,0] = ccd_on_sky[:,0] + x
            ccd_on_sky[:,1] = ccd_on_sky[:,1] + y
            
            ax.add_patch( patches.Polygon( ccd_on_sky, fill=False ) )

        return ax

def plot_wfi_footprint(x=0.0,y=0.0):
    """Function to plot the footprint of WFIRST's WFI"""

    fig = plt.figure()

    ax = fig.add_subplot(111, aspect='equal')
    
    wfirst = WFIRSTFootprint()
    
    wfirst.draw_footprint(ax,x,y)

    rot_field = rotate_fov()
    
    rot_field.draw_footprint(ax,x,y)
    
    plt.axis([-1.0,1.0,-1.0,1.0])
    
    plt.grid()
    
    fig.savefig('wfirst_wfi_footpring.png', dpi=90, bbox_inches='tight')

def rotate_fov():
    
    rot = -30.0
    ref_theta = 210.0 * np.pi/180.0
    dtheta = rot * np.pi/180.0
    
    r = 0.491
    wfirst = WFIRSTFootprint()
    
    detector = []
    
    xr = [1e6,-1e6]
    yr = [1e6,-1e6]
    
    for ccd in wfirst.footprint:

        xprime = ccd[:,0]*np.cos(dtheta) - ccd[:,1]*np.sin(dtheta)
        yprime = ccd[:,0]*np.sin(dtheta) + ccd[:,1]*np.cos(dtheta)
        
        rot_ccd = np.zeros(ccd.shape)
        
        rot_ccd[:,0] = xprime
        rot_ccd[:,1] = yprime
        
        detector.append(rot_ccd)
        
        if rot_ccd[:,0].min() < xr[0]:
            xr[0] = rot_ccd[:,0].min()
        if rot_ccd[:,0].max() > xr[1]:
            xr[1] = rot_ccd[:,0].max()
        if rot_ccd[:,1].min() < yr[0]:
            yr[0] = rot_ccd[:,1].min()
        if rot_ccd[:,1].max() > yr[1]:
            yr[1] = rot_ccd[:,1].max()
            
    wfirst.footprint = detector
    
    print('X range ',xr,(xr[1]-xr[0]))
    print('Y range ',yr,(yr[1]-yr[0]))
    
    return wfirst
    
if __name__ == '__main__':
    
    if len(argv) == 1:
        x = 0.0
        y = 0.0
    else:
        x = float(argv[1])
        y = float(argv[2])
    
    plot_wfi_footprint(x=x,y=y)
    