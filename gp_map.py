import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from os import path
from astropy import units as u
from astropy_healpix import HEALPix
from astropy.coordinates import Galactic, TETE, SkyCoord

# Configuration
STAR_MAP_DIR = '/Users/rstreet1/software/LSST-TVSSC_software_tools/star_density_maps'
STAR_MAP_FILE = 'TRIstarDensity_r_nside_64.npz'

def load_star_density_map(limiting_mag=28.0):

    data_file = path.join(STAR_MAP_DIR, STAR_MAP_FILE)
    if path.isfile(data_file):
        npz_file = np.load(data_file)
        with np.load(data_file) as npz_file:
            star_map = npz_file['starDensity']
            mag_bins = npz_file['bins']

            dmag = abs(mag_bins - limiting_mag)
            idx = np.where(dmag == dmag.min())[0]

            star_density_map = np.copy(star_map[:,idx]).flatten()
            star_density_map = hp.reorder(star_density_map, n2r=True)

        return star_density_map

    else:
        raise IOError('Cannot find star density map data file at '+data_file)

    return None

def rotateHealpix(hpmap, transf=['C','G'], phideg=0., thetadeg=0.):
    """Rotates healpix map from one system to the other. Returns reordered healpy map.
    Healpy coord transformations are used, or you can specify your own angles in degrees.
    To specify your own angles, ensure that transf has length != 2.
    Original code by Xiaolong Li
    """

    # For reasons I don't understand, entering in ['C', 'G'] seems to do the
    # transformation FROM galactic TO equatorial. Possibly something buried in
    # the conventions used by healpy.

    # Heavily influenced by stack overflow solution here:
    # https://stackoverflow.com/questions/24636372/apply-rotation-to-healpix-map-in-healpy

    nside = hp.npix2nside(len(hpmap))

    # Get theta, phi for non-rotated map
    t,p = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))

    # Define a rotator
    if len(transf) == 2:
        r = hp.Rotator(coord=transf)
    else:
        r = hp.Rotator(deg=True, rot=[phideg,thetadeg])

    # Get theta, phi under rotated co-ordinates
    trot, prot = r(t,p)

    # Interpolate map onto these co-ordinates
    rot_map = hp.get_interp_val(hpmap, trot, prot)

    return rot_map

def calc_hp_pixels_for_region(l_center, b_center, l_width, b_height, n_points, ahp):

    halfwidth_l = l_width / 2.0
    halfheight_b = b_height / 2.0
    l = np.linspace(l_center-halfwidth_l, l_center+halfwidth_l, n_points) * u.deg
    b = np.linspace(b_center-halfheight_b, b_center+halfheight_b, n_points) * u.deg

    LL,BB = np.meshgrid(l, b)

    coords = SkyCoord(LL, BB, frame=Galactic())

    pixels = ahp.skycoord_to_healpix(coords)

    return pixels

def bono_survey_regions():

    n_points = 500
    l = np.linspace(-20.0, 20.0, n_points) * u.deg
    b = np.linspace(-15.0, 10.0, n_points) * u.deg
    LL,BB = np.meshgrid(l, b)
    coords = SkyCoord(LL, BB, frame=Galactic())
    shallow_pix = ahp.skycoord_to_healpix(coords)


    n_points = 100
    l = np.linspace(-20.0, 20.0, n_points) * u.deg
    b = np.linspace(-3.0, 3.0, n_points) * u.deg
    LL,BB = np.meshgrid(l, b)
    coords = SkyCoord(LL, BB, frame=Galactic())
    deep_pix = ahp.skycoord_to_healpix(coords)

    return shallow_pix, deep_pix

NSIDE = 64
NPIX = hp.nside2npix(NSIDE)
m = np.zeros(NPIX)

# Trilegal data is in Galactic coordinates, so this needs to be rotated in
# order to map it to healpix
star_density_map = load_star_density_map(limiting_mag=24.7)
hp_star_density = rotateHealpix(star_density_map)
hp_log_star_density = np.log10(hp_star_density)

ahp = HEALPix(nside=NSIDE, order='ring', frame=TETE())

# Galactic Plane survey regions, -85.0 < l (phi) <+85.0◦, -10.0 < b (theta) <+10.0◦
gp_region_pix1 = calc_hp_pixels_for_region(43.5, 0.0, 85.0, 20.0, 500, ahp)
gp_region_pix2 = calc_hp_pixels_for_region(317.5, 0.0, 85.0, 20.0, 500, ahp)
gp_region_pix = np.concatenate((gp_region_pix1.flatten(),gp_region_pix2.flatten()))

gp_region_pix1 = calc_hp_pixels_for_region(7.5, 0.0, 15.0, 20.0, 500, ahp)
gp_region_pix2 = calc_hp_pixels_for_region(352.5, 0.0, 15.0, 20.0, 500, ahp)
Gonzalez_gp_pix = np.concatenate((gp_region_pix1.flatten(),gp_region_pix2.flatten()))

(Bono_shallow_pix, Bono_deep_pix) = bono_survey_regions()

# Magellenic Clouds regions
# LMC  277.77 - 283.155, -35.17815 - -30.59865
LMC_pix = calc_hp_pixels_for_region(280.4652, -32.888443, (322.827/60), (274.770/60), 100, ahp)

# SMC 301.4908 - 304.126, -45.1036 - -43.5518
SMC_pix = calc_hp_pixels_for_region(302.8084, -44.3277, (158.113/60), (93.105/60), 100, ahp)

bulge_pix = calc_hp_pixels_for_region(2.216, -3.14, 3.5, 3.5, 50, ahp)

# Clementini survey regions
M54_pix = calc_hp_pixels_for_region(5.60703,-14.08715, 3.5, 3.5, 20, ahp)
Sculptor_pix = calc_hp_pixels_for_region(287.5334, -83.1568, 3.5, 3.5, 20, ahp)
Carina_pix = calc_hp_pixels_for_region(260.1124, -22.2235, 3.5, 3.5, 20, ahp)
Fornax_pix = calc_hp_pixels_for_region(237.1038, -65.6515, 3.5, 3.5, 20, ahp)
Phoenix_pix = calc_hp_pixels_for_region(272.1591, -68.9494, 3.5, 3.5, 20, ahp)
Antlia2_pix = calc_hp_pixels_for_region(264.8955, 11.2479, 3.5, 3.5, 20, ahp)
Clementini_regions = np.concatenate((M54_pix.flatten(), Sculptor_pix.flatten()))
for cluster in [Carina_pix, Fornax_pix, Phoenix_pix, Antlia2_pix]:
    Clementini_regions = np.concatenate((Clementini_regions, cluster.flatten()))

# Bonito survey regions
EtaCarina_pix = calc_hp_pixels_for_region(287.5967884538, -0.6295111793, 3.5, 3.5, 20, ahp)
OrionNebula_pix = calc_hp_pixels_for_region(209.0137, -19.3816, 3.5, 3.5, 20, ahp)
NGC2264_pix = calc_hp_pixels_for_region(202.9358, 2.1957, 3.5, 3.5, 20, ahp)
NGC6530_pix = calc_hp_pixels_for_region(6.0828, -01.3313, 3.5, 3.5, 20, ahp)
NGC6611_pix = calc_hp_pixels_for_region(16.9540, 0.7934, 3.5, 3.5, 20, ahp)
Bonito_regions = np.concatenate((EtaCarina_pix.flatten(), OrionNebula_pix.flatten()))
for cluster in [NGC2264_pix, NGC6530_pix, NGC6611_pix]:
    Bonito_regions = np.concatenate((Bonito_regions, cluster.flatten()))

# List of high-priority survey regions, in the form of HP pixels:
high_priority_regions = {'Galactic_Plane': gp_region_pix,
                         'Gonzalez_Plane_region': Gonzalez_gp_pix,
                         'Bono_shallow_survey': Bono_shallow_pix,
                         'Bono_deep_survey': Bono_deep_pix,
                         'Large_Magellenic_Cloud': LMC_pix,
                         'Small_Magellenic_Cloud': SMC_pix,
                         'Galactic_Bulge': bulge_pix,
                         'Clementini_regions': Clementini_regions,
                         'Bonito_regions': Bonito_regions}
regions_outside_plane = [LMC_pix, SMC_pix, Clementini_regions, Bonito_regions]

# Plotting scale range:
plot_min = 0.598
plot_max = 7.545

# Plot all regions separately for reference:
for name,region in high_priority_regions.items():
    map = np.zeros(NPIX)
    map[region] = hp_star_density[region]

    fig = plt.figure(1,(10,10))
    hp.mollview(np.log10(map), title=name,
                min=plot_min, max=plot_max)
    hp.graticule()
    plt.tight_layout()
    plt.savefig(name+'_footprint_map'+'.png')
    plt.close(1)

# Maximum footprint map:
max_footprint_map = np.zeros(NPIX)
for region_name, region_pix in high_priority_regions.items():
    max_footprint_map[region_pix] = hp_star_density[region_pix]

fig = plt.figure(1,(10,10))
hp.mollview(np.log10(max_footprint_map), title="Priority regions of the Galactic Plane - maximum footprint",
            min=plot_min, max=plot_max)
hp.graticule()
plt.tight_layout()
plt.savefig('max_GalPlane_footprint_map.png')
plt.close(1)

# Medium footprint map
density_thresholds = [ 0.75, 0.80 ]
for threshold in density_thresholds:
    medium_footprint_map = np.zeros(NPIX)
    idx = np.where(hp_log_star_density >= threshold*hp_log_star_density.max())[0]
    medium_footprint_map[idx] = hp_star_density[idx]
    for region_pix in regions_outside_plane:
        medium_footprint_map[region_pix] = hp_star_density[region_pix]

    plot_min = threshold * plot_max
    fig = plt.figure(2,(10,10))
    hp.mollview(np.log(medium_footprint_map),
                title="Priority regions of the Galactic Plane - medium footprint, "+\
                    str(round(threshold*100.0,0))+'% of max threshold')
    hp.graticule()
    plt.tight_layout()
    plt.savefig('medium_GalPlane_footprint_map_'+str(round(threshold*100.0,0))+'.png')
    plt.close(2)

# Minimum footprint map
minimum_footprint_map = np.zeros(NPIX)
regions = [Bono_deep_pix,
LMC_pix, SMC_pix, bulge_pix,
Clementini_regions, Bonito_regions]
for region_pix in regions:
    minimum_footprint_map[region_pix] = hp_star_density[region_pix]

fig = plt.figure(3,(10,10))
hp.mollview(np.log10(minimum_footprint_map), title="Priority regions of the Galactic Plane - minimum footprint")
hp.graticule()
plt.tight_layout()
plt.savefig('min_GalPlane_footprint_map.png')
plt.close(3)

plot_min = 0.598
plot_max = 7.545
fig = plt.figure(4,(10,10))
hp.mollview(np.log10(hp_star_density), title="Density of stars within Rubin viewing zone",
            min=plot_min, max=plot_max)
hp.graticule()
plt.tight_layout()
plt.savefig('rubin_star_density_map.png')
plt.close(4)
