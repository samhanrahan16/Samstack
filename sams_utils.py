import numpy as np
from sams_skymaps import sams_skymap
from astropy.io import ascii
from catalog_binning import Field_catalog
from astropy.constants import c,h,k_B

def pack_stacked_fluxes(stacked_fluxes,stacked_noise,stacked_stats,bin_names,ngals,map_rms,names_to_mids):

    packed_stacks  ={}
    for name in bin_names:
        packed_stacks[name] = {}
        packed_stacks[name]['bin_midpoint'] = names_to_mids[name]
        packed_stacks[name]['value'] = stacked_fluxes[name]
        packed_stacks[name]['ngals_bin'] = ngals[name]
        packed_stacks[name]['psnerr'] = map_rms/np.sqrt(float(ngals[name]))
        packed_stacks[name]['mean'] = stacked_stats[name]['mean']
        packed_stacks[name]['median'] = stacked_stats[name]['median']
        packed_stacks[name]['stdev'] = stacked_stats[name]['stdev']
        packed_stacks[name]['noise'] = stacked_noise[name]

    return packed_stacks


def circle_around_galaxy(pixmap,rad,pixsize,ind_gals):
    radius = rad/pixsize
    x_shifts = range(9)
    y_shifts = range(9)
    co_ords = []
    for x in x_shifts:
        for y in y_shifts:
            if x**2 + y**2 <= radius**2:
                co_ords.append((x,y))

    ind_around = ([],[])
    vals_to_add = []

    for x,y in co_ords:
        for x_gal_loc,y_gal_loc in zip(ind_gals[0],ind_gals[1]):
            add_x = x_gal_loc + x
            add_y = y_gal_loc + y
            val_to_add = pixmap[x_gal_loc,y_gal_loc]
            vals_to_add.append(val_to_add)
            ind_around[0].append(add_x)
            ind_around[1].append(add_y)

    for x,y,val_to_add in zip(ind_around[0],ind_around[1],vals_to_add):
        if pixmap[x,y] == 0:
            pixmap[x,y] += val_to_add

    return pixmap,ind_around


def get_maps(params):
    '''Read the maps and store them in dictionaries'''
    sky_library = {}

    for t in params['library_keys']:
        sky = sams_skymap(params['map_files'][t],params['noise_files'][t],params['psfs'][t+'_fwhm'],colour_correction=params['colour_correction'][t],wavelength=params['wavelength'][t],fwhm = params['psfs'][t+'_fwhm'])
        sky_library[t] = sky

    return sky_library

def get_catalogs(params):

    catalog_path = params['catalogs']['catalog_path'] + params['catalogs']['catalog_file']

    tbl = ascii.read(catalog_path)

    zkey = params['zkey']
    mkey = params['mkey']
    rkey = params['ra_key']
    dkey = params['dec_key']

    catout = Field_catalog(tbl,params,zkey,mkey,rkey,dkey)

    return catout

def get_bin_ids(params,cats):

    catalog = get_catalogs(params)

    bin_names = catalog.make_bin_names()

    return bin_names

def get_bin_radec(params,cats):

    catalog = get_catalogs(params)

    binned_ra_dec = catalog.split_catalog_by_mass()

    return binned_ra_dec

def get_bin_names_to_midpoints(params,cats):

    catalog = get_catalogs(params)
    catalog.split_catalog_by_mass()

    names_to_mids = catalog.bin_names_to_midpoints
    return names_to_mids


def beam_solid_angle(fwhm):
    fwhm_deg = fwhm/3600.
    fwhm_rad = fwhm_deg*(2*np.pi)/360
    solid_angle = (fwhm_rad**2)*np.pi/(4*np.log(2))

    return solid_angle

def specific_intensity(flux_density,solid_angle):
    '''Flux density in Jy'''
    intensity = flux_density*(10**(-26))/solid_angle

    return intensity 

def blackbody_temp(wavelength,flux,beam_fwhm):
    '''Takes wavelength in micrometers and flux in Jy'''
    beam_solid_angle = beam_solid_angle(beam_fwhm)
    intensity = specific_intensity(flux,beam_solid_angle)
    wavelength_m = wavelength * (10**(-6))

    lnf = (2*h*c**2)/(intensity*wavelength_m**5)

    ln = np.log(1 + lnf)
    A = h*c/(wavelength_m*k_B)

    T = A /ln

    return T
    
    

    
