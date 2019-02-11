import numpy as np
import ConfigParser
import os
import logging

def get_params(param_file_path):

    'Turn parameter file into parameter dictionary'

    config = ConfigParser.SafeConfigParser()
    config.read(param_file_path)

    # Get dictionaries from config object
    raw_params = dict(config.items('general'))
    raw_cats_params = dict(config.items('catalogs_to_stack'))
    raw_cats_path_params = dict(config.items('catalog_path'))
    raw_cats_file_params = dict(config.items('catalog_file'))
    raw_binning_params = dict(config.items('binning'))
    raw_maps_params = dict(config.items('maps_to_stack'))
    raw_map_path_params = dict(config.items('map_path'))
    raw_map_file_params = dict(config.items('map_file'))
    raw_noise_file_params = dict(config.items('noise_file'))
    raw_beam_params = dict(config.items('beams'))
    raw_colour_correction_params = dict(config.items('color_correction'))
    # Convert raw config dictionary to organisex dictionary params

    params = get_general_params(raw_params)
    params['map_files'] = get_maps_params(raw_maps_params,raw_map_path_params,raw_map_file_params)
    params['noise_files'] = get_maps_params(raw_maps_params,raw_map_path_params,raw_noise_file_params)
    params['wavelength'] = get_wavelength_params(raw_maps_params)
    params['psfs'] = get_beams_params(raw_maps_params,raw_beam_params)
    params['colour_correction'] = get_colour_correction_params(raw_maps_params,raw_colour_correction_params)
    params['catalogs'] = get_catalogs_params(raw_cats_params,raw_cats_path_params,raw_cats_file_params)
    params['bins'] = get_binning_params(raw_binning_params)
    params['library_keys'] = params['map_files'].keys()

    logging.info('------PARAMETER VALUES------')

    return params

def get_general_params(raw_params):

    params = {}
    params['zkey'] = raw_params['zkey']
    params['mkey'] = raw_params['mkey']
    params['ra_key'] = raw_params['ra_key']
    params['dec_key'] = raw_params['dec_key']
    try:
        params['save_bin_ids'] = string_is_true(raw_params['save_bin_ids'])
    except:
        params['save_bin_ids'] = True

    return params

def get_wavelength_params(raw_maps_params):

    wavelengths = {}
    for imap in raw_maps_params:
        if string_is_true(raw_maps_params[imap].split()[1]) == True:
            wavelengths[imap] = float(raw_maps_params[imap].split()[0])
    return wavelengths


def get_binning_params(raw_params):

    binning = {}
    z_nodes = []
    m_nodes = []
    for i in raw_params['redshift_nodes'].split():
        z_nodes.append(float(i))
    for j in raw_params['mass_nodes'].split():
        m_nodes.append(float(j))

    binning['z_nodes'] = z_nodes
    binning['m_nodes'] = m_nodes

    return binning

def get_maps_params(raw_maps_params,raw_map_path_params,raw_map_file_params):

    maps = {}

    for imap in raw_maps_params:
        if string_is_true(raw_maps_params[imap].split()[1]) == True:
            maps[imap] = raw_map_path_params[imap].split()[0] + raw_map_file_params[imap]

    return maps

def get_beams_params(raw_maps_params,raw_beam_params):
    psfs = {}

    for imap in raw_maps_params:
        if string_is_true(raw_maps_params[imap].split()[1]) == True:
            psfs[imap + '_beam_area'] = float(raw_beam_params[imap].split()[1])
            if is_float(raw_beam_params[imap].split()[0]) == True:
                psfs[imap+'_fwhm'] = float(raw_beam_params[imap].split()[0])
            else:
                psfs[imap+'_beam_file'] = raw_beam_params[imap].split()[0]
    return psfs


def get_colour_correction_params(raw_maps_params,raw_colour_correction_params):
    colour_correction = {}

    for imap in raw_maps_params:
        if string_is_true(raw_maps_params[imap].split()[1]) == True:
            colour_correction[imap+''] = float(raw_colour_correction_params[imap])

    return colour_correction


def get_catalogs_params(raw_cats_params,raw_cats_path_params,raw_cats_file_params):
    catalog = {}

    for cat in raw_cats_params:
        if string_is_true(raw_cats_params[cat].split()[0]) == True:
            catalog['catalog_path'] = str(raw_cats_path_params[cat].split()[0])
            catalog['catalog_file'] = str(raw_cats_file_params[cat])
            catalog['cat_name'] = str(cat)

    return catalog

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def string_is_true(sraw):

    s = sraw.lower()
    true_strings = ['true','t','yes','y','1']
    false_strings = ['false','f','no','n','0']
    if s in true_strings:
        return True
    elif s in false_strings:
        return False
    else:
        logging.warning('Input not recognised for parameter %s' % (key))
        logging.warning('You providedL %s' % (sraw))
        raise


    
    
