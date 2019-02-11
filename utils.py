import pdb
import gc
import numpy as np
from numpy import zeros
from numpy import shape
#from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import Planck15 as cosmo
import astropy.units as u
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
import scipy.io
import pylab as plt
import pandas as pd
from astropy.table import Table
from itertools import combinations_with_replacement

pi=3.141592653589793


## A



## B



def clean_args(dirty_args):
    return dirty_args.replace('.','p').replace('-','_')

def clean_arrays(x_array, y_array, z_array=None):
    xout = []
    yout = []
    if z_array != None:
        zout = []
    for i in range(len(y_array)):
        if y_array[i] != 0:
            if np.sum(np.isnan(x_array[i])) > 0:
                #print 'nan!'
                pass
            else:
                yout.append(y_array[i])
                xout.append(x_array[i])
                if z_array != None:
                    zout.append(z_array[i])
    if z_array != None:
        return np.array(xout),np.array(yout),np.array(zout)
    else:
        return np.array(xout),np.array(yout)

def clean_nans(dirty_array, replacement_char=0.0):
  clean_array = dirty_array
  clean_array[np.isnan(dirty_array)] = replacement_char
  clean_array[np.isinf(dirty_array)] = replacement_char

  return clean_array








def gauss(x, x0, y0, sigma):
    p = [x0, y0, sigma]
    return p[1]* np.exp(-((x-p[0])/p[2])**2)

def gauss_kern(fwhm, side, pixsize):
  ''' Create a 2D Gaussian (size= side x side)'''

  sig = fwhm / 2.355 / pixsize
  delt = zeros([int(side),int(side)])
  delt[0,0]=1.0
  ms = shape(delt)
  delt = shift_twod(delt, ms[0] / 2, ms[1] / 2)
  kern = delt
  gaussian_filter(delt, sig, output= kern)
  kern /= np.max(kern)

  #pdb.set_trace()
  return kern










## M

def map_rms(map,header=None,mask=None,silent=True):
    if mask != None:
         ind = np.where((mask == 1) & (clean_nans(map) != 0))
         print 'using mask'
    else:
         ind = clean_nans(map) != 0
    map /= np.max(map)

    #hist, bin_edges = np.histogram(map[ind], density=True)
    #hist, bin_edges = np.histogram(map[ind],range=(np.min(map),0),bins=50)
    #hist, bin_edges = np.histogram(map[ind],range=(np.min(map),abs(np.min(map))),bins=50,density=True)
    #x0 = 0.9*np.min(map)
    x0 = abs(np.percentile(map,99))
    #hist, bin_edges = np.histogram(np.unique(map),range=(np.min(map),abs(np.min(map))),bins=50,density=True)
    hist, bin_edges = np.histogram(np.unique(map),range=(-x0,x0),bins=30,density=True)

    p0 = [0., 1., x0/3]
    x = .5 * (bin_edges[:-1] + bin_edges[1:])
    #x_peak = x[hist == max(hist)][0]
    x_peak = 1+np.where((hist - max(hist))**2 < 0.01)[0][0]
    #x_peak = find_nearest_index(hist, max(hist)[0])

    # Fit the data with the function
    #fit, tmp = curve_fit(gauss, x, hist/max(hist), p0=p0)
    fit, tmp = curve_fit(gauss, x[:x_peak], hist[:x_peak]/max(hist), p0=p0)
    #sig_rad = fit[2] * pixsize_deg * (3.14159 / 180)
    #fwhm = fit[2] * pixsize_deg * 3600. * 2.355
    rms_1sig = abs(fit[2])
    if silent == False:
        print('1sigma rms=%.2e' % rms_1sig)
        plt.plot(x,hist)
        plt.plot(x[:x_peak],hist[:x_peak])
        plt.plot(np.linspace(-abs(x0),abs(x0),121),
                max(hist)*gauss(np.linspace(-abs(x0),abs(x0),121),*fit),'m--')
        plt.show()
    #pdb.set_trace()

    return rms_1sig


## P
def pad_and_smooth_psf(mapin, psfin):

  s = np.shape(mapin)
  mnx = s[0]
  mny = s[1]

  s = np.shape(psfin)
  pnx = s[0]
  pny = s[1]

  psf_x0 = pnx/2
  psf_y0 = pny/2
  psf = psfin
  px0 = psf_x0
  py0 = psf_y0

  # pad psf
  psfpad = np.zeros([mnx, mny])
  psfpad[0:pnx,0:pny] = psf

  # shift psf so that centre is at (0,0)
  psfpad = shift_twod(psfpad, -px0, -py0)
  smmap = np.real( np.fft.ifft2( np.fft.fft2(zero_pad(mapin) ) *
    np.fft.fft2(zero_pad(psfpad)) ) )

  return smmap[0:mnx,0:mny]



def shift(seq, x):
  from numpy import roll
  out = roll(seq, int(x))
  return out

def shift_twod(seq, x, y):
  from numpy import roll
  out = roll(roll(seq, int(x), axis = 1), int(y), axis = 0)
  return out

def shift_bit_length(x):
  return 1<<(x-1).bit_length()

def smooth_psf(mapin, psfin):

  s = np.shape(mapin)
  mnx = s[0]
  mny = s[1]

  s = np.shape(psfin)
  pnx = s[0]
  pny = s[1]

  psf_x0 = pnx/2
  psf_y0 = pny/2
  psf = psfin
  px0 = psf_x0
  py0 = psf_y0

  # pad psf
  psfpad = np.zeros([mnx, mny])
  psfpad[0:pnx,0:pny] = psf

  # shift psf so that centre is at (0,0)
  psfpad = shift_twod(psfpad, -px0, -py0)
  smmap = np.real( np.fft.ifft2( np.fft.fft2(mapin) *
    np.fft.fft2(psfpad))
    )

  return smmap


def string_is_true(sraw):
    """Is string true? Returns boolean value.
    """
    s       = sraw.lower() # Make case-insensitive

    # Lists of acceptable 'True' and 'False' strings
    true_strings    = ['true', 't', 'yes', 'y', '1']
    false_strings    = ['false', 'f', 'no', 'n', '0']
    if s in true_strings:
        return True
    elif s in false_strings:
        return False
    else:
        logging.warning("Input not recognized for parameter: %s" % (key))
        logging.warning("You provided: %s" % (sraw))
        raise

def solid_angle_from_fwhm(fwhm_arcsec):
  sa = np.pi*(fwhm_arcsec / 3600.0 * np.pi / 180.0)**2.0 / (4.0 * np.log(2.0))
  return sa






def shift_ra_dec_values(ras,decs,ra_shift,dec_shift):
    'Shifts the ra and dec values to create lists for plot'
    ras_new = [ra + ra_shift for ra in ras]
    decs_new = [dec + dec_shift for dec in decs]
    return ras_new,decs_new

def shift_RAs(ras,ra_shift):
    ras_new = [ra + ra_shift for ra in ras]
    return ras_new

def shift_DECs(decs,dec_shift):
    decs_new = [dec + dec_shift for dec in decs]
    return decs_new


def find_image_indexes(shifted_list,pixsize,side_length):
    matrix_indexes = {}
    for ra_shift,dec_shift in shifted_list:
        ra_initial = ra_shift/pixsize
        dec_initial = dec_shift/pixsize
        ra_index = (side_length + 1)/2 - ra_initial - 1 # flipped ra and dec index to test
        dec_index = dec_initial + (side_length + 1)/2  -1
        matrix_indexes[(ra_shift,dec_shift)] = (int(ra_index),int(dec_index))
    return matrix_indexes

def create_shifted_list(side_length,pixsize):
    co_ords = []
    for i in range(side_length):
        if len(co_ords) > 0:
            pos = 0 + i
            neg = 0 - i
            co_ords.append(pos)
            co_ords.append(neg)
        else:
            co_ords.append(0)
        if len(co_ords) == side_length:
            break
    shift_co_ords = [co_ord*pixsize for co_ord in co_ords]
    shifted_list = list(combinations_with_replacement(shift_co_ords,2))
    new_list = []
    for x,y in shifted_list:
        new_list.append((x,y))
        new_list.append((y,x))
    cleaned_list = list(set(new_list))
    return cleaned_list


def get_stacked_statistics(imap):
    mean = np.mean(imap)
    median = np.median(imap)
    var = np.var(imap)
    std = np.std(imap)
    summed = np.sum(imap)
    return mean,median,var,std,summed

def get_gaussian_prob(mean,sigma,x):

  A = 1/(sigma*np.sqrt(2*np.pi))
  b = -(x - mean)**2
  f = A * np.exp(b/(2*sigma**2))
  return f


def get_stacked_image_data_table(data_dict,shifted_list,param_names):

  
  for iparam in param_names:
    means = np.array([data_dict[(ra_shift,dec_shift)][iparam]['mean'] for ra_shift,dec_shift in shifted_list])
    medians = np.array([data_dict[(ra_shift,dec_shift)][iparam]['median'] for ra_shift,dec_shift in shifted_list])
    sums = np.array([data_dict[(ra_shift,dec_shift)][iparam]['sum'] for ra_shift,dec_shift in shifted_list])
    stdevs = np.array([data_dict[(ra_shift,dec_shift)][iparam]['stdev'] for ra_shift,dec_shift in shifted_list])
  x_positions = np.array([ra_shift for ra_shift,dec_shift in shifted_list])
  y_positions = np.array([dec_shift for ra_shift,dec_shift in shifted_list])

  table = Table([x_positions,y_positions,means,medians,sums,stdevs],names=('Ra_shift','dec_shift','mean','median','sum','Stdev'))
  print table
  return table


def get_stacked_image_data_table_sams_stacking(data_dict,shifted_list,bin_names):

  for name in bin_names:
    stacked_fluxes = np.array([data_dict[(ra_shift,dec_shift)][name]['value'] for ra_shift,dec_shift in shifted_list])

  x_positions = np.array([ra_shift for ra_shift,dec_shift in shifted_list])
  y_positions = np.array([dec_shift for ra_shift,dec_shift in shifted_list])

  table = Table([x_positions,y_positions,stacked_fluxes],names=('Ra shift','dec shift','stacked flux'))

  return table
