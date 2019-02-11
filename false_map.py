#!/usr/bin/env python

import numpy as np
import pandas as pd
import parameters
import os
import os.path
import random
from astropy.io import fits
from astropy import wcs
from astropy.table import Table
from astropy.io import ascii
from utils import get_gaussian_prob
from utils import create_shifted_list
from bincatalogs import Field_catalogs

def get_catalogs(params):

    # Formatting no longer needed as
    tbl = pd.read_table(params['catalogs']['catalog_path']+params['catalogs']['catalog_file'],sep=',')
    #print tbl

    tbl['ID'] = range(len(tbl))
    if 'sfg' in tbl.keys():
        pass
    elif 'CLASS' in tbl.keys():
        tbl['sfg']=tbl['CLASS']

    zkey = params['zkey']
    mkey = params['mkey']
    rkey = params['ra_key']
    dkey = params['dec_key']
    #print list(tbl[rkey])
    catout = Field_catalogs(tbl,zkey=zkey,mkey=mkey,rkey=rkey,dkey=dkey)

    return tbl

def gaussian_spread(distance,beam_sigma):
    exp_term = distance**2/(2*beam_sigma**2)
    f = np.exp(-exp_term)
    return f


class fake_map:
    '''Class to create a fake map with know flux values at each location to test the stacking code'''

    def __init__(self,file_map,file_noise,params,output_file_location,random_cat_location):

        self.original_map = file_map
        self.original_noise = file_noise
        self.rkey = params['ra_key']
        self.dkey = params['dec_key']
        self.zkey = params['zkey']
        self.mkey = params['mkey']
        self.cat = get_catalogs(params)
        #print self.cat
        self.RAs = list(self.cat[self.rkey])
        self.DECs = list(self.cat[self.dkey])
        self.masses = list(self.cat[self.mkey])
        self.zs = list(self.cat[self.zkey])
        self.new_map_location = output_file_location
        self.random_cat_location = random_cat_location
        

    def create_fake_matrix_positions(self,cat_number):

        hdu_list = fits.open(self.original_map)
        original_header = hdu_list[0].header
        hdu_list_noise = fits.open(self.original_noise)
        map_data = hdu_list[0].data
        noise_data = hdu_list_noise[0].data
        map_dimensions = map_data.shape
        flattened_noise = np.ndarray.flatten(noise_data)
        sigma_noise = np.absolute(np.mean(flattened_noise))
        random_map_numbers = np.random.normal(0,sigma_noise,len(flattened_noise))
        pixsize = hdu_list[0].header['CDELT2'] * 3600
        new_image_data = np.matrix(random_map_numbers.reshape(map_dimensions),dtype=np.float32)
        wcs_1 = wcs.WCS(self.original_map)
        pys,pxs = wcs_1.wcs_world2pix(self.RAs,self.DECs,0)
        beam_FWHM = 17.6
        beam_sigma = beam_FWHM/(2*np.sqrt(2*np.log(2)))


        # CREATE SHIFTED LIST FOR PIXEL POSITIONS
        shifted_list = create_shifted_list(5,pixsize)
        pixel_shifted_list = []
        fluxes_from_centre = {}
        for x_shift,y_shift in shifted_list:
            x_pixel = x_shift/pixsize
            y_pixel = y_shift/pixsize
            pixel_shifted_list.append((x_pixel,y_pixel))
            distance = np.sqrt(x_shift**2 + y_shift**2)
            f = gaussian_spread(distance,beam_sigma)
            flux = f*sigma_noise
            fluxes_from_centre[(x_pixel,y_pixel)] = flux
        print sigma_noise, 'mean noise value'
       # print new_image_data
       # print new_image_data[500,600]
        

        #CHANGE THE FLUX OF THE PIXEL LOCATIONS WHERE GALAXIES ARE

        for py,px in zip(pys,pxs):
            for x_pixel,y_pixel in pixel_shifted_list:
                x_pixel_location = int(x_pixel + px)
                y_pixel_location = int(y_pixel + py)
               # print x_pixel_location
               # print px, 'px'
               # print x_pixel, 'x_pixel'
                new_image_data[x_pixel_location,y_pixel_location] += fluxes_from_centre[(x_pixel,y_pixel)]

        # think I want to make the original array of noise into fits file then can add in fluxes of galaxy locations

        #CREATE THE IMAGE FILE
        w = wcs.WCS(naxis=2)
        #print original_header
        w.wcs.crpix = [original_header['CRPIX1'],original_header['CRPIX2']]
        w.wcs.cdelt = np.array([-pixsize/3600,pixsize/3600])
        w.wcs.crval = [original_header['CRVAL1'],original_header['CRVAL2']]
        w.wcs.ctype = [original_header['CTYPE1'],original_header['CTYPE2']]
       # w.wcs.crpix = [-1134.0,-189.0]
        #w.wcs.cdelt = np.array([-pixsize/3600.0,pixsize/3600.0])
        #w.wcs.crval = [197.99294599559,26.090059794022]
        #w.wcs.ctype = ['RA---TAN','DEC---TAN']
        header_w = w.to_header()
        hdu = fits.PrimaryHDU(data = new_image_data,header = header_w)
        #header['BITPIX'] = -32
        #header['CTYPE1'] = 'RA---TAN'
        #header['CTYPE2'] = 'DEC---TAN'
        #header['CRVAL1'] = 197.99294599559
        #header['CRVAL2'] = 26.090059794022
        #header['CRPIX1'] = -1134.0
        #header['CRPIX2'] = -189.0
        #header['CDELT1'] = -pixsize/3600
        #header['CDELT2'] = pixsize/3600
        #header['CUNIT1'] = 'deg'
        #header['CUNIT2'] = 'deg'
        #header['LONPOLE'] = 180.0
        #header['LATPOLE'] = 26.090059794022
        #header['WCSAXES'] = 2
        #header['RADESYS'] = 'ICRS'
        file_path = self.new_map_location + '/FAKE_MAP_POSITIONS_' + str(cat_number) + '.FITS'
        if os.path.exists(file_path) == True:
            os.remove(file_path)
            hdu.writeto(file_path)
        else:
            hdu.writeto(file_path)


    def create_fake_matrix_positions_random(self,cat_number):

        hdu_list = fits.open(self.original_map)
        original_header = hdu_list[0].header
        hdu_list_noise = fits.open(self.original_noise)
        map_data = hdu_list[0].data
        noise_data = hdu_list_noise[0].data
        map_dimensions = map_data.shape
        flattened_noise = np.ndarray.flatten(noise_data)
        sigma_noise = np.absolute(np.mean(flattened_noise))
        random_map_numbers = np.random.normal(0,sigma_noise,len(flattened_noise))
        pixsize = hdu_list[0].header['CDELT2'] * 3600
        new_image_data = np.matrix(random_map_numbers.reshape(map_dimensions),dtype=np.float32)
        wcs_1 = wcs.WCS(self.original_map)
        random_pxs = random.sample(range(0,map_dimensions[0]-20),len(self.RAs))
        random_pys = random.sample(range(0,map_dimensions[1]-20),len(self.RAs))
        beam_FWHM = 17.6
        beam_sigma = beam_FWHM/(2*np.sqrt(2*np.log(2)))

        # CREATE CATALOG OF RANDOM POSITIONS
        ras,decs = wcs_1.wcs_pix2world(random_pys,random_pxs,0)
        classes = [1] * len(ras)
        reffs = [4.0] * len(ras)
        data_table = Table([ras,decs,reffs,self.masses,classes,self.zs],names=[self.rkey,self.dkey,'Reff[arcsec]',self.mkey,'CLASS',self.zkey])
        print data_table
        output_cats_path = self.random_cat_location + '/cat_of_random_positions_' + str(cat_number) + '.csv'
        if os.path.exists(output_cats_path) == True:
            os.remove(output_cats_path)
            ascii.write(data_table,output_cats_path,format = 'csv')
        else:
            ascii.write(data_table,output_cats_path,format = 'csv')
            
        


        # CREATE SHIFTED LIST FOR PIXEL POSITIONS
        shifted_list = create_shifted_list(5,pixsize)
        pixel_shifted_list = []
        fluxes_from_centre = {}
        for x_shift,y_shift in shifted_list:
            x_pixel = x_shift/pixsize
            y_pixel = y_shift/pixsize
            pixel_shifted_list.append((x_pixel,y_pixel))
            distance = np.sqrt(x_shift**2 + y_shift**2)
            f = gaussian_spread(distance,beam_sigma)
            flux = f*sigma_noise
            fluxes_from_centre[(x_pixel,y_pixel)] = flux
        print sigma_noise, 'mean noise value'
       # print new_image_data
       # print new_image_data[500,600]
        

        #CHANGE THE FLUX OF THE PIXEL LOCATIONS WHERE GALAXIES ARE

        for py,px in zip(random_pys,random_pxs):
            for x_pixel,y_pixel in pixel_shifted_list:
                x_pixel_location = int(x_pixel + px)
                y_pixel_location = int(y_pixel + py)
               # print x_pixel_location
               # print px, 'px'
               # print x_pixel, 'x_pixel'
                new_image_data[x_pixel_location,y_pixel_location] += fluxes_from_centre[(x_pixel,y_pixel)]

        # think I want to make the original array of noise into fits file then can add in fluxes of galaxy locations

        #CREATE THE IMAGE FILE
        w = wcs.WCS(naxis=2)
        #print original_header
        w.wcs.crpix = [original_header['CRPIX1'],original_header['CRPIX2']]
        w.wcs.cdelt = np.array([-pixsize/3600,pixsize/3600])
        w.wcs.crval = [original_header['CRVAL1'],original_header['CRVAL2']]
        w.wcs.ctype = [original_header['CTYPE1'],original_header['CTYPE2']]
       # w.wcs.crpix = [-1134.0,-189.0]
        #w.wcs.cdelt = np.array([-pixsize/3600.0,pixsize/3600.0])
        #w.wcs.crval = [197.99294599559,26.090059794022]
        #w.wcs.ctype = ['RA---TAN','DEC---TAN']
        header_w = w.to_header()
        hdu = fits.PrimaryHDU(data = new_image_data,header = header_w)
        #header['BITPIX'] = -32
        #header['CTYPE1'] = 'RA---TAN'
        #header['CTYPE2'] = 'DEC---TAN'
        #header['CRVAL1'] = 197.99294599559
        #header['CRVAL2'] = 26.090059794022
        #header['CRPIX1'] = -1134.0
        #header['CRPIX2'] = -189.0
        #header['CDELT1'] = -pixsize/3600
        #header['CDELT2'] = pixsize/3600
        #header['CUNIT1'] = 'deg'
        #header['CUNIT2'] = 'deg'
        #header['LONPOLE'] = 180.0
        #header['LATPOLE'] = 26.090059794022
        #header['WCSAXES'] = 2
        #header['RADESYS'] = 'ICRS'
        file_path = self.new_map_location + '/FAKE_MAP_RANDOM_POSITIONS_' + str(cat_number) + '.FITS'
        if os.path.exists(file_path) == True:
            os.remove(file_path)
            hdu.writeto(file_path)
        else:
            hdu.writeto(file_path)

    def create_fake_matrix_no_sources(self):

        hdu_list = fits.open(self.original_map)
        original_header = hdu_list[0].header
        hdu_list_noise = fits.open(self.original_noise)
        map_data = hdu_list[0].data
        noise_data = hdu_list_noise[0].data
        map_dimensions = map_data.shape
        flattened_noise = np.ndarray.flatten(noise_data)
        sigma_noise = np.absolute(np.mean(flattened_noise))
        random_map_numbers = np.random.normal(0,sigma_noise,len(flattened_noise))
        pixsize = hdu_list[0].header['CDELT2'] * 3600
        new_image_data = np.matrix(random_map_numbers.reshape(map_dimensions),dtype=np.float32)
        wcs_1 = wcs.WCS(self.original_map)
        pys,pxs = wcs_1.wcs_world2pix(self.RAs,self.DECs,0)
        beam_FWHM = 17.6
        beam_sigma = beam_FWHM/(2*np.sqrt(2*np.log(2)))
      

        #CREATE THE IMAGE FILE
        w = wcs.WCS(naxis=2)
        #print original_header
        w.wcs.crpix = [original_header['CRPIX1'],original_header['CRPIX2']]
        w.wcs.cdelt = np.array([-pixsize/3600,pixsize/3600])
        w.wcs.crval = [original_header['CRVAL1'],original_header['CRVAL2']]
        w.wcs.ctype = [original_header['CTYPE1'],original_header['CTYPE2']]
        header_w = w.to_header()
        hdu = fits.PrimaryHDU(data = new_image_data,header = header_w)
        file_path = self.new_map_location + '/FAKE_MAP_NO_SOURCES.FITS'
        if os.path.exists(file_path) == True:
            os.remove(file_path)
            hdu.writeto(file_path)
        else:
            hdu.writeto(file_path)


map_path_original = '/home/sam/Desktop/LSB/data/maps/HATLAS_NGP_DR2_RAW250_CUTOUT.FITS'
noise_path_original = '/home/sam/Desktop/LSB/data/maps/HATLAS_NGP_DR2_SIGMA250_CUTOUT.FITS'
param_file_path = '/home/sam/Desktop/LSB/simstack-master/Parameters1.cfg'
params = parameters.get_params(param_file_path)


    
#map_basis = false_map(map_path_original,noise_path_original,params)
#false_map.create_fake_matrix()
            
        
