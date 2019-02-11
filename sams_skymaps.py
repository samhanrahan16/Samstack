from utils import gauss_kern
from utils import map_rms
from astropy.io import fits
import numpy as np



class sams_skymap:

    def __init__(self,map_file,map_noise,psf,colour_correction =1.0,wavelength=None,fwhm=None):

        cmap , header = fits.getdata(map_file, 0 ,header=True)

        cnoise, noise_header = fits.getdata(map_noise, 0,header=True)

        if header['CDELT2'] == 1.0:
            pixsize = header['PC2_2'] * 3600.
        else:
            pixsize = header['CDELT2'] * 3600

        #kern = gauss_kern(psf,np.floor(psf*8.)/pixsize,pixsize)

        self.cmap = cmap * colour_correction
        self.noise = cnoise * colour_correction
        self.header = header
        #self.psf = kern
        self.rms = map_rms(self.cmap.copy(),silent=True)
        self.pixel_size = pixsize
        if fwhm != None:
            self.fwhm = fwhm
        if wavelength != None:
            self.wavelength = wavelength

        
