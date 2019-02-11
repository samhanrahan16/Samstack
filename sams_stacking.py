import os
import os.path
import sys
import time
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import ascii
from astropy.table import Table
from utils import *
from sams_utils import *
import sams_parameters
from sams_skymaps import sams_skymap


class samstack:
    '''Class contains methods to do some simple stacking, to plot the results, get images of the stacked maps'''

    def __init__(self,param_file_path):

        params = sams_parameters.get_params(param_file_path)
        self.zkey = params['zkey']
        self.mkey = params['mkey']
        self.rkey = params['ra_key']
        self.dkey = params['dec_key']

        self.m_nodes = params['bins']['m_nodes']
        self.z_nodes = params['bins']['z_nodes']
        self.cat_name = params['catalogs']['cat_name']

        self.n_slices = len(params['bins']['z_nodes']) - 1

        self.sky_library = get_maps(params)
        cats = get_catalogs(params)
        self.binned_ra_dec = get_bin_radec(params,cats)
        self.bin_names_to_midpoints = get_bin_names_to_midpoints(params,cats)
        self.bin_ids = get_bin_ids(params,cats)

    def simple_stack(self):
        t0 = time.time()

        map_library = self.sky_library
        subcatalog_library = self.binned_ra_dec

        map_names = [i for i in map_library.keys()]
        self.map_names = map_names
        cwavelengths = [map_library[i].wavelength for i in map_names]
        self.wavelengths = cwavelengths
        uwavelengths = np.sort(np.unique(cwavelengths))
        nwv = len(uwavelengths)

        bin_names = subcatalog_library.keys()
        self.bin_names = bin_names
        nbins = len(bin_names)

        stacks_all_wvs = {}
        radius = 1.0

        for iwv in range(nwv):
            print 'stacking '+ self.cat_name + ' ' + map_library.keys()[iwv]

            map_key = map_library.keys()[iwv]
            cmap = map_library[map_key].cmap
            cnoise = map_library[map_key].noise
            cwv = map_library[map_key].wavelength
            crms = map_library[map_key].rms
            chd = map_library[map_key].header
            pixsize = map_library[map_key].pixel_size
            fwhm = map_library[map_key].fwhm
            cw = WCS(chd)
            map_shape = np.shape(cmap)
            kern = gauss_kern(fwhm,np.floor(fwhm*10.)/pixsize,pixsize)

            # Create map with number of galaxies in each pixel

            ngals_per_pix_all = {}
            ngals_layer = {}
            stacked_fluxes = {}
            stacked_noises = {}
            stacked_stats = {}
            for name in bin_names:
                if len(subcatalog_library[name][0]) >0:
                    ngals_per_pix = np.zeros([map_shape[0],map_shape[1]])
                    ra = subcatalog_library[name][0]
                    dec = subcatalog_library[name][1]
                    ty,tx = cw.wcs_world2pix(ra,dec,0)
                    ind_keep = np.where((np.round(tx) >= 0) & (np.round(tx) < map_shape[0]) & (np.round(ty) >=0) & (np.round(ty) < map_shape[1]))
                    real_x = np.round(tx[ind_keep]).astype(int)
                    real_y = np.round(ty[ind_keep]).astype(int)

                    ind_nz = np.where(cmap[real_x,real_y] != 0)
                    nt = np.shape(ind_nz)[1]

                    ngals_layer[name] = nt
                    if nt >0:
                        real_x = real_x[ind_nz]
                        real_y = real_y[ind_nz]

                        for ni in range(nt):
                            ngals_per_pix[real_x[ni],real_y[ni]] +=1.0

                    ngals_per_pix_all[name] = ngals_per_pix
                else:
                    ngals_layer[name] = 1
                    ngals_per_pix_all[name] = 'empty'

            for name in bin_names:
                if ngals_per_pix_all[name] == 'empty':
                    stacked_fluxes[name] = 0
                    stacked_noises[name] = 0
                    stacked_stats[name] = {}
                    stacked_stats[name]['mean'] = 0
                    stacked_stats[name]['median'] = 0
                    stacked_stats[name]['stdev'] = 0
                else:
                
                    background = np.mean(np.ndarray.flatten(cmap))
                    gals_loc_map = ngals_per_pix_all[name]
                    ind_gals = np.where(ngals_per_pix_all[name] >=1)
                    gals_loc_map_around,ind_around = circle_around_galaxy(gals_loc_map,fwhm,pixsize,ind_gals)
                    gals_loc_map_convovled = smooth_psf(gals_loc_map,kern)
                    N_alpha = np.ndarray.flatten(gals_loc_map_convovled[ind_around])
                    fluxes = np.ndarray.flatten(cmap[ind_around])
                    fluxes_backsub = [flux - background for flux in fluxes]
                    fluxes_stats = np.ndarray.flatten(cmap[ind_gals])
                    fluxes_stats_backsub = [flux - background for flux in fluxes_stats]
                    noises = np.ndarray.flatten(cnoise[ind_around])
                    noises_backsub = [noise - background for noise in noises]
                    vals_to_sum = [N * flux for N,flux in zip(N_alpha,fluxes)]
                    noises_to_sum = [N*noise for N,noise in zip(N_alpha,noises)]
                    summed_vals = np.sum(vals_to_sum)
                    summed_noises = np.sum(noises_to_sum)
                    mean_gals_per_pix = np.mean(list(np.ndarray.flatten(gals_loc_map_around)))
                    N_pix = np.float(map_shape[0]) * np.float(map_shape[1])

                    stacked_flux = (1/(N_pix * mean_gals_per_pix))* summed_vals
                    stacked_noise = (1/(N_pix * mean_gals_per_pix)) * summed_noises
                    mean,median,var,stdev,summed = get_stacked_statistics(fluxes_stats)
                    stacked_fluxes[name] = stacked_flux
                    stacked_noises[name] = stacked_noise
                    stacked_stats[name] = {}
                    stacked_stats[name]['mean'] = mean
                    stacked_stats[name]['median'] = median
                    stacked_stats[name]['stdev'] = stdev
                names_to_mids = self.bin_names_to_midpoints

            packed_stack_outputs = pack_stacked_fluxes(stacked_fluxes,stacked_noises,stacked_stats,bin_names,ngals_layer,crms,names_to_mids)

            stacks_all_wvs[map_key] = packed_stack_outputs

        #print stacks_all_wvs
        t1 = time.time()
        tpass = t1 - t0

        logging.info('Done!')
        logging.info('')
        logging.info('Total time :{:.4f} minutes\n'.format(tpass/60.))

        return stacks_all_wvs

    def save_stacked_fluxes(self,save_path):
        '''Saves the output of the stacked fluxes in a csv format'''

        stack_output = self.simple_stack()
        tbl_maps = []
        bin_names_all = []
        stacked_fluxes_all = []
        stacked_noises_all = []
        means_all = []
        medians_all = []
        stdevs_all = []
        ngals_bin_all = []
        for map_name in self.map_names:
            mid_vals = []
            stacked_fluxes = []
            stacked_noises = []
            means = []
            medians = []
            stdevs = []
            ngals_bin = []
            bin_names = []
            
            for name in self.bin_names:
                mid_vals.append(stack_output[map_name][name]['bin_midpoint'])
                tbl_maps.append(map_name)
                bin_names.append(name)
                ngals_bin.append(stack_output[map_name][name]['ngals_bin'])
                stacked_fluxes.append(stack_output[map_name][name]['value'] * (10**3))
                stacked_noises.append(stack_output[map_name][name]['noise']*(10**3))
                means.append(stack_output[map_name][name]['mean']*(10**3))
                medians.append(stack_output[map_name][name]['median']*(10**3))
                stdevs.append(stack_output[map_name][name]['stdev']*(10**3))

            mid_vals_sorted = sorted(mid_vals)
            for val in mid_vals_sorted:
                i = mid_vals.index(val)
                bin_names_all.append(bin_names[i])
                stacked_fluxes_all.append(stacked_fluxes[i])
                stacked_noises_all.append(stacked_noises[i])
                means_all.append(means[i])
                medians_all.append(medians[i])
                stdevs_all.append(stdevs[i])
                ngals_bin_all.append(ngals_bin[i])
                

        

        results_table = Table([tbl_maps,bin_names_all,ngals_bin_all,stacked_fluxes_all,stacked_noises_all,means_all,medians_all,stdevs_all],names=['map','bin_name','ngals_bin','stacked_flux/mJy','stacked_noise/mJy','mean/mJy','median/mJy','stdev/mJy'])

        tbl_name = self.cat_name + '_stacked_output_' + 'bins_' + str(len(self.bin_names)) + '.csv'
        file_path = save_path + '/' + tbl_name

        if os.path.exists(save_path) == True:
            if os.path.exists(file_path) == True:
                os.remove(file_path)
                ascii.write(results_table,file_path,format='csv')
            else:
                ascii.write(results_table,file_path,format='csv')
        else:
            os.makedirs(save_path)
            ascii.write(results_table,file_path,format='csv')

    def plot_the_stacked_fluxes(self,save_path):
        '''This method plots the stacked fluxes for the different wavelengths and bins'''

        stack_output = self.simple_stack()
        plot_masses = []

        for map_name in self.map_names:
            mid_vals = []
            stacked_fluxes = []
            for name in self.bin_names:
                mid_val = stack_output[map_name][name]['bin_midpoint']
                stacked_flux = stack_output[map_name][name]['value'] * (10**3)
                mid_vals.append(mid_val)
                stacked_fluxes.append(stacked_flux)

            mid_vals_sorted = sorted(mid_vals)
            fluxes_sorted = []
            for val in mid_vals_sorted:
                index = mid_vals.index(val)
                flux = stacked_fluxes[index]
                fluxes_sorted.append(flux)

            plt.figure()
            #plt.plot(mid_vals_sorted,stacked_fluxes_sorted)
            plt.scatter(mid_vals,stacked_fluxes)
            savename = self.cat_name + '_' + map_name
            plt.title(savename)
            plt.grid()
            plt.xlabel('PetroMag_i')
            plt.ylabel('Flux/mJy')
            if os.path.exists(save_path) == True:
                plt.savefig(save_path + '/' + savename)
            else:
                os.makedirs(save_path)
                plt.savefig(save_path + '/' + savename)
            plt.show()

    def get_image_of_stack(self,side_length,table_path):
        '''Stacks the pixels and produces an image of each stacked output'''
        t0 = time.time()
        map_library = self.sky_library
        subcatalog_library = self.binned_ra_dec

        map_names = [i for i in map_library.keys()]
        cwavelengths = [map_library[i].wavelength for i in map_names]
        uwavelengths = np.sort(np.unique(cwavelengths))
        nwv = len(uwavelengths)

        bin_names = subcatalog_library.keys()
        nbins = len(bin_names)

        stacks_all_wvs = {}
        plot_matrices_all_wvs = {}
        plot_matrices_all_wvs_noise = {}

        for iwv in range(nwv):
            print 'stacking ' + self.cat_name + ' ' + map_library.keys()[iwv]

            map_key = map_library.keys()[iwv]
            cmap = map_library[map_key].cmap
            cnoise = map_library[map_key].noise
            cwv = map_library[map_key].wavelength
            crms = map_library[map_key].rms
            chd = map_library[map_key].header
            pixsize = map_library[map_key].pixel_size
            fwhm = map_library[map_key].fwhm
            cw = WCS(chd)
            map_shape = np.shape(cmap)

            shifted_list = create_shifted_list(side_length,pixsize)
            stacked_fluxes_all_shifts = {}
            stacked_noises_all_shifts = {}
            stacked_means_all_shifts = {}
            for ra_shift,dec_shift in shifted_list:
                binned_fluxes = {}
                binned_noises = {}
                binned_means = {}
                ngals_per_pix_all = {}
                ngals_layer = {}
                stacked_fluxes = {}
                stacked_noises = {}
                stacked_stats = {}
                for name in bin_names:
                    if len(subcatalog_library[name][0]) > 0:
                        ngals_per_pix = np.zeros([map_shape[0],map_shape[1]])
                        ra_shift_deg = ra_shift/3600.
                        ra = subcatalog_library[name][0]
                        ra_shifted = shift_RAs(ra,ra_shift_deg)
                        dec_shift_deg = dec_shift/3600.
                        dec = subcatalog_library[name][1]
                        dec_shifted = shift_DECs(dec,dec_shift_deg)
                        ty,tx = cw.wcs_world2pix(ra_shifted,dec_shifted,0)
                        ind_keep = np.where((np.round(tx) >=0) & (np.round(tx) < map_shape[0]) & (np.round(ty) >= 0) & (np.round(ty) < map_shape[1]))
                        real_x = np.round(tx[ind_keep]).astype(int)
                        real_y = np.round(ty[ind_keep]).astype(int)

                        ind_nz = np.where(cmap[real_x,real_y] != 0)
                        nt = np.shape(ind_nz)[1]

                        ngals_layer[name] = nt
                        if nt > 0:
                            real_x = real_x[ind_nz]
                            real_y = real_y[ind_nz]

                            for ni in range(nt):
                                ngals_per_pix[real_x[ni],real_y[ni]] += 1.0
                            ngals_per_pix_all[name] = ngals_per_pix
                        else: ngals_layer[name] = 1

                for name in bin_names:
                    background = np.mean(np.ndarray.flatten(cmap))
                    gals_loc_map = ngals_per_pix_all[name]
                    ind_gals = np.where(ngals_per_pix_all[name] >=1)
                    N_alpha = np.ndarray.flatten(gals_loc_map[ind_gals])
                    fluxes = np.ndarray.flatten(cmap[ind_gals])
                    fluxes_backsub = [flux - background for flux in fluxes]
                    noises = np.ndarray.flatten(cnoise[ind_gals])
                    noises_backsub = [noise - background for noise in noises]
                    vals_to_sum = [N * flux for N,flux in zip(N_alpha,fluxes)]
                    noises_to_sum = [N * noise for N,noise in zip(N_alpha,noises)]
                    summed_vals = np.sum(vals_to_sum)
                    summed_noise = np.sum(noises_to_sum)

                    mean_gals_per_pix = np.mean(list(np.ndarray.flatten(gals_loc_map)))
                    N_pix = np.float(map_shape[0]) * np.float(map_shape[1])
                    stacked_flux = (1/(N_pix * mean_gals_per_pix)) * summed_vals
                    stacked_noise = (1/(N_pix * mean_gals_per_pix)) * summed_noise
                    mean,median,var,stdev,summed = get_stacked_statistics(fluxes)
                    
                    stacked_fluxes[name] = stacked_flux
                    stacked_noises[name] = stacked_noise
                    stacked_stats[name] = {}
                    stacked_stats[name]['mean'] = mean
                    stacked_stats[name]['median'] = median
                    stacked_stats[name]['stdev'] = stdev
                    names_to_mids = self.bin_names_to_midpoints

                packed_stack_outputs = pack_stacked_fluxes(stacked_fluxes,stacked_noises,stacked_stats,bin_names,ngals_layer,crms,names_to_mids)
                for name in bin_names:
                    binned_fluxes[name] = packed_stack_outputs[name]['value']
                    binned_noises[name] = packed_stack_outputs[name]['noise']
                    binned_means[name] = packed_stack_outputs[name]['mean']
                stacked_fluxes_all_shifts[(ra_shift,dec_shift)] = binned_fluxes
                stacked_noises_all_shifts[(ra_shift,dec_shift)] = binned_noises
                stacked_means_all_shifts[(ra_shift,dec_shift)] = binned_means

            #CREATE TABLE OF VALS FOR PLOT
            for name in bin_names:
                ras = [ra_shift for ra_shift,dec_shift in shifted_list]
                decs = [dec_shift for ra_shift,dec_shift in shifted_list]
                fs = [stacked_fluxes_all_shifts[(ra_shift,dec_shift)][name] for ra_shift,dec_shift in shifted_list]
                ns = [stacked_noises_all_shifts[(ra_shift,dec_shift)][name] for ra_shift,dec_shift in shifted_list]
                ms = [stacked_means_all_shifts[(ra_shift,dec_shift)][name] for ra_shift,dec_shift in shifted_list]

                im_tbl = Table([ras,decs,fs,ns,ms],names=['ra_shift','dec_shift','flux','noise','mean'])
                tbl_name = self.cat_name + ':_' + map_key + '_' + name + '.csv'
                file_path = table_path + '/' + tbl_name
                if os.path.exists(table_path) == True:
                    if os.path.exists(file_path) == True:
                        os.remove(file_path)
                        ascii.write(im_tbl,file_path)
                    else:
                        ascii.write(im_tbl,file_path)
                else:
                    os.makedirs(table_path)
                    if os.path.exists(file_path) == True:
                        os.remove(file_path)
                        ascii.write(im_tbl,file_path)
                    else:
                        ascii.write(im_tbl,file_path)
            

            plot_matrices = {}
            plot_matrices_noise = {}
            for name in bin_names:
                values = [stacked_fluxes_all_shifts[(ra_shift,dec_shift)][name] for ra_shift,dec_shift in shifted_list]
                side_length = int(np.sqrt(float(len(values))))
                plot_matrix = np.zeros((side_length,side_length),dtype=float)
                plot_matrix_noise = np.zeros((side_length,side_length),dtype=float)
                matrix_indexes = find_image_indexes(shifted_list,pixsize,side_length)
                for ra_shift,dec_shift in shifted_list:
                    ra_index,dec_index = matrix_indexes[(ra_shift,dec_shift)]
                    plot_matrix[ra_index,dec_index] += stacked_fluxes_all_shifts[(ra_shift,dec_shift)][name]
                    plot_matrix_noise[ra_index,dec_index] += stacked_noises_all_shifts[(ra_shift,dec_shift)][name]
                plot_matrices[name] = plot_matrix
                plot_matrices_noise[name] = plot_matrix_noise

            plot_matrices_all_wvs[map_key] = plot_matrices
            plot_matrices_all_wvs_noise[map_key] = plot_matrices_noise

        t1 = time.time()
        tpass = t1 - t0

        logging.info('Done!')
        logging.info('')
        logging.info('Total time :{:.4f} minutes\n'.format(tpass/60.))

        return plot_matrices_all_wvs,plot_matrices_all_wvs_noise


    def plot_the_stacked_images(self,side_length,save_path):
        '''Plots the stacked images'''

        subcatalog_library = self.binned_ra_dec
        map_library = self.sky_library
        bin_names = subcatalog_library.keys()
        map_names = map_library.keys()
        plot_matrices_all_wvs,plot_matrices_all_wvs_noise = self.get_image_of_stack(side_length,save_path)

        if os.path.exists(save_path) == True:
            for wv in map_names:
                for name in bin_names:
                    image_name = self.cat_name + ':_' + str(wv) + '_' + str(name)
                    save_name = image_name + '.jpg'
                    file_path = save_path + '/' + save_name
                    plot_matrix = plot_matrices_all_wvs[wv][name]
                    plt.figure()
                    plt.imshow(plot_matrix)
                    plt.colorbar()
                    plt.title(image_name)
                    plt.savefig(file_path)
                    plt.show()
                    
        else:
            os.makedirs(save_path)
            for wv in map_names:
                for name in bin_names:
                    image_name = self.cat_name + ':_' + str(wv) + '_' + str(name) 
                    save_name = image_name + '.jpg'
                    file_path = save_path + '/' + save_name
                    plot_matrix = plot_matrices_all_wvs[wv][name]
                    plt.figure()
                    plt.imshow(plot_matrix)
                    plt.colorbar()
                    plt.title(image_name)
                    plt.savefig(file_path)
                    plt.show()

                    

        

    def plot_stacked_flux_increasing_gals(self,save_path):
        '''Plots a graph of the stacked flux increasing the number of galaxies'''
        t0 = time.time()

        map_library = self.sky_library
        subcatalog_library = self.binned_ra_dec

        map_names = [i for i in map_library.keys()]
        cwavelengths = [map_library[i].wavelength for i in map_names]
        uwavelengths = np.sort(np.unique(cwavelengths))
        nwv = len(uwavelengths)

        bin_names = subcatalog_library.keys()
        nbins = len(bin_names)

        stacks_all_wvs = {}

        for iwv in range(nwv):
            print 'stacking ' + self.cat_name + '_' + map_library.keys()[iwv]

            map_key = map_library.keys()[iwv]
            cmap = map_library[map_key].cmap
            cnoise = map_library[map_key].noise
            cwv = map_library[map_key].wavelength
            crms = map_library[map_key].rms
            chd = map_library[map_key].header
            pixsize = map_library[map_key].pixel_size
            fwhm = map_library[map_key].fwhm
            cw = WCS(chd)
            map_shape = np.shape(cmap)
            kern = gauss_kern(fwhm,np.floor(fwhm*10.)/pixsize,pixsize)

            ngals_per_pix_all = {}
            ngals_layer = {}
            stacked_fluxes = {}
            for name in bin_names:
                ngals_diff_lengths = {}
                length_list = range(10,len(subcatalog_library[name][0])+1,4)
                if len(subcatalog_library[name][0]) > 0:
                    for length in length_list:
                        ngals_per_pix = np.zeros([map_shape[0],map_shape[1]])
                        ra = subcatalog_library[name][0][:length]
                        dec = subcatalog_library[name][1][:length]
                        ty,tx = cw.wcs_world2pix(ra,dec,0)

                        ind_keep = np.where((np.round(tx) >= 0) & (np.round(tx) < map_shape[0]) & (np.round(ty) >= 0) & (np.round(ty) < map_shape[1]))

                        real_x = np.round(tx[ind_keep]).astype(int)
                        real_y = np.round(ty[ind_keep]).astype(int)

                        ind_nz = np.where(cmap[real_x,real_y] != 0)
                        nt = np.shape(ind_nz)[1]

                        ngals_layer[name] = nt
                        if nt >0:
                            real_x = real_x[ind_nz]
                            real_y = real_y[ind_nz]

                            for  ni in range(nt):
                                ngals_per_pix[real_x[ni],real_y[ni]] += 1.0

                        ngals_diff_lengths[length] = ngals_per_pix
                    ngals_per_pix_all[name] = ngals_diff_lengths
                else: ngals_layer[name] = 1

            for name in bin_names:
                length_list = range(10,len(subcatalog_library[name][0])+1,4)
                stacked_flux_diff_lengths = {}
                for length in length_list:
                    gals_loc_map = ngals_per_pix_all[name][length]
                    ind_gals = np.where(gals_loc_map >= 1)
                    gals_loc_map_around,ind_around = circle_around_galaxy(gals_loc_map,fwhm,pixsize,ind_gals)
                    gals_loc_map_convolved = smooth_psf(gals_loc_map,kern)
                    N_alpha = np.ndarray.flatten(gals_loc_map_convolved[ind_around])
                    fluxes = np.ndarray.flatten(cmap[ind_around])
                    noises = np.ndarray.flatten(cnoise[ind_around])
                    fluxes_stats = np.ndarray.flatten(cmap[ind_gals])
                    vals_to_sum = [N * flux for N,flux in zip(N_alpha,fluxes)]
                    noises_to_sum = [N*noise for N,noise in zip(N_alpha,noises)]
                    summed_vals = np.sum(vals_to_sum)

                    mean_gals_per_pix = np.mean(list(np.ndarray.flatten(gals_loc_map_around)))
                    N_pix = np.float(map_shape[0]) * np.float(map_shape[1])
                    stacked_flux = (1/(N_pix * mean_gals_per_pix)) * summed_vals
                    stacked_flux_diff_lengths[length] = stacked_flux
                stacked_fluxes[name] = stacked_flux_diff_lengths

        for name in bin_names:
            length_list = range(10,len(subcatalog_library[name][0]) +1,4)
            stacked_flux_data = [stacked_fluxes[name][length] for length in length_list]
            stacked_flux_mjy = [flux*(10**3) for flux in stacked_flux_data]
            save_name = self.cat_name + ':_' + map_key + '_' + name
            file_path = save_path + '/' + save_name
            plt.figure()
            plt.plot(length_list,stacked_flux_mjy)
            plt.grid()
            plt.xlabel('N')
            plt.ylabel('stacked_flux/mJy')
            plt.savefig(file_path)
            plt.show()

        t1 = time.time()
        tpass = t1 - t0

        logging.info('Done!')
        logging.info('')
        logging.info('Total time :{:.4f} minutes\n'.format(tpass/60.))
            
                                
                    
                    
                    
                        
                        
                


    

            
                    
                        

        


        
            

        

   

        

        

            
