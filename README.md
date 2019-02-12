# Samstack
A stacking code based on Marco Viero's simstack code. It is a simplifed and more user friendly version developed in Python.

The samstack class will run stacking of a single catalog on several different skymaps. It will produce files of the stacked fluxes 
and images of 'stacked galaxy'. It also contains a SED fitter for analysis of the stacked fluxes.

# Setup
Ensure the following repositories are installed:
* Python 2
* numpy
* scipy
* matplotlib
* astropy
* emcee
* corner

The catalogs need to be in a csv format and must contain columns with the RA,DEC,redshift,Mass of the galaxies (in our example we 
use i band magnitude as a substitute for mass).

The catalog names, catalog paths and file name of the catalogs should be specified in the parameter file, as should the map 
names,paths and file names. Follow the example file provided 'Sams_parameters_example.cfg' as a reference. 

The data can be separated into mass bins, by specyfing the bin edges in the parameter file.

The samstack class should be initialised with the path to the parameter file. The seperate methods can then be run to produce the
stacking output. Follow the example notebook 'Run_samstack_example' to see how it works, using the map and catalog in the 
'example_data' folder.

# Results

The stacked fluxes for each map and mass bin are saved into a csv file. While the stacked fluxes are plotted against the mid value
of their mass bin.

Images of the 'stacked galaxy' are also saved into the path provided.

# Additional Features

The modified_bb class in the 'SED_fitter' file can be used to fit a modified blackcckbody spectrum to the flux densities from the 
stacking, using MCMC from emcee.

In order to test the stacking is working correctly, a fake map where the fluxes are all known can be created using the 
'false_map' file. Use the 'create_fake_map' notebook as an example on how to do this.


# Future Work
The catalog can currently only be binned by mass and not redshift. We are working on including this.

We are currently working on adding a bootstrap technique to the stacking, to provide an error on the stacked flux.

# People

* Sam Hanrahan(Imperial College London)
* James Harrison(Imperial College London)

# Credits

This code is based upon the work of Marco Viero's 'simstack' see viero et. al. (2013).




