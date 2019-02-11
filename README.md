# Samstack
A stacking code based on Marco Viero's simstack code. It is a simplifed and more user friendly version developed in Python.

The samstack class can be used to perform stacking of a number of galaxies within a catalog. The maps and catalogs should be specified 
within the parameter file. An example Parameter file Sams_paramters_example.cfg can be found. The catalog can be split into bins of its mass or
redshift and multiple maps can be stacked at once. It is important that only one catalog can be stacked at a time and this should be 
specified in the parameter file. See the run_samstack_example notebook as an example on how to run the outputs of samstack.

To run samsack, it needs to be initialised with the path to the parameter file. The methods then require the paths of the file where you
would like to save the stacked outputs.

The output will be saved in the form of a csv file, while the images will be saved in image format.

For analysis of the stacked results the file contains a SED fitter. This can be used to fit a modified blackbody function to your stacked 
fluxes for further analysis of them.



