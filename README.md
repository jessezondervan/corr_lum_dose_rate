# corr_lum_dose_rate
Calculate the dose rate of luminescence samples, and correct for changes in the water table. 
The model assumes saturation water content until the water table drops below the samples at an age specified by the user,
after which the wet water content measured in the lab is used to correct dose rates.

# Requirements
This code is written in Python 3

# Usage
1. Calculate the external gamma dose rate added to the samples by cosmic rays. This depends on the altitude and geographical 
location of the samples, as well as the depth of samples. This can be done in a standard luminescence spreadsheet, for example
the one provided by Riso
2. Complete a csv file with sample name, a priori age (optional), end of saturation age (time of water table lowering below the sample),
saturation water content, wet water content, dry gamma dose rate, dry beta dose rate, conglomerate correction (if applicable), 
external gamma dose rate and the dose. See an example csv file, which you can download and use to fill in your own data
3. run the code for quartz and/or feldspar luminescence ages

The output is a print out of the corrected luminescence ages and dose rates based on the water content model, and also gives the
precision only error (leaving out systematic errors related to measurement uncertainties shared by all samples)

[![DOI](https://zenodo.org/badge/193535834.svg)](https://zenodo.org/badge/latestdoi/193535834)
