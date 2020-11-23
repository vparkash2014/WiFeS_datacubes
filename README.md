
# Implementing pPXF on WiFeS DataCubes
Hello fellow WiFeS User/ student. Congrats on a successful observing run! Now time to have some fun. If you are trying to extract the spectra from WiFeS IFU and fit a stellar spectrum then you are in the right spot.

Below provide a description of the main fitting script and all the other files and functions included in this "package." If you have any questions, you can contact me at vparkash2014@gmail.com. Remember NO question is a dumb question.

Main script:
The main function of this script is to extract stellar kinematic, gas kinematic and emission line fluxes from WiFeS IFU cubes using pPXF. This script will only work on reduced WiFeS cubes. If you have not already done so, reduce your observations using [PyWiFeS](http://www.mso.anu.edu.au/pywifes/doku.php).

Warning: this code was designed to specifically for analysing WiFeS cubes of local galaxies (z < 0.01). It can be adopted for different types of observations but with some caution. Please contact me if you need assistance.

More information about [WiFeS](https://rsaa.anu.edu.au/observatories/instruments/wide-field-spectrograph-wifes) and [pPXF](https://www-astro.physics.ox.ac.uk/~mxc/software/) can be found their respected webpages.

Below is a list of required the packages that you will need:
* Python 3.8
* pPXF : This script has been tested on the version of pPXF released on 11 September 2020. Please contact me if newer versions of pPXF breaks this script and I will de-bug it.
* Numpy 1.18
* Scipy 1.4
* Matplotlib 3.1
* Astropy
* Seaborn
* [VorBin](https://pypi.org/project/vorbin/)


## Summary of the script/ methods:
1. The script implements the Voronoi binning method (Cappellari & Copin,2003, MNRAS, 342, 345) "optimally solves the problem of preserving the maximum spatial resolution of general two-dimensional data (or higher dimensions), given a constraint on the minimum signal-to-noise ratio."  You will need to set the target signal-to-noise ratio ("targetSN"). To figure out the ideal "targetSN" for your data, I suggest that you measure the signal to noise ratio of your observation and determine the perfect balance between resolution and S/N that you need for your science. I have written a separate jupyter notebook that only runs the Voronoi bin method to test a range of S/N (Please contact me if you are interested in that script). Note that the script uses the red arm (data cube) to determine the bins because of my science goals. If you don't want to run the Voronoi bin method, set "targetSN" to a very small number so each spaxel is in its bin.

2. The main part of this script it to implments pPXF on the Voronio bin spectra to extract kinematics and gas emission lines. Using the function "retrieveTemplates", the script will look for [INDO-US_Valdes stellar spectra](https://www.noao.edu/cflib/). I have included the templates. You will need to change the path (line 206 in WiFeS_function.py) to these templates in the "retrieveTemplates" function for this to work.

3. Once we have the kinematics for each bin (and hence pixel), the script will velocity-shift all the spaxels to the centre of the galaxy and combine all the spaxels to create a single integrated galaxy spectrum. Then, it will re-run pPXF on the integrated spectrum. Note, I roughly implement the same method when I extract a nuclear spectrum for in a 3" diameter aperture. That script is currently not included in this script, but if you are interested, please contact me.


## Syntax:
```console
$ python Main_script.py input_file
```

## Input:
The only input variable is a text file that lists/ defines the variables to run this script. Below is the list of variables that you should define in this text file. Each line should have one variable (no need to include the name of the variable-- i.e. "galaxy = ". The variables should be ordered in the following way. Please look at example input file for clarification.

Input file:
```console
destdir (the directory where the datacubes are)
fileR (The R arm datacube)
fileB (The B arm datacube)
vel (Velocity of the galaxy in km/s. I usually use the values reported from NED)
PA
ba
targetSN (The target signal-to-noise ratio for Voronoi bin)
```

## Output:
The output of this script includes various of plots that you can look at to see if the script is working the way you want it to and fits files containing information of the kinematic and emission lines:

destdir+'/'+hdrB['OBJECT']+'_voronoi_2d_binning_pp.csv' -- csv that save the results of the Voronoi bin method. Column: x pixel, y pixel, signal, noise, binNum

destdir+'/'+hdrB['OBJECT']+'.pdf' -- Shows results from the Voronoi binning. The left and right panels shows the results for the red and blue arm respectively. The 1st row is the 2D image of the datacubes, the 2nd row shows the location and size of each bin, the 3rd row shows the 2D image of the datacubes after the spaxels are combined with respect to the Voronio bins and 4th row is the difference between the 1st and 3rd row. The main purpose of these plots is to see how the spatial resolution is affected by the Voronoi bin method.

destdir+'/'+hdrB['OBJECT']+'_voronoi_spectra.pdf' -- This pdf/plots compares 10 random individual spaxel spectrum to the Voronoi bin spectrum. This is a visual test to see if we are obtaining the S/N that we want and to see if the Voronoi bin spectrum are combined correctly.

destdir+'/bins/'+'pPXFoutput_'+str(bin)+'.csv' -- The output csv file of pPXF for each bin. Columns: wavelength, observed spectrum, bestfit spectrum, stellar component, & gas component.

destdir+'/bins/'+'pPXFoutput_'+str(bin)+'.pdf' -- Plots of the spectra, the bestfit model, the stellar component & the gas component for each bin. I also include zoom-ins to the key emission line features to visual inspect how well the emission-lines were fitted.

destdir+'/bins/'+hdrB['OBJECT']+'_ppxf_velgas_output.csv' -- csv file that stores the kinematic and emission fluxes for all the bins. Columns: bin number, gas velocity, stellar velocity, gas sigma, stellar sigma, emission flux of H delta, H gamma, H beta, H alpha, SII 6717, SII 6731, OIII 5007d, OI 6300d, NII

destdir+'/'+hdrB['OBJECT']+'_stellar_vel.fits' -- fits file that assigns a stellar velocity each spaxel according the results of pPXF

destdir+'/pPXFoutput_integral_spec.csv' --  The output csv file of pPXF for the integrated spectra. Columns: wavelength, observed spectrum, bestfit spectrum, stellar component, & gas component.

destdir+'/pPXFoutput_integral_spec.pdf' -- Plots of the spectra, the bestfit model, the stellar component & the gas component for the intergrated spectra.

destdir+'/integral_spec_ppxf_velgas_output.csv' -- csv file that stores the kinematic and emission fluxes for the integrated spectra. Columns: bin number, gas velocity, stellar velocity, gas sigma, stellar sigma, emission flux of H delta, H gamma, H beta, H alpha, SII 6717, SII 6731, OIII 5007d, OI 6300d, NII


## Important Note:
* The cosmology applied in this script is H0 =70kms, ΩM =0.3, and ΩΛ =0.7. This can be fixed on line ## of the Main_script.py.
* Please cite the proper papers to [Voronoi bins](https://ui.adsabs.harvard.edu/abs/2003MNRAS.342..345C/abstract), pPXF ([Paper 1](https://ui.adsabs.harvard.edu/abs/2017MNRAS.466..798C/abstract), [Paper 2](https://ui.adsabs.harvard.edu/abs/2004PASP..116..138C/abstract)), WiFeS ([Paper 1](https://ui.adsabs.harvard.edu/abs/2007Ap%26SS.310..255D/abstract), [Paper 2](https://ui.adsabs.harvard.edu/abs/2010Ap%26SS.327..245D/abstract)) and [PyWiFeS](https://ui.adsabs.harvard.edu/abs/2014Ap%26SS.349..617C/abstract) in your publication. 
* Main_script.py calls upon the functions written in WiFeS_function.py. Make sure that this file is in the same directory. To ensure that everything is running properly, I add the WiFeS_function.py to my Python path

## Example Code:
The 'Example' directory contains example input files, and some of the output files.
I was not able to uplaod the datacubes to GitHub as they are too large. However you can access the datacube and the rest of output files from my [dropbox](https://www.dropbox.com/sh/39fhc3dyh56r0lm/AACZxqE9pa97py_9q4gplPn3a?dl=0)


In this directory, there is a C-shell script called "runscript". This script will run the "Main_script.py." Just make sure all the paths are defined properly.
To execute, make sure you are in csh, tcsh or zcsh and type the following command into a terminal:
```console
$ ./runscript
```

## Extra code: 
I have written various of visualisation scripts that are currently not available at this time. Please check-out my paper ([Parkash et al., 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.3169P/abstract)) to see those visualistions/ plots. If you like what you see and interested in those scripts, please free to contact me. 
