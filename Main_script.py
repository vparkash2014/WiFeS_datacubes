#!/usr/bin/env python
##############################################################################
#
# this script implements the voronoi bin method and pPXF method
# Please read the README file for more details
##############################################################################


# Loading all the needed libraries
from time import perf_counter as clock
from os import path

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import optimize

from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from astropy.io import fits
import astropy.units as u

from ppxf import ppxf
from vorbin import voronoi_2d_binning
from WiFeS_function import *


##### Testing input Parameters

try:
    input = sys.argv[1]
    print('The input file is ', input)
except:
    print("No input file was provide")
    sys.exit()

f = open(str(input), 'r')
Lines = f.readlines()

if len(Lines) != 7:
    print("There is not information in the input file.")
    sys.exit()

destdir = str(Lines[0].strip())
fileR, fileB = str(str(Lines[1].strip())), str(str(Lines[2].strip()))
vel = float(Lines[3].strip())
PA = float(Lines[4].strip())
ba = float(Lines[5].strip())
targetSN_R = float(Lines[6].strip())


##### Calculating redshift
c = 299792.458 #Km/s
z = vel/c


##### Read datacube

crop = [2,35,0,26] # since the edges of the cubes are messy, I like to crop the datacube. You can uncomment this if you want.

try:
    cubeR, hdrR, varR, lambR, lamRangeR = read_fits(destdir + "/" + fileR, limits = crop)
    cubeB, hdrB, varB, lambB, lamRangeB = read_fits(destdir + "/" + fileB, limits = crop)
except:
    print("Can't find the datacubes. Check the path")
    sys.exit()

print('Working with Galaxy: ', hdrR['OBJECT'])


##### Finding the centre galaxy (based of Sextractor)
center = centroid(cubeR, PA, ba)
print("The center of the galaxy is ", center[0][0], center[0][1] )



##### Voronoi binning the data by using the var information in the data cube as the noise
x_1d, y_1d = np.array([x for x in range( cubeR.shape[2] ) for y in range( cubeR.shape[1] )]),\
    np.array([y for x in range( cubeR.shape[2] ) for y in range( cubeR.shape[1] )])

# we don't use the full spectrum to determine the median of each spaxel because of emission lines. Below, I have defined the region that is flat and has no spectrum features
bin_red = np.array([6000.,6200.])*(1+z)
bin_blue = np.array([5310, 5680])*(1+z)

signal_R, noise_R = singal_noise_var_sec(cubeR, varR, lamRangeR, lambR, bin_red)
signal_B, noise_B = singal_noise_var_sec(cubeB, varB, lamRangeB, lambB, bin_blue)

try:
    print('Voronoi binning the R data')
    ind_neg_R, = np.where(signal_R <= 0)
    ind_pos_R, = np.where(signal_R > 0)
    binNum, xNode_R, yNode_R, xBar_R, yBar_R, sn_R, nPixels_R, scale_R\
        = voronoi_2d_binning(x_1d[np.where(signal_R > 0)], y_1d[np.where(signal_R > 0)], \
            signal_R[np.where(signal_R > 0)], noise_R[np.where(signal_R > 0)], targetSN_R,\
                             plot=False, quiet=True)

except ValueError:
    print("voronoi binning R won't work")
    print("")
    binNum, xNode_R, yNode_R, xBar_R, yBar_R, sn_R, nPixels_R, scale_R = \
    np.arange(len(x_1d)),np.zeros(len(x_1d)),np.zeros(len(x_1d)),np.zeros(len(x_1d)),np.zeros(len(x_1d)),np.zeros(len(x_1d)),np.zeros(len(x_1d)), np.zeros(len(x_1d))
    pass

binNum_R = np.zeros(len(x_1d))
binNum_R[ind_pos_R] = binNum
binNum_R[ind_neg_R] = -1
binNum_R = [int(i) for i in binNum_R]

filename = destdir+'/'+hdrB['OBJECT']+'_voronoi_2d_binning_pp.csv'
dat = Table([x_1d, y_1d, signal_R, noise_R, binNum_R], names=('x', 'y','signal', 'noise','bin'))
dat.write(filename, format='csv', overwrite=True)

fig = voronoi_outputs(cubeR, signal_R, noise_R, signal_B, noise_B, xNode_R, xNode_R, binNum_R, binNum_R,destdir+'/'+hdrB['OBJECT']+'.pdf')

fig = voronoi_spectra_outputs(cubeR, cubeB, lamRangeR, lambR, binNum_R, lamRangeB, lambB, binNum_R, destdir+'/'+hdrB['OBJECT']+'_voronoi_spectra.pdf')



##### Setup the input parameters 

spectrum_B = single_spectra(cubeB, results[0][0],results[0][1])
spectrum_R = single_spectra(cubeR, results[0][0],results[0][1])
flux_t, lamb_t = glue_BR(spectrum_B, lambB, spectrum_R, lambR)

specNew, logLam, velscale = util.log_rebin([lamb_t[0], lamb_t[-1]], flux_t, oversample=False,flux=False)
wave_limit = [3800,6950] 
FWHM_gal = np.mean(np.exp(logLam)/3000)
wave = np.exp(logLam)[np.where((np.exp(logLam) >= wave_limit[0]) & (np.exp(logLam) <= wave_limit[1]))]
velscale = c*np.log(wave[1]/wave[0])
lam_range_gal = np.array([np.min(wave), np.max(wave)])/(1 + z) #redshift corrected
print('the un-redshift wavelength is:', lam_range_gal)
print('the redshift wavelength is:', [wave[0], wave[-1]])

nagalaxy, logLam1, navelscale = util.log_rebin(lam_range_gal, flux_t, oversample=False, flux=False) # the rest frame velocity


##### Setup templates 


try:
    logLam2, stars_templates = retrieveTemplates(velscale, FWHM_gal, range_sel = [3400,7200])#, ID_template='INDO_US')
except:
    print "the stellar templates are not found. Please fix the path to the stellar templates in the retrieveTemplates function and try again."
    sys.exit(1)

print("")
gas_templates, line_names, line_wave = util.emission_lines(logLam2, lam_range_gal,FWHM_gal)
print('the gas lines are:', line_names)
templates = np.hstack([stars_templates, gas_templates])

dv = c*(logLam2[0] - np.log(wave[0])) #the difference between the templetes vel and the start of the redshift observed velocity...

start = [vel, 180.]  # (km/s), starting guess for [V, sigma]

nTemps = stars_templates.shape[1]
nLines = gas_templates.shape[1]
component = [0]*nTemps + [1]*nLines
moments = [4, 2]  # fit (V, sig, h3, h4) for the stars and (V, sig) for the gas
start = [start, start]  # adopt the same starting value for both gas and stars
logLam1 = np.log(wave/(1 + z)) #redshift corrected wavelengths
goodpixels = util.determine_goodpixels(logLam1, [np.exp(logLam2[0]), np.exp(logLam2[-1])], z )  #Generates a list of goodpixels to mask a given set of gas emissionlines.



##### ------------------- Running loop over all the bins -----------------------

x_1d, y_1d = np.array([x for x in range( cubeR.shape[2] ) for y in range( cubeR.shape[1] )]),\
                 np.array([y for x in range( cubeR.shape[2] ) for y in range( cubeR.shape[1] )])

spec_vor, solutions, velF, sigmaF, V_gas_vector, V_stellar_vector, gas_sigma, stellar_sigma = [], [], [], [], [], [], [], []
Hdelta, Hgamma, Hbeta, Halpha, SII_6717, SII_6731, OIII_5007d, OI6300d, NII, fluxes =  [], [], [], [], [], [], [], [], [], []
EW_Ha_vector, SFR_vector, Ha_area1, NII_area_22, SII_area1 = [], [], [], [], []

errvelF, errsigmaF, errh3F, errh4F, all_data = [], [], [], [], []

for bin in np.unique(binNum_R):
    x, y = x_1d[np.where(binNum_R == bin)], y_1d[np.where(binNum_R == bin)]
    print("pPXF is currently working on bin," bin,len(x))
    spectrum_R = np.sum(np.array([cubeR[:, int(ypair), int(xpair)] for ypair,xpair in zip(y,x)]), axis = 0)
    spectrum_B = np.sum(np.array([cubeB[:, int(ypair), int(xpair)] for ypair,xpair in zip(y,x)]), axis = 0)
    spectrum_B = spectrum_B[0:2830]
    lambB = lambB[0:2830]
    flux_t, lamb_t = glue_BR(spectrum_B, lambB, spectrum_R, lambR)
    spec_vor.append(flux_t)

    specNew, logLam, velscale = util.log_rebin([lamb_t[0], lamb_t[-1]], flux_t, oversample=False,flux=False)

    galaxy = specNew[np.where((np.exp(logLam) >= wave_limit[0]) & (np.exp(logLam) <= wave_limit[1]))]


    ##### pPXF: Here the actual fit starts. The best fit is plotted on the screen.

    noise = np.sqrt(np.var(galaxy))
    noise = np.full_like(galaxy,noise)

    fig = plt.figure(num=0, figsize=(12.6, 5))
    pp = ppxf(templates, galaxy, noise, velscale, start, plot=False, moments=moments, degree=-1,\
                  mdegree=10, vsyst=dv, clean=False, component=component, quiet=True)

    gas = np.array(component) == 1  # Select weights of gas emissions only
    stars = pp.matrix[:, ~gas].dot(pp.weights[~gas])
    gas = pp.matrix[:, gas].dot(pp.weights[gas])
    w = np.where(np.array(component) == 1)[0]

    filename = destdir+'/bins/'+'pPXFoutput_'+str(bin)+'.csv'
    dat = Table([wave, pp.galaxy, pp.bestfit, stars, gas], names=('wave', 'galaxy','bestfit', 'stars', 'gas'))
    dat.write(filename, format='csv', overwrite=True)

    fig, ax = plt.subplots(figsize = (17,15))
    plt.subplot(211)
    plt.xlabel("Observed Wavelength ($\AA$)")
    plt.ylabel("Flux")
    plt.xlim([np.min(wave), np.max(wave)])

    plt.plot(wave, pp.galaxy, 'k')
    plt.axhline(y=-0, linestyle='--', color='k', linewidth=2)
    plt.plot(wave, (stars), 'r', linewidth=2)  # overplot stellar templates alone
    plt.plot(wave, (pp.galaxy-stars),'o', ms=1, color='LimeGreen', mec='LimeGreen')
    plt.plot(wave, (gas), 'b', linewidth=1)  # overplot emission lines alone

    plt.savefig(destdir+'/bins/'+'pPXFoutput_'+str(bin)+'.pdf', bbox_inches='tight')
    plt.clf()

    
    ##### Defining all the outputs
    solutions.append(pp.sol)
    velF.append(pp.sol[0])
    sigmaF.append(pp.sol[1])
    velF.append(pp.sol[0])
    sigmaF.append(pp.sol[1])
    V_gas_vector = np.append(V_gas_vector,pp.sol[1][0])
    V_stellar_vector = np.append(V_stellar_vector,pp.sol[0][0])
    gas_sigma = np.append(gas_sigma,pp.sol[1][1])
    stellar_sigma = np.append(stellar_sigma,pp.sol[0][1])

    for name, weight, line in zip(line_names, pp.weights[w], pp.matrix[:,w].T):
        flux = weight*np.max(line)
        fluxes = np.append(fluxes, flux)

    Hdelta, Hgamma, Hbeta, Halpha, SII_6717, SII_6731, OIII_5007d, OI6300d, NII, fluxes =\
    np.append(Hdelta,fluxes[0]), np.append(Hgamma,fluxes[1]), np.append(Hbeta,fluxes[2]), np.append(Halpha,fluxes[3]),\
    np.append(SII_6717,fluxes[4]), np.append(SII_6731,fluxes[5]), np.append(OIII_5007d,fluxes[6]),\
    np.append(OI6300d,fluxes[7]), np.append(NII,fluxes[8]), []

print("pPXF has ran through all the bins- all done. You should check the results in the bin directory")
filename = destdir+'/bins/'+hdrB['OBJECT']+'_ppxf_velgas_output.csv'
dat = Table([np.unique(binNum_R),V_gas_vector,V_stellar_vector, gas_sigma, stellar_sigma, \
                            Hdelta, Hgamma, Hbeta, Halpha, SII_6717, SII_6731, OIII_5007d, OI6300d, NII], \
                names=("bin","V_gas_vector","V_stellar_vector", "gas_sigma", "stellar_sigma", \
                            "Hdelta", "Hgamma", "Hbeta", "Halpha","SII_6717", "SII_6731", "OIII_5007d", "OI6300d", "NII"))
dat.write(filename, format='csv', overwrite=True)




##### Each spaxels will be assigned a stellar velocity based of ppxf
print("Next we are going to create an integrated spectrum and run pPXF on that")
filename = destdir+'/'+hdrB['OBJECT']+'_voronoi_2d_binning_pp.csv'
bin_data = Table.read(filename, format = 'csv')

filename = destdir+'/bins/'+hdrB['OBJECT']+'_ppxf_velgas_output.csv'
fits_v = Table.read(filename, format = 'csv')
star_vel_cube = np.zeros( ( len(cubeR[0]), len(cubeR[0,0,:]) ) )

for ind, bin in enumerate(fits_v['bin']):
    coord = bin_data['x','y'][np.where(bin_data['bin'] == bin)]
    if bin >= 0.0:
        for spaxel in coord:
            x, y = np.array(spaxel[0]), np.array(spaxel[1])
            star_vel_cube[ y, x] = fits_v['V_stellar_vector'][ind]
hdu = fits.PrimaryHDU()
hdu.data = star_vel_cube
hdu.header
filename = destdir+'/'+hdrB['OBJECT']+'_stellar_vel.fits'
hdu.writeto(filename,  overwrite = True)

central_vel = star_vel_cube[int(results[0][1]), int(results[0][0])]



##### Each voronoi bin will be velocity-shifted to that of the central spaxel and an integrated spectra is created

fig, ax = plt.subplots(figsize = (17,7))
integral_spec = np.zeros(len(wave))

for ind, bin in enumerate(fits_v['bin']):
    if bin >= 0.0:
        filename = destdir+'/bins/'+'pPXFoutput_'+str(bin)+'.csv'
        fits_data = Table.read(filename, format = 'ascii.csv')
        bin_vel = fits_v['V_stellar_vector'][ind]

        diff = np.round( (central_vel - bin_vel) / velscale )
        shift_wave = np.round( np.log(wave) + ( diff*(velscale/c) ) , 10)
        for index,value in enumerate( np.round( np.log(wave), 10) ) :
            mask1 = (shift_wave == value)
            if np.sum(mask1) >= 1:
                integral_spec[index] += fits_data['galaxy'][mask1]


##### ------------------- Running ppxf on the integral spectrum -----------------------

spec_vor, solutions, velF, sigmaF, V_gas_vector, V_stellar_vector, gas_sigma, stellar_sigma = [], [], [], [], [], [], [], []
Hdelta, Hgamma, Hbeta, Halpha, SII_6717, SII_6731, OIII_5007d, OI6300d, NII, fluxes =  [], [], [], [], [], [], [], [], [], []
EW_Ha_vector, SFR_vector, Ha_area1, NII_area_22, SII_area1 = [], [], [], [], []

errvelF, errsigmaF, errh3F, errh4F, all_data = [], [], [], [], []



##### pPxf Here the actual fit starts. The best fit is plotted on the screen.

noise = np.sqrt(np.var(integral_spec))
noise = np.full_like(integral_spec,noise)

fig = plt.figure(num=0, figsize=(12.6, 5))
pp = ppxf(templates, integral_spec, noise, velscale, start, plot=False, moments=moments, degree=-1,\
                  mdegree=10, vsyst=dv, clean=False, component=component, quiet=True)

gas = np.array(component) == 1  # Select weights of gas emissions only
stars = pp.matrix[:, ~gas].dot(pp.weights[~gas])
gas = pp.matrix[:, gas].dot(pp.weights[gas])
w = np.where(np.array(component) == 1)[0]

filename = destdir+'/pPXFoutput_integral_spec.csv'
dat = Table([wave, pp.galaxy, pp.bestfit, stars, gas], names=('wave', 'galaxy','bestfit', 'stars', 'gas'))
dat.write(filename, format='csv', overwrite=True)

fig, ax = plt.subplots(figsize = (17,15))
plt.subplot(211)
plt.xlabel("Observed Wavelength ($\AA$)")
plt.ylabel("Flux")
plt.xlim([np.min(wave), np.max(wave)])

plt.plot(wave, pp.galaxy, 'k')
plt.axhline(y=-0, linestyle='--', color='k', linewidth=2)
plt.plot(wave, (stars), 'r', linewidth=2)  # overplot stellar templates alone
plt.plot(wave, (pp.galaxy-stars),'o', ms=1, color='LimeGreen', mec='LimeGreen')
plt.plot(wave, (gas), 'b', linewidth=1)  # overplot emission lines alone

plt.savefig(destdir+'/pPXFoutput_integral_spec.pdf', bbox_inches='tight')
plt.clf()

##### Defining all the outputs
solutions.append(pp.sol)
velF.append(pp.sol[0])
sigmaF.append(pp.sol[1])
velF.append(pp.sol[0])
sigmaF.append(pp.sol[1])
V_gas_vector = np.append(V_gas_vector,pp.sol[1][0])
V_stellar_vector = np.append(V_stellar_vector,pp.sol[0][0])
gas_sigma = np.append(gas_sigma,pp.sol[1][1])
stellar_sigma = np.append(stellar_sigma,pp.sol[0][1])

for name, weight, line in zip(line_names, pp.weights[w], pp.matrix[:,w].T):
    flux = weight*np.max(line)
    fluxes = np.append(fluxes, flux)

Hdelta, Hgamma, Hbeta, Halpha, SII_6717, SII_6731, OIII_5007d, OI6300d, NII, fluxes =\
np.append(Hdelta,fluxes[0]), np.append(Hgamma,fluxes[1]), np.append(Hbeta,fluxes[2]), np.append(Halpha,fluxes[3]),\
np.append(SII_6717,fluxes[4]), np.append(SII_6731,fluxes[5]), np.append(OIII_5007d,fluxes[6]),\
np.append(OI6300d,fluxes[7]), np.append(NII,fluxes[8]), []

print("pPXF is done running on the integrated spectra.")
filename = destdir+'/integral_spec_ppxf_velgas_output.csv'
dat = Table([V_gas_vector,V_stellar_vector, gas_sigma, stellar_sigma, \
                            Hdelta, Hgamma, Hbeta, Halpha, SII_6717, SII_6731, OIII_5007d, OI6300d, NII], \
                names=("V_gas_vector","V_stellar_vector", "gas_sigma", "stellar_sigma", \
                            "Hdelta", "Hgamma", "Hbeta", "Halpha","SII_6717", "SII_6731", "OIII_5007d", "OI6300d", "NII"))
dat.write(filename, format='csv', overwrite=True)
