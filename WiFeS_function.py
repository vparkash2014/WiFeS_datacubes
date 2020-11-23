import numpy as np
import ppxf.ppxf_util as util
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from matplotlib import gridspec
import os, glob, time, random, math
from astropy.table import Table
from scipy import ndimage, optimize
from sys import stdout
import scipy.interpolate
import math
from matplotlib.colors import ListedColormap


def read_fits(filename, limits = [0,0,0,0]):
    """
    Read datacube
    Input: filname and the y and x boundaries: limits = [ymin,ymax,xmin,xmax]
    """
    fits_cube = fits.open(filename)
    cube, hdr, var = fits_cube[0].data, fits_cube[0].header, fits_cube[1].data
    if sum(limits) > 0:
        cube, var = cube[:,limits[0]:limits[1],limits[2]:limits[3]], var[:,limits[0]:limits[1],limits[2]:limits[3]]
    fits_cube.close()

    """
    Calibrate in wavelength from header keywords
    """
    cdelt3 = hdr['CDELT3'] # the scale
    crval3 = hdr['CRVAL3'] # the starting wavelength
    naxis3 = hdr['NAXIS3'] # the length of the axis
    pixel = np.arange(cube.shape[0])
    lamb = crval3 + (pixel) * cdelt3 #the wavelength ranges
    lamRange =  [crval3 + 0,crval3 + cdelt3*(naxis3-1)] # the first and last wavelength of the spectrum

    return (cube, hdr, var, lamb, lamRange)

def centroid(cube, PA, ba):
    """
    This function finds the centre galaxy (based of Sextractor)
    """
    frame = np.median(cube, axis=0)
    x_cent,y_cent  = np.where(frame == np.amax(frame))[1], np.where(frame == np.amax(frame))[0]

    center = optimize.fmin(callMeasureVariation, [x_cent,y_cent],args = (PA, ba, frame), xtol=1./x_cent, ftol=0.1, full_output=1, disp=1)

    return center

def callMeasureVariation(pos, PA, ba, cubeR):
    """
    This function see how best the centre of a galaxy was estimated.
    pos = [x,y], map= is the 2-D image of the galaxy
    """
    totStd = measureVariation(pos, cubeR, PA, ba, binSize=0.2, lim=1.)
    return totStd

def measureVariation(posMax, cubeR, PA, ba, binSize=0.1, lim=1.2):
    """
    This function was written by Amelia Fraser-McKelvie and is based off Sextractor
    """
    listR, listInten = [], []
    for ii in np.arange(len(cubeR)):
        for jj in np.arange(len(cubeR[ii])):
            listInten.append(cubeR[ii][jj])
            angleRot = (np.pi/180.)*(PA-90.)
            xRot = ((float(ii)-posMax[0])*np.cos(angleRot)-(float(jj)-posMax[1])*np.sin(angleRot))
            yRot = ((float(ii)-posMax[0])*np.sin(angleRot)+(float(jj)-posMax[1])*np.cos(angleRot))
            listR.append(np.sqrt(ba*(xRot**2.)+(yRot**2.)/ba))

    indices = permutation_indices(listR)
    listR_log = np.log10(np.array(listR)[indices])
    listInten_log = np.log10(np.array(listInten)[indices])
    #Binning
    listStd = 0
    listR_log[np.nonzero(np.isinf(listR_log))] = -2
    for ii in np.arange(np.min(listR_log), listR_log[-1], binSize):
        listTMP = []
        for jj in np.arange(len(listR_log)):
            if ii <= listR_log[jj] < (ii + binSize):
                listTMP.append(listInten_log[jj])
        listStd += np.std(listTMP)
    return listStd

def permutation_indices(data):
    return sorted(range(len(data)), key = data.__getitem__)

def singal_noise_var_sec(cube, var, lamRange, lamb, wavelength_range):
    """
    This function calculates a mean signal and noise of each spaxel using the var
    output: 2 1D array listing all the signal and noise for the spaxels
    """
    signal, noise = [],[]

    x_1d, y_1d = np.array([x for x in range( cube.shape[2] ) for y in range( cube.shape[1] )]), np.array([y for x in range( cube.shape[2] ) for y in range( cube.shape[1] )])

    for x,y in zip(x_1d, y_1d) :
        specNew, logLam, velscale = rebin_spectra(cube, x, y, lamRange, lamb)
        var_int = var[:, int(y), int(x)]
        lamNew = np.exp(logLam)
        specNew, var_int  = \
        specNew[np.where( (lamNew >= wavelength_range[0]) &( lamNew <= wavelength_range[1]) )], \
        var_int[np.where( (lamNew >= wavelength_range[0]) &( lamNew <= wavelength_range[1]) )]

        signal_range, noise_range = np.mean(specNew), np.mean(np.sqrt(var_int))

        signal.extend([signal_range])
        noise.extend([noise_range])

    return(np.array(signal), np.array(noise))

def voronoi_outputs(cube1, signal_1, noise_1, signal_2, noise_2, xNode_1, xNode_2, binNum_1, binNum_2, output, gas = 0):
    """
    Display the results from the voronoi binning.
    Update: 14 June 2018: I added a key to deal with the gas data cube that doesn't use two arms, so the first two plots will
    show the signal map and S/N mar
    """
    x_1d, y_1d = np.array([x for x in range( cube1.shape[2] ) for y in range( cube1.shape[1] )]), np.array([y for x in range( cube1.shape[2] ) for y in range( cube1.shape[1] )])

    "setting up the plot"
    fig, axs = plt.subplots(4,2, figsize=(4, 10))
    fig.subplots_adjust(hspace=.5)
    plt.setp(axs, xticklabels=[])
    plt.setp(axs, yticklabels=[])
    "the raw galaxy"
    img_1, img_2 = np.full(( cube1.shape[1], cube1.shape[2]), np.nan), np.full((cube1.shape[1], cube1.shape[2]), np.nan)

    img_1[y_1d, x_1d] = signal_1
    img_2[y_1d, x_1d] = signal_2
    if gas == 1:
        img_1[y_1d, x_1d] = signal_1
        SN = np.full(( cube1.shape[1], cube1.shape[2]), np.nan)
        SN[y_1d, x_1d] = signal_1/noise_1
        img_2 = SN

    img0 = axs[0,0].imshow(img_1, interpolation='nearest', origin='upper', cmap='rainbow', vmin=np.min(img_1), vmax=np.max(img_1) )
    fig.colorbar(img0, ax=axs[0, 0])
    img1 = axs[0,1].imshow(img_2, interpolation='nearest', origin='upper', cmap='rainbow', vmin=np.min(img_2), vmax=np.max(img_2) )
    fig.colorbar(img1, ax=axs[0, 1])

    "The voronoi bins"

    rnd_1, rnd_2 = np.argsort(np.random.random(xNode_1.size)), np.argsort(np.random.random(xNode_2.size))
    counts_1, counts_2 = rnd_1[binNum_1], rnd_2[binNum_2]
    img_1, img_2 = np.full(( cube1.shape[1], cube1.shape[2]), np.nan),\
    np.full((cube1.shape[1], cube1.shape[2]), np.nan)
    img_1[y_1d, x_1d] = counts_1
    img_2[y_1d, x_1d] = counts_2

    if gas == 0:
        im2 = axs[1,0].imshow(img_1, interpolation='nearest', cmap='prism')
        im3 = axs[1,1].imshow(img_2, interpolation='nearest', cmap='prism')
    if gas == 1:
        im2 = axs[1,0].imshow(img_1, interpolation='nearest', cmap='prism')

    # After binning
    bin_signal_1 = [np.mean(signal_1[np.where(binNum_1 == ind)]) for num in binNum_1 for ind in np.unique(binNum_1) if num == ind]

    bin_signal_2 = [np.mean(signal_2[np.where(binNum_2 == ind)]) for num in binNum_2 for ind in np.unique(binNum_2) if num == ind]

    img_1, img_2 = np.full(( cube1.shape[1], cube1.shape[2]), np.nan), np.full((cube1.shape[1], cube1.shape[2]), np.nan)

    img_1[y_1d, x_1d], img_2[y_1d, x_1d] = bin_signal_1, bin_signal_2

    img4 = axs[2,0].imshow(img_1, interpolation='nearest', origin='upper', cmap='rainbow', vmin=np.min(img_1), vmax=np.max(img_1) )
    fig.colorbar(img4, ax=axs[2, 0])
    if gas == 0:
        img5 = axs[2,1].imshow(img_2, interpolation='nearest', origin='upper', cmap='rainbow', vmin=np.min(img_2), vmax=np.max(img_2) )
        fig.colorbar(img5, ax=axs[2, 1])

    "plotting the difference"
    diff_1, diff_2 = signal_1 - bin_signal_1, signal_1 - bin_signal_1

    img_1, img_2 = np.full(( cube1.shape[1], cube1.shape[2]), np.nan),\
    np.full((cube1.shape[1], cube1.shape[2]), np.nan)

    img_1[y_1d, x_1d], img_2[y_1d, x_1d] = diff_1, diff_2

    img6 = axs[3,0].imshow(img_1, interpolation='nearest', origin='upper', cmap='rainbow', vmin=np.min(img_1), vmax=np.max(img_1) )
    fig.colorbar(img6, ax=axs[3, 0])
    if gas == 0:
        img7 = axs[3,1].imshow(img_2, interpolation='nearest', origin='upper', cmap='rainbow', vmin=np.min(img_2), vmax=np.max(img_2))
        fig.colorbar(img7, ax=axs[3, 1])

    plt.savefig(output, bbox_inches='tight')
    return fig

def voronoi_spectra_outputs(cubeR, cubeB, lamRangeR, lambR, binNum_R, lamRangeB, lambB, binNum_B, output):
    """
    This function is going to compare the individual spaxel spectrum to the voronoi bin spectrum
    """

    "First we are going to pick 10 random spaxels in the data cube"
    x_1d, y_1d = np.array([x for x in range( cubeR.shape[2] ) for y in range( cubeR.shape[1] )]),\
                 np.array([y for x in range( cubeR.shape[2] ) for y in range( cubeR.shape[1] )])
    num = 10
    y_rand, x_rand = np.random.randint(0, cubeR.shape[1],num), np.random.randint(0,cubeR.shape[2], num)

    "Setting up the Figure"
    fig = plt.subplots(figsize = (17,70))
    gs = gridspec.GridSpec(11, 2, width_ratios=[3, 1])

    frame_R = np.mean(cubeR, axis=0)

    for i in range(num):
        ax0 = plt.subplot(gs[2*i])
        ax1 = plt.subplot(gs[2*i+1])

        specNewR, logLamR, velscaleR = rebin_spectra(cubeR, x_rand[i], y_rand[i], lamRangeR, lambR)
        specNewB, logLamB, velscaleB = rebin_spectra(cubeB, x_rand[i], y_rand[i], lamRangeB, lambB)

        lamNewR, lamNewB = np.exp(logLamR[np.where(specNewR > 0)]), np.exp(logLamB[np.where(specNewB > 0)])
        specNewR, specNewB = specNewR[np.where(specNewR > 0)], specNewB[np.where(specNewB > 0)]
        ax0.plot(lamNewR, specNewR, sns.xkcd_rgb["pale red"], linewidth = 2)
        ax0.plot(lamNewB, specNewB, sns.xkcd_rgb["denim blue"], linewidth = 2)
        ax1.imshow(frame_R, interpolation='nearest', origin='upper', cmap='rainbow', vmin=np.min(frame_R), vmax=np.max(frame_R) )
        ind_R, ind_B = binNum_R[int(np.where( (x_1d == x_rand[i]) & (y_1d == y_rand[i]) )[0])], \
              binNum_B[int(np.where( (x_1d == x_rand[i]) & (y_1d == y_rand[i]) )[0])]

        x_ind_R, y_ind_R, x_ind_B, y_ind_B = x_1d[np.where(np.array(binNum_R) == ind_R)],\
                                                  y_1d[np.where(np.array(binNum_R) == ind_R)], \
                                                  x_1d[np.where(np.array(binNum_B) == ind_B)],\
                                                  y_1d[np.where(np.array(binNum_B) == ind_B)]

        specNewR, logLamR, velscaleR = rebin_spectra(cubeR, x_ind_R, y_ind_R, lamRangeR, lambR)
        specNewB, logLamB, velscaleB = rebin_spectra(cubeB, x_ind_B, y_ind_B, lamRangeB, lambB)
        lamNewR, lamNewB = np.exp(logLamR[np.where(specNewR > 0)]), np.exp(logLamB[np.where(specNewB > 0)])
        specNewR, specNewB = specNewR[np.where(specNewR > 0)], specNewB[np.where(specNewB > 0)]

        ax0.plot(lamNewR, specNewR, "k", linestyle= '--', linewidth = 1)
        ax0.plot(lamNewB, specNewB, "k", linestyle= '--', linewidth = 1)
        ax1.scatter(x_ind_R, y_ind_R,c= sns.xkcd_rgb["pale red"], marker = 's', alpha = 1)
        ax1.scatter(x_ind_B, y_ind_B,c= sns.xkcd_rgb["denim blue"], marker = 's', alpha = 1)
        ax1.scatter(x_rand[i],y_rand[i],c = 'k')
    plt.savefig(output, bbox_inches='tight')
    return fig

def rebin_spectra(cube, x, y, lamRange, lamb):
    """
    Extracts the a single spectrum for a data given the desired x,y by binning the data using pPXF util.log_rebin
    """
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    spectrum = np.median(np.array([cube[:, int(ypair), int(xpair)] for ypair,xpair in zip(y,x)]), axis = 0)
    c = 299792.458 # speed of light in km/s
    specNew, logLam, velscale = util.log_rebin(lamRange, spectrum, oversample=False, velscale=c*np.log(lamb[1]/lamb[0]), flux=False)
    return (specNew, logLam, velscale)

def single_spectra(cube, x, y):
    """
    Extracts the a single spectrum for a data given the desired x,y without rebinning
    """
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    spectrum = np.median(np.array([cube[:, int(ypair), int(xpair)] for ypair,xpair in zip(y,x)]), axis = 0)
    return (spectrum)

def glue_BR(fluxB, lambB, fluxR, lambR):
    """
    'glues' together the Blue and red arm by scaling the red arm according to the blue velscaleR
    """

    """First we are figuring out the difference between the two arm. aka how much does the red arm need to be shifted to align with the blue arm"""
    lim = [lambR[1],lambB[-1]]
    waveB, galB, waveR, galR = \
        lambB[np.where(lambB > lim[0])], fluxB[np.where(lambB > lim[0])],\
        lambR[np.where(lambR < lim[-1])], fluxR[np.where(lambR < lim[-1])]
    diff = []
    for wave in waveR:
        x = galR[np.where(waveR == wave)] - np.mean(galB[np.where(abs(waveB - wave) <= 2.0)])
        if math.isnan(x) == False:
            diff = np.append( diff,galR[np.where(waveR == wave)] - np.mean(galB[np.where(abs(waveB - wave) <= 2.0)]) )
    diff = np.mean(diff)
    fluxR[:] = [x - diff for x in fluxR]

    """Increasing the resolution of the red arm to the same resolution of the blue arm to create on spectra"""
    lambscale = np.mean(np.diff(lambB))
    new_wave = np.arange(lambB[0], lambR[-1], lambscale)
    new_flux = fluxB
    diff_lim = 3.0
    for wave in new_wave:
        if wave > new_wave[len(lambB)-1]:
            flux, lamb = fluxR[np.where(abs(lambR-wave) <= diff_lim)], lambR[np.where(abs(lambR-wave) <= diff_lim)]
            flux_interp = scipy.interpolate.interp1d(lamb, flux)
            new_flux = np.append( new_flux, flux_interp(wave))
    return new_flux, new_wave


def retrieveTemplates(velscale, FWHM_WiFeS_ang, range_sel = None):
    """
    This function reads and sets up the stellar libraries for pPXF
    This function was written by Amelia Fraser-McKelvie
    """
    templatesInput = glob.glob('/Users/vaishalp/Dropbox/Amelia/IC1059/Archive/StellarLibraries/INDO-US_Valdes/*.txt')
    FWHM_tem_ang = 1.35 #In \AA

    #SINGLE TEMPLATE
    wv_temp1, inten_temp1 = readAsciiSpec(templatesInput[0])
    if range_sel:
        range_sel_rest = range_sel
    else:
        range_sel_rest = [np.min(wv_temp1), np.max(wv_temp1)]

    #Selection Chunk of template
    sel_temp = np.nonzero((wv_temp1 >= range_sel_rest[0]) & (wv_temp1 <= range_sel_rest[1]))
    wv_temp, inten_temp = wv_temp1[sel_temp], inten_temp1[sel_temp]
    sspNew, logLam2, velscale = util.log_rebin([wv_temp[0], wv_temp[-1]], inten_temp, velscale=velscale)

    #ALL TEMPLATES
    templates = np.empty((sspNew.size, len(templatesInput)))
    #convolution
    FWHM_dif = np.sqrt(FWHM_WiFeS_ang**2 - FWHM_tem_ang**2)
    sigma = FWHM_dif/2.355/(wv_temp[1] - wv_temp[0]) # Sigma difference in pixels
    for ii in range(len(templatesInput)):
        wv_temp2, inten_temp2 = readAsciiSpec(templatesInput[ii])
        wv_temp_ii, inten_temp_ii = wv_temp2[sel_temp], inten_temp2[sel_temp]
        ssp = ndimage.gaussian_filter1d(inten_temp_ii, sigma) #One-dimensional Gaussian filter.
        sspNew, logLam2, velscale = util.log_rebin([wv_temp[0], wv_temp[-1]], ssp, velscale=velscale)

        templates[:,ii] = sspNew/np.median(sspNew)
        stdout.write("\rDONE %i/%i" % (ii+1, len(templatesInput)))
        stdout.flush()

    return logLam2, templates
