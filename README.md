# How to Model Optically-thin Continuum Emisison from Debris Disks

This repository contains everything necessary to model your favorite debris disk with a Markov Chain Monte Carlo (MCMC) method using emcee (Foreman-Mackey et al. 2012, https://github.com/dfm/emcee) and code originally written by Noel Moldvai and Angelo Ricarte (Ricarte et al. 2013, https://github.com/noeleo/diskmodel). I edited these codes (disk_spl.py and visgen_master.py) in order to account for multiple disk components and added additional functions for analyzing the mid-infrared. In addition, I've appended codes for analyzing each MCMC run, for averaging together an IRS spectrum (if there is one) to be used in fitting to the Spectral Energy Distribution (SED), for de-redding photometry to line up nicely with a model photosphere, and for basic imaging. Both visibility-only and simultaneous visibility and SED fitting is possible. All the information is contained within. 

## Getting Started

To get started, you need the .uvf files of the data, a clone of this repository, and MIRIAD installed on your machine. There was a CASA bug when I was working on this project related to exporting multiple spectral windows (spws) in one .uvf file, so it's currently set up to handle three separate spws called DISKNAMEspw#.uvf (in this repo, DISKNAME == HIP79516 as an example, and # goes from 1-3).

#### Making an Image

To make an image, simply call imageMaker.py in vis/ and enter the disk name:

```bash
python imageMaker.py
What's the name of this disk? HIP79516
```

The various MIRIAD commands are talkative and will tell you some information about the imaging process, but things to note that are hard coded that you may want to change:
  * Cutoff = 100 microJy
  * Cell = 0.00625
  * Cell = 1024
  * Contours are [-3, 3, 6, 9, etc.] x 45 microJy (the RMS noise) in the first image that is displayed
The final MIRIAD command called by imageMaker.py calls imstat on a region of the image off the disk. The RMS noise reported here should be the contour level. The command below, filling in appropriately for DISKNAME and NOISE (in microJy), will give more accurate contour levels. 

```bash
cgdisp in=DISKNAMEspw123.cm,DISKNAMEspw123.cm type=pix,con slev=a,NOISEe-6 levs1=-3,3,6,9,12,15,18 labtyp=arcsec device=/xs beamtyp=b,l,4
```

#### Prepping your data/ gathering stellar info

##### Stellar Parameters

For stellar parameters, go to [this version of the Hipparcos catalog](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/MNRAS/427/343). Look up your star and note the the luminosity, effective temperature, and distance. We can derive the star's radius from the Stefan-Boltzmann law. Wolfram alpha is a big help for natural language inqueries and dealing with units, so after you prove to yourself what's going on, I highly recommend going to wolframalpha.com and entering (where the luminosity and temperature are the star's temperature, of course):

```bash
(99999999 solar luminosities / (4 pi * stefan boltzmann constant * (1000000000K)^4))^.5 / solar radius
```

These values need to be input into disk_spl.py. See lines 34 (effective temperature), 36 (stellar radius), 38 (luminosity), and 40 (distance). 

##### Model Stellar Spectrum/ Photometry

You'll need a stellar spectrum and photometry at minimum, and if there is an IRS spectrum, you'll want that too. 

Go to [this Kurucz-Lejeune stellar photosphere model catalog](http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+AS/125/229) and input the star's effective temperature rounded to the nearest 250K. [M/H] = 0.0 and log(g) = 4.5 are fine assumptions if you don't know the specifics about the star. Edit these choices if more information is available of course! Select the .corr file, then click "values as a table." Copy these data into a text file called Model_Photosphere.txt where the fluxes are non-zero-- usually starting at about 50 nanometers! This text file is called by disk_spl.py on line 74. 

If you've got an IRS spectrum, make sure it's in three columns in a text file (wavelength of obs, flux, and error), and use it as input for avgIRS.py (line 14). The spectrum is generally ~360 points. This file averages them to 9 for sake of computational efficiency when fitting the SED (the main use of processing power when making models) and adds a 10% systematic flux uncertainty. Set the final range in "ranges" on line 34 to correspond to the last data point in the IRS spectrum. See the length of the text file to deduce what number point this is. avgIRS.py will spit out a text file of the averaged points (named on line 52) and save a plot (line 62) that confirms that the data were averaged appropriately. [Here](https://github.com/jliemansifry/DustDiskModeling/blob/master/HIP79516_IRS_spectrum_avg.png) is the example in this repo.

Gather your photometry. Insert the IRS points into this file by hand ordered by wavelength. If there is a reddening factor, you should correct for this so that your photometry will happily agree with your model photosphere. Open up ISM_extinction.py, set your aV (line 6), input and output filenames (lines 11 and 18), and run it! Make sure that you've already added any systematic errors before doing this. 







