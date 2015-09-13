# How to Model Optically-Thin Continuum Emisison From Debris Disks with an Affine-Invariant MCMC Method

This repository contains everything necessary to model your favorite debris disk with a Markov Chain Monte Carlo (MCMC) method using emcee (Foreman-Mackey et al. 2012, https://github.com/dfm/emcee) and code originally written by Noel Moldvai and Angelo Ricarte (Ricarte et al. 2013, https://github.com/noeleo/diskmodel). I edited these codes (disk_spl.py and visgen_master.py) in order to account for multiple disk components and added additional functions for analyzing the mid-infrared. In addition, I've appended codes for analyzing each MCMC run, for averaging together an IRS spectrum (if there is one) to be used in fitting to the Spectral Energy Distribution (SED), for de-redding photometry to line up nicely with a model photosphere, and for basic imaging. Both visibility-only and simultaneous visibility and SED fitting is possible. All the information is contained within. 

## Getting Started

To get started, you need the .uvf files of the data, a clone of this repository, and MIRIAD, astropy, and emcee installed on your machine. There was a CASA bug when I was working on this project related to exporting multiple spectral windows (spws) in one .uvf file, so it's currently set up to handle three separate spws called DISKNAMEspw#.uvf (in this repo, DISKNAME == HIP79516 as an example, and # goes from 1-3).

#### Making an Image

To make an image, simply call imageMaker.py in vis/ and enter the disk name:

```bash
python imageMaker.py
```

The various MIRIAD commands are talkative and will tell you some information about the imaging process, but things to note that are hard coded that you may want to change:
  * Cutoff = 100 microJy
  * Cell = 0.00625
  * Cell = 1024
  * Contours are [-3, 3, 6, 9, etc.] x 45 microJy (the RMS noise) in the first image that is displayed
The final MIRIAD command called by imageMaker.py calls imstat on a region of the image off the disk. The RMS noise reported here should be the contour level. The command below, filling in appropriately for DISKNAME and NOISE (in microJy), will give more accurate contour levels. 

```
cgdisp in=DISKNAMEspw123.cm,DISKNAMEspw123.cm type=pix,con slev=a,NOISEe-6 levs1=-3,3,6,9,12,15,18 labtyp=arcsec device=/xs beamtyp=b,l,4
```

#### Prepping Your Data/ Gathering Stellar Info

##### Stellar Parameters

For stellar parameters, go to [this version of the Hipparcos catalog](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/MNRAS/427/343). Look up your star and note the the luminosity, effective temperature, and distance. We can derive the star's radius from the Stefan-Boltzmann law. Wolfram alpha is a big help for natural language inqueries and dealing with units, so after you prove to yourself what's going on, I highly recommend going to wolframalpha.com and entering (where the luminosity and temperature are the star's temperature, of course):

```
(99999999 solar luminosities / (4 pi * stefan boltzmann constant * (1000000000K)^4))^.5 / solar radius
```

These values need to be input into disk_spl.py. See lines 34 (effective temperature), 36 (stellar radius), 38 (luminosity), and 40 (distance). 

##### Model Stellar Spectrum/ Photometry

You'll need a stellar spectrum and photometry at minimum, and if there is an IRS spectrum, you'll want that too. 

Go to [this Kurucz-Lejeune stellar photosphere model catalog](http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+AS/125/229) and input the star's effective temperature rounded to the nearest 250K. [M/H] = 0.0 and log(g) = 4.5 are fine assumptions if you don't know the specifics about the star. Edit these choices if more information is available of course! Select the .corr file, then click "values as a table." Copy these data into a text file called Model_Photosphere.txt where the fluxes are non-zero-- usually starting at about 50 nanometers! This text file is called by disk_spl.py on line 74. 

If you've got an IRS spectrum, make sure it's in three columns in a text file (wavelength of obs, flux, and error), and use it as input for avgIRS.py (line 14). The spectrum is generally ~360 points. This file averages them to 9 for sake of computational efficiency when fitting the SED (the main use of processing power when making models) and adds a 10% systematic flux uncertainty. Set the final range in "ranges" on line 34 to correspond to the last data point in the IRS spectrum. See the length of the text file to deduce what number point this is. avgIRS.py will spit out a text file of the averaged points (named on line 52) and save a plot (line 62) that confirms that the data were averaged appropriately. [Here](https://github.com/jliemansifry/DustDiskModeling/blob/master/HIP79516_IRS_spectrum_avg.png) is the example in this repo.

Gather your photometry. Insert the IRS points into this file by hand ordered by wavelength. If there is a reddening factor, you should correct for this so that your photometry will happily agree with your model photosphere. Open up ISM_extinction.py, set your aV (line 6), input and output filenames (lines 11 and 18), and run it! Make sure that you've already added any systematic errors before doing this. 

## MCMC Modeling

Time to do some MCMC modeling! If you want/need more details than those offered here, please see Ricarte et al. 2013 and Lieman-Sifry et al. (in prep). Email me if you want a copy of the paper draft.

I highly recommend doing a visibility-only fit the first time around in order to deduce what the radio data say about your disk before attempting to simultaneously recreate the fluxes of the SED. Open up doModel_vis_only.py. We describe a blowout grain size (essentially the smallest grains that can live on stable orbits around the star without being blown away) in doing these visibility-only fits, as grain size is degenerate with disk mass. Fill in appropriate values for your star. Define starting center positions and widths for the gaussian ball of the parameters you're going to vary on lines 16 and 17. 

```python
pC = np.array([60, 50, -2, 0, 0]) # center positions for the gaussian ball for each parameter to be varied
pW = np.array([20, 30, 0.5, 30, 30]) # gaussian widths for each parameter
```

In this example, the 5 parameters being varied are:
  1. Inner edge of the disk (innerRadius)
  2. Width of the disk (delR)
  3. log(diskMass)
  4. Inclination (i)
  5. Position angle (PA)

In the log likelihood function, these values are called by the Disk instance on line 28 and VisibilityGenerator instance on line 30. The 0.0001 refers to the inner edge of an inner belt that we may or may not need yet, with the 0 being its corresponding mass. 0.8 is the long wavelength power law index of grain emission efficiency (beta), and 1 is the power law (p) slope of the surface density (decreasing). If your data are good enough to warrant varing beta and p, your calls on lines 28 and 30 would look something like:

```python
disk = Disk(0.0001,p[0],p[1],10**p[2], 0, p[3], blowoutSize*1e6, p[4])
visGen = VisibilityGenerator(1024,p[5],p[6],fits_file) # renumbering which values correspond to inc and PA because it makes sense to keep them at the end
```

If you add parameters, make sure you edit lines 50 and 51, adjust the boundaries for the walkers on line 22, and set start positions/widths on lines 16 and 17. I recommend starting the first MCMC run really wide, then narrowing in after you've got an idea where the best fit is. Set the number of walkers/steps on line 18. Start small(ish)... however, you'll want enough walkers that they'll converge quickly and enough steps that the run will burn in. I can't provide exacts here-- you'll want to play around with this based on how much time you've got/how many cores you can run your code on. Speaking of running the code...

```bash
mpirun -np NUMBEROFCORES python doModel_vis_only.py
```

Watch the output go for a little bit... Make sure that nothing complains. Also, make sure it is printing chi-squareds on the order of the number of visibilities you have. If it is printing inf and continues to print inf, your visibility files are probably corrupted. Check back in occasionally, especially if you've changed things. Sometimes, if the boundaries not appropriate and a parameter wonders off into la-la-land (power laws of 30, huge disk masses, etc.), the code will freeze, but only on one core, and everything else will try to keep going. If you just glance at the output from the code, it might *look* like it's all good, but really it is frozen. Scroll up a bit to make sure everyone is chugging along happily, or wait until it starts to print some more stuff. If it is frozen, look carefully at the values of the parameters that made it freeze and the error message! Hopefully you'll be able to track down what is going wrong. 
