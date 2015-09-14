# How to Model Optically-Thin Continuum Emisison From Debris Disks with an Affine-Invariant MCMC Method

This repository contains everything necessary to model your favorite debris disk with a Markov Chain Monte Carlo (MCMC) method using emcee (Foreman-Mackey et al. 2012, https://github.com/dfm/emcee) and code originally written by Noel Moldvai and Angelo Ricarte (Ricarte et al. 2013, https://github.com/noeleo/diskmodel). I edited these codes (disk_as.py and visgen.py in that repo) in order to account for multiple disk components and added additional functions for analyzing the mid-infrared. In addition, I've appended codes for analyzing each MCMC run, for averaging together an IRS spectrum (if there is one) to be used in fitting to the Spectral Energy Distribution (SED), for de-redding photometry to line up nicely with a model photosphere, and for basic imaging. Both visibility-only and simultaneous visibility and SED fitting is possible. All the information is contained within. 

## Getting Started

To get started, you need the .uvf files of the data, a clone of this repository, and MIRIAD, astropy, and emcee installed on your machine. There was a CASA bug when I was working on this project related to exporting multiple spectral windows (spws) in one .uvf file, so it's currently set up to handle three separate spws called DISKNAMEspw#.uvf (in this repo, DISKNAME = HIP79516 as an example, and # goes from 1-3).

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
The final MIRIAD command called by imageMaker.py calls imstat on a region of the image off the disk. The RMS noise reported here should be the contour level. The command below, filling in appropriately for DISKNAME and NOISE (in microJy), will give appropriate contour levels. 

```
cgdisp in=DISKNAMEspw123.cm,DISKNAMEspw123.cm type=pix,con slev=a,NOISEe-6 levs1=-3,3,6,9,12,15,18 labtyp=arcsec device=/xs beamtyp=b,l,4
```

#### Prepping Your Data/ Gathering Stellar Info

##### Stellar Parameters

For stellar parameters, go to [this version of the Hipparcos catalog](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/MNRAS/427/343). Look up your star and note the the luminosity, effective temperature, and distance. We can derive the star's radius from the Stefan-Boltzmann law. Wolfram alpha is a big help for natural language inqueries and dealing with units, so after you prove to yourself what's going on, I highly recommend going to wolframalpha.com and entering (where the luminosity and temperature are your star's characteristics, of course):

```
(4.04 solar luminosities / (4 pi * stefan boltzmann constant * (6500K)^4))^.5 / solar radius
```

These values need to be input into disk_spl.py. See lines 34 (effective temperature), 36 (stellar radius), 38 (luminosity), and 40 (distance). 

##### Model Stellar Spectrum/ Photometry

You'll need a stellar spectrum and photometry at minimum, and if there is an IRS spectrum, you'll want that too. 

Go to [this Kurucz-Lejeune stellar photosphere model catalog](http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+AS/125/229) and input the star's effective temperature rounded to the nearest 250K. [M/H] = 0.0 and log(g) = 4.5 are fine assumptions if you don't know the specifics about the star. Edit these choices if more information is available though, of course! Select the .corr file, then click "values as a table." Copy these data into a text file called Model_Photosphere.txt where the fluxes are non-zero-- usually starting at about 50 nanometers. This text file is called by disk_spl.py on line 74. 

If you've got an IRS spectrum, make sure it's in three columns in a text file (wavelength of obs, flux, and error), and use it as input for avgIRS.py (line 14). The spectrum is generally ~360 points. This file averages them to 9 for sake of computational efficiency when fitting the SED (the main use of processing power when making models) and adds a 10% systematic flux uncertainty. Set the final range in "ranges" on line 34 to correspond to the last data point in the IRS spectrum. See the length of the text file to deduce what number point this is. avgIRS.py will spit out a text file of the averaged points (named on line 52) and save a plot (line 62) that confirms that the data were averaged appropriately. [Here](https://github.com/jliemansifry/DustDiskModeling/blob/master/HIP79516_IRS_spectrum_avg.png) is the example for HIP79516's IRS spectrum.

Gather your photometry from the literature. Insert the IRS points into this file by hand ordered by wavelength. If there is a reddening factor, you should correct for this so that your photometry will happily agree with your model photosphere. Open up ISM_extinction.py, set your aV (line 6), input and output filenames (lines 11 and 18), and run it! Make sure that you've already added any systematic errors before doing this. 

## MCMC Modeling

Time to do some MCMC modeling! If you want/need more details about the modeling process itself, please see Ricarte et al. 2013 and Lieman-Sifry et al. (in prep). Email me if you want a copy of the paper draft. If you want more details about emcee, see Foreman-Mackey et al. 2012 or [the online documentation](http://dan.iel.fm/emcee/current/), which is really good. If you want more details about the style of MCMC used by emcee, see Goodman & Weare (2010). 

#### Visibility-Only Fitting

I highly recommend doing a visibility-only fit the first time around in order to deduce what the radio data say about your disk before attempting to simultaneously recreate the fluxes of the SED. disk_spl.py describes the surface density distribution of the dust disk with a single power law (thus the _spl... more complex descrptions can be found in [my thesis](https://www.dropbox.com/s/suvi5gqk2bdmg4s/JesseLiemanSifry_Thesis.pdf?dl=0). See any of the references mentioned thus far for a more complete description of the specifics of the disk model, and see my thesis for an excessively complete description of the details.

Open up doModel_vis_only.py. We describe a blowout grain size (essentially the smallest grains that can live on stable orbits around the star without being blown away) in doing these visibility-only fits, as grain size is degenerate with disk mass when not simultaneously fitting to the SED. Fill in appropriate values for your star. Define starting center positions and widths for the gaussian ball of the parameters you're going to vary on lines 16 and 17. 

```python
pC = np.array([60, 50, -2, 0, 0]) # center positions for the gaussian ball for each parameter to be varied
pW = np.array([20, 30, 0.5, 30, 30]) # gaussian widths for each parameter
```

In this example, the 5 parameters being varied are:
  1. Inner edge of the disk (*rIn*)
  2. Width of the disk (*delR*)
  3. log(*mDisk*)
  4. Inclination (*i*)
  5. Position angle (*PA*)

In the log likelihood function, these values are called by the Disk instance on line 28 and VisibilityGenerator instance on line 30. The 0.0001 refers to the inner edge of an inner belt that we may or may not need yet, with the 0 being its corresponding mass. 0.8 is the long wavelength power law index of grain emission efficiency (beta), and 1 is the power law (p) slope of the surface density (decreasing). If your data are good enough to warrant varing beta and p, your calls on lines 28 and 30 would look something like:

```python
disk = Disk(0.0001, p[0], p[1], 10**p[2], 0, p[3], blowoutSize*1e6, p[4])
visGen = VisibilityGenerator(1024, p[5], p[6], fits_file) # renumbering which values correspond to inc and PA because it makes sense to keep them at the end
```

If you add parameters, make sure you edit lines 50 and 51, adjust the boundaries for the walkers on line 22, and set start positions/widths on lines 16 and 17. I recommend starting the first MCMC run really wide, then narrowing in after you've got an idea where the best fit is. Set the number of walkers/steps on line 18. Start small(ish)... however, you'll want enough walkers that they'll converge quickly and enough steps that the run will burn in. I can't provide exacts here-- you'll want to play around with this based on how much time you've got/how many cores you can run your code on. To run the code, run this on the command line (with NUMBEROFCORES = anything greater than 1). 

```bash
mpirun -np NUMBEROFCORES python doModel_vis_only.py
```

Watch the output go for a little bit... Make sure that nothing complains. Be sure that your computer settings are such that only your display goes to sleep, and not your processors! Make sure it is printing chi-squareds on the order of the number of visibilities you have in between all the MIRIAD stuff. If it is printing inf and continues to print inf, your visibility files are probably corrupted. Check back in occasionally, especially if you've changed things. Sometimes, if the boundaries not appropriate and a parameter wonders off into la-la-land (power laws of 30, huge disk masses, etc.), the code will freeze, but only on one core, and everything else will try to keep going. If you just glance at the output from the code, it might *look* like it's all good, but really it is frozen. Scroll up a bit to make sure everyone is chugging along happily, or wait until it starts to print some more stuff. If it is frozen, look carefully at the values of the parameters that made it freeze and the error message! Hopefully you'll be able to track down what is going wrong. Remember that if your parameters are wandering into la-la-land, even though you might believe that they're constrained, they don't believe they're constrained. Check your code. 

#### Analyzing the Results of Your Visibility-Only MCMC Run

mcmcPlaytime_vis_only.py provides a bunch of handy tools for analyzing how well your model fitting is going. Make sure your blowoutSize is the same as described in doModel_vis_only.py. mcmcPlaytime will ask you "MCMC what by what?" to which you should simply respond NUMWALKERSxNUMSTEPS, i.e. "20x100" or whatever you asked for in doModel_vis_only.py on line 18 (this is helping it find the files that doModel created). It will ask you for the disk name, which will be used in output file naming, and the RMS noise, which you should know from when you first imaged the disk (quick note on imaging: make sure the param you used for imaging in imageMaker.py are the same as the ones in the functions miriad_residual and model_maker). 

You will then be asked if you want a residual image, followed by a model image. These will be displayed to screen rather than saved, as the postscript files that MIRIAD generates are upside-down and in negative colors. Take a screenshot if you please. The images are all saved in MCMCRUNS/vis_only/NUMWALKERSxNUMSTEPS/ and could be made into pretty model/residual plots with just slight edits of [this code](https://github.com/jliemansifry/49_Ceti/blob/master/49Ceti_dataPlots.py). You'll just need to create .fits versions of the model and residual, which can be done quickly with the code below. Make sure you're in the appropriate directory (i.e. MCMCRUNS/vis_only/NUMWALKERS/NUMSTEPS). 

```bash
fits in=NUMWALKERSxNUMSTEPS_BestFit_image.cm op=xyout out=NUMWALKERSxNUMSTEPS_BestFit_image.fits
fits in=NUMWALKERSxNUMSTEPS_BestFit_Residual.cm op=xyout out=NUMWALKERSxNUMSTEPS_BestFit_Residual.fits
```

The code will also prompt you as to whether you'd like an SED generated with the best-fit parameters from the vis-only fit, and while this may sound silly, it's actually quite enlightening. You can quickly see whether or not you'll need a secondary inner belt component (see page 32 of my thesis for an example).

Finally, we get to make a triangle plot with emcee's awesome [triangle](https://github.com/dfm/triangle.py) module and print out some numbers. Take a look at the raw chi-squared vs. step number for all the walkers in the plot that appears. The run is "burned in" when the raw chi-squareds have stopped falling precipitiously. See page 23 of my thesis for an example- in that case, the run is burned in after ~150 steps. If you have multiple horizontal lines around different raw chi-squareds, some of your walkers have gotten stuck in a local minimum, and you should make sure to start your next run with smaller gaussian widths for each parameter around the best fit. This will also show up on the triangle plot, which will be generated after you specify the burn-in. The best fit values for each parameter will be highlighted by green crosshairs. 

A "table" of median values plus/minus one sigma error bars will be saved, along with everything else, in the MCMCRUNS/vis_only/NUMWALKERSxNUMSTEPS/ directory. I say "table" because really it's only half a table, and you'll need to provide the rest of the formatting when you copy it into latex. It should save you some time though. Don't forget to deal with sigfigs! One sigma upper and lower limits are also provided in case the distribution is one tailed. See the triangle plot to dechiper whether or not this is the case. When you get to the point of wanting to report values from your MCMC runs, a good rule of thumb is that you want the number of steps you're formulating the triangle plot/tables off of (ie. NUMSTEPS - burn-in) to be ~10x the autocorrelation time for the reported uncertainties to be robust. Side-note: autocorrelation time goes down as you increase the number of walkers, but you could also just do more steps for a given number of walkers. Either will work. Once you've got a visibility-only fit you're happy with, we can move on to simultaneous visibility and SED fitting. 
