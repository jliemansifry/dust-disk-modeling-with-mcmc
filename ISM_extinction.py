from numpy import *
from scipy.interpolate import interp1d

mwl, ext1, ext2 = loadtxt('Mathis90.dat', unpack=True)
interpol = interp1d(mwl, ext1)
aV= 0.1 # input the extinction value here

def extinct(lambd):
    return aV * .282 * interpol(lambd)

lambdObs, fluxObs, errorObs = loadtxt('HIP79516_phot_correctUncertainties.txt', usecols=(0,1,2), unpack=True) # photometry text file as input with wavelength of observations, flux, and error as the 3 columns

newFluxObs = fluxObs * 10**(.4*extinct(lambdObs))
newErrorObs = errorObs * 10**(.4*extinct(lambdObs))

obsDAT = column_stack((lambdObs,newFluxObs,newErrorObs))

savetxt('HIP79516_phot_red_correctUncertainties_test.txt', obsDAT, delimiter=" ", fmt='%s') #use this output file as the photometry in disk_MASTER.py
