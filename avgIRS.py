import math
import numpy as np
import matplotlib # necessary to save plots remotely; comment out if local
matplotlib.use('Agg') # comment out if local
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
from scipy.special import gamma
import warnings
import pyfits
import time

f = open('HIP79516_spectrum.txt', 'r') # the full IRS spectrum
IRS_lambda = []
IRS_flux = []
IRS_error = []
line = f.readline()
while line != '':
    info = line.split()
    lam = float(info[0])
    IRS_lambda.append(lam)
    IRS_flux.append(float(info[1]))
    IRS_error.append(float(info[2]))
    line = f.readline()
log_IRS_lambda = np.log10(IRS_lambda)
log_IRS_flux = np.log10(IRS_flux)
plt.errorbar(10**log_IRS_lambda, 10**log_IRS_flux, yerr=IRS_error, color = 'b', alpha = .5, label='IRS spectrum \n(only absolute errors plotted)')

log_IRS_avgLambda = []
log_IRS_avgFlux = []
log_IRS_avgError = []

ranges = [0,40,80,120,160,200,240,280,320,363] # ranges of the indexes you want to average over... generally between 360 and 365 points, make last index correspond to last data point in spectrum

for i in range(9):
    log_IRS_avgLambda.append(np.mean(log_IRS_lambda[ranges[i]:ranges[i+1]]))
    log_IRS_avgFlux.append(np.mean(log_IRS_flux[ranges[i]:ranges[i+1]]))

IRS_avgLambda = np.power(10,log_IRS_avgLambda)
IRS_avgFlux = np.power(10,log_IRS_avgFlux)

sumVar = np.zeros(9)
fluxVar = np.zeros(9)

for i in range(9):
    for j in range(ranges[i],ranges[i+1]):
        sumVar[i] += (IRS_flux[j]/40) * IRS_error[j]**2
        fluxVar[i] = (2.3025850929940 * IRS_avgFlux[i])**2 * sumVar[i]
fluxStd = (fluxVar**2 +(0.1*IRS_avgFlux)**2)**(.5) # assuming 10% IRS systematic flux uncertainty! (added in quadrature)

f2 = open('HIP79516_IRSavg.txt' ,'w')

for i in range(9):
    f2.write(str(IRS_avgLambda[i])+' '+str(IRS_avgFlux[i])+' '+str(fluxStd[i])+'\n')
f2.close()
plt.errorbar(IRS_avgLambda, IRS_avgFlux, yerr=fluxStd,color='r', label = 'Averaged IRS points \n(assumed 10% systematic flux uncertainty)')
plt.plot(IRS_avgLambda, IRS_avgFlux,color='r')
plt.xlabel('Wavelength (microns)')
plt.ylabel('Flux (Jy)')
plt.legend()
plt.savefig('HIP79516_IRS_spectrum_avg.png')
#plt.show()
f.close()
