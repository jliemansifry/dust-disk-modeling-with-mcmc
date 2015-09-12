import math
import matplotlib # necessary to save plots remotely; comment out if local
matplotlib.use('Agg') # comment out if local
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pylab
import numpy
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
from scipy.special import gamma
import warnings
from astropy.io import fits
import time

class Disk:

    stefbolt_const = 5.67e-8
    h_const = 6.6260695729e-34
    k_const = 1.380648813e-23
    c_const = 2.99792458e8
    
    """
    inner and outer radius in units of AU
    grain size in units of microns
    disk mass in units of earth masses
    grain efficiency, ie. beta, unitless
    power law (positive values mean decreasing surface density with radius), unitless
    """
    def __init__(self, ininR, innerRadius, delR, diskMass, beltMass, grainEfficiency, grainSize, powerLaw):
        # set the grain density to 2.7 grams/cm^3 in kilograms/m^3
        self.grainDensity = 2.7e3
        # set the star temperature to solar temperature in Kelvin
        self.starTemperature = 6500.
        # set the star radius to solar radius in meters
        self.starRadius = 6.955e8*1.59
        # set the star luminosity to solar luminosity in Watts
        self.starLuminosity = 3.839e26*4.04
        # set the star distance in meters
        self.starDistance = 133.69*3.08568e16
        # disable overflow warnings
        warnings.filterwarnings("ignore", "overflow encountered in double_scalars")
        # convert the radii to meters
        self.ininR = float(ininR)*1.496e11 # inner edge of inner belt
        self.ininW = 10 * 1.496e11 # width of inner belt
        self.innerRadius = float(innerRadius)*1.496e11
        self.outerRadius = (float(innerRadius)+float(delR))*1.496e11
        # convert grain size to meters
        self.grainSize = float(grainSize)/1e6 #2.46/1e6
        # convert to kilograms
        self.diskMass = float(diskMass)*5.9742e24
        self.grainEfficiency = float(grainEfficiency)
        self.grainEfficiency2 = float(grainEfficiency)
        self.beltMass = float(beltMass)*5.9742e24
        # set power law for outer disk
        self.powerLaw = float(powerLaw)
        # calculate zeta(grainEfficency + 4)
        argument = self.grainEfficiency + 4
        self.zeta = 0
        for n in range(1, 100):
            self.zeta += 1.0/(n**argument)
        # calculate the surface density at 100 AU (sigma_100)
        j = 2 - self.powerLaw
        m = (self.outerRadius**j) - (self.innerRadius**j)
        m2 = self.ininR+self.ininW - self.ininR # inner belt power law = 1 (decreasing sigma with radius... could fit for a new powerlaw if desired)
        self.surfaceSigma = self.diskMass*j/(2*math.pi*((100*1.496e11)**self.powerLaw)*m)
        self.surfaceSigma2 = self.beltMass/(2*math.pi*100*1.496e11*m2)
        self.grainMass = self.grainDensity*(4.0/3.0*math.pi*(self.grainSize**3))
        self.grainMass2 = self.grainDensity*(4.0/3.0*math.pi*(self.grainSize**3))
        #self.grainMass2 = self.grainDensity*(4.0/3.0*math.pi*(self.asteroidGrain**3)) # could describe a different characteristic grain size for inner belt if desired
        
        print 'Reading in SED data...'
        # read and store data from the Kurucz-Lejeune model 
        f = open('Model_Spectrum.txt','r')
        radius = self.starRadius
        dist = self.starDistance
        self.data_lambda = []
        self.data_flux = []
        self.data_flux_raw = []
        self.data_lambda_raw = []
        line = f.readline()
        while line != '':
            info = line.split()
            lamma = float(info[0])*1e-9 
	    flux = float(info[2])*4*3.14159*(radius/dist)**2*1e23
	    flux_raw = flux
            flux *= self.c_const/lamma
            self.data_lambda.append(lamma)
            self.data_flux.append(flux)
            self.data_flux_raw.append(flux_raw)
            self.data_lambda_raw.append(lamma)
            line = f.readline()
        f.close()
        self.data_lambda.append(1e4*1e-6)
        self.data_flux.append(2.914e-13*4*3.14159*(radius/dist)**2*1e23*self.c_const/1e-2)
        self.data_lambda = self.convertToMicrons(self.data_lambda)
        self.data_lambda_raw = self.convertToMicrons(self.data_lambda_raw)

        # read and store observed fluxes (including IRS fluxes if available) after dereddening observed photometry
        f = open('HIP79516_phot_correctUncertainties.txt', 'r')
        self.sample_lambda = []
        self.sample_flux = []
        self.sample_error = []
        line = f.readline()
        while line != '':
            # read wavelength, flux, and error
            info = line.split()
            # convert microns to meters, but use microns
            lam = float(info[0])*1e-6
            self.sample_lambda.append(float(info[0]))
            # make flux in Janksys into Jansky*Hz by muliplying with nu
            flu = float(info[1])*self.c_const/lam
            self.sample_flux.append(flu)
            # do the same thing with error
            err = float(info[2])*self.c_const/lam
            self.sample_error.append(err)
            line = f.readline()
        f.close()
        
        #read and store IRAS fluxes, which will only be used for display purposes, not fitting (if available)
        # comment this section out if no IRS spectrum is present
        f = open('HIP79516_spectrum.txt', 'r')
        self.IRS_lambda = []
        self.IRS_flux = []
        self.IRS_error = []
        self.IRS_lambda_raw = []
        self.IRS_flux_raw = []
        self.IRS_error_raw = []

        line = f.readline()
        while line != '':
            info = line.split()
            lam = float(info[0])
            self.IRS_lambda.append(lam)
            self.IRS_lambda_raw.append(lam)
            flu = float(info[1])*self.c_const/(lam*1e-6) #Want nu*F_nu
            self.IRS_flux.append(flu)
            self.IRS_flux_raw.append(float(info[1]))
            err = float(info[2])/flu*self.c_const/(lam*1e-6)*100 # could be plotted with SED but is quite large; not plotted for cleanliness
            self.IRS_error.append(err)
            self.IRS_error_raw.append(float(info[2]))
            line = f.readline()
        f.close()
        print 'Done!'
        
        print 'Reading in tabulated integrals...'
        #Read in Q*B table and associated arrays.
        compiled_temp = [float(x) for x in fits.open('./dust/compiled_temperature.fits')[0].data]
        compiled_grain_sizes = [float(x) for x in fits.open('./dust/compiled_grain_sizes.fits')[0].data]
        compiled_QBintegrals = fits.open('./dust/compiled_QBintegrals.fits')[0].data
        compiled_QFintegrals = fits.open('./dust/compiled_QFintegrals.fits')[0].data
        
        #Sorting, since the arrays were generated in an arbitrary order.
        self.sorted_lambda = numpy.array([float(x) for x in fits.open('./dust/compiled_lambda.fits')[0].data])
        self.sorted_QBintegrals = numpy.array([y for (x,y) in sorted(zip(compiled_grain_sizes,compiled_QBintegrals))])
        self.sorted_QFintegrals = numpy.array([y for (x,y) in sorted(zip(compiled_grain_sizes,compiled_QFintegrals))])
        self.sorted_grain_sizes = numpy.array(sorted(compiled_grain_sizes))
        self.sorted_temp = numpy.array(compiled_temp)
        
        print 'Calculating T(r)...'

        #Interpolate between 2 closest grain sizes. Need a second QBintegral_list and all related variables if a second grain size is used.
        self.QB_interp_funcs = [interp1d(self.sorted_grain_sizes,self.sorted_QBintegrals[:,temp], bounds_error=False, fill_value = 0) for temp in range(len(self.sorted_temp))]
        self.QBintegral_list = [f(self.grainSize) for f in self.QB_interp_funcs]
        
    #Generate an estimate of T(r) based on the tabulated integrals.
        self.rad_steps = numpy.arange(self.innerRadius,self.outerRadius,1.496e10)
        self.rad_steps2 = numpy.arange(self.ininR,self.ininR+self.ininW,1.496e10)
        self.grain_temps = []
        self.grain_temps2 = []
        raw_integral = interp1d(self.sorted_grain_sizes, self.sorted_QFintegrals)(self.grainSize)
    
        for i in range(len(self.rad_steps)):
            #Scale integrals as 1/r^2.  The integrals were calculated with r = 50 AU.
            lhs = raw_integral*(50*1.496e11/self.rad_steps[i])**2/4/math.pi
            QBintegral_close = min(self.QBintegral_list, key=lambda y: math.fabs(y-lhs))
            QBintegral_index = numpy.where(self.QBintegral_list==QBintegral_close)[0][0]
            temperature = self.sorted_temp[QBintegral_index]
            self.grain_temps.append(temperature)
        if len(self.rad_steps) > 1:
            self.temp_function = interp1d(self.rad_steps, self.grain_temps, bounds_error=False, fill_value = self.grain_temps[len(self.grain_temps)-1])
        if len(self.rad_steps) == 1:
            def temp_function(self):
                return self.grain_temps[0]        
        
        for i in range(len(self.rad_steps2)):
            lhs2 = raw_integral*(50*1.496e11/self.rad_steps2[i])**2/4/math.pi
            QBintegral_close = min(self.QBintegral_list, key=lambda y: math.fabs(y-lhs2))
            QBintegral_index = numpy.where(self.QBintegral_list==QBintegral_close)[0][0]
            temperature2 = self.sorted_temp[QBintegral_index]
            self.grain_temps2.append(temperature2)
        if len(self.rad_steps2) > 1:
            self.temp_function2 = interp1d(self.rad_steps2, self.grain_temps2, bounds_error=False, fill_value = self.grain_temps2[len(self.grain_temps2)-1])
        if len(self.rad_steps2) == 1:
            def temp_function2(self):
                return self.grain_temps2[0]        
        print 'Done!'
       
        # generate interpolation function
        loglamb = map(math.log10, self.data_lambda)
        logflux = map(math.log10, self.data_flux)
        # if out of bounds, interpolates to 0
        logFunct = interp1d(loglamb, logflux, bounds_error=False, fill_value=0)
        self.interpol_funct = lambda x: 10**(logFunct(math.log10(x)))
    
        # Add in an asteroid belt with a fixed temperature and mass
    # this simple blackbody version of the asteroid belt isn't used in this code, but could be reinstated if desired. leaving these variables defined for simplicity
        self.asteroid_radius = 1.0e-6 #arbitrary
        self.M_aster = 4/3*math.pi*self.asteroid_radius**3*self.grainDensity 
        self.n_asteroids = self.beltMass/self.M_aster
        self.Temp_a = 0.001
        
    def getOuterRadius(self):
        """ gets the outer radius in AU """
        return self.outerRadius/1.496e11
    
    def getInnerRadius(self):
        """ gets the inner radius in AU """
        return self.innerRadius/1.496e11
    
    def getStarDistance(self):
        """ gets the star distance in parsecs """
        return self.starDistance/3.08568e16
    
    def convertToMicrons(self, lst):
        """ converts a list of meters into microns """
        return [x*1e6 for x in lst]
    
    def calculatePointFlux(self, radius, frequency):
        """ Takes a radius (in AU) and frequency (in GHz), returns the point flux at that radius and frequency in Janskys"""
        radius = float(radius)*1.496e11
        # make sure that we are inside of the ring
        if radius > self.outerRadius or radius < self.innerRadius:
            return 0
        # convert frequency to GHz
        lamma = self.c_const/(frequency*1e9)
        Q = self.qFunction(lamma)
        flux = 1e26*Q*(self.grainSize**2)*self.calculateGrainDistribution(radius)*self.calculateGrainBlackbody(radius, lamma)/(self.starDistance**2)
        return flux
   
    def calculatePointFlux2(self, radius, frequency):
        radius = float(radius)*1.496e11
        if radius > self.ininR+self.ininW or radius < self.ininR:
            return 0
        lamma = self.c_const/(frequency*1e9)
        Q2 = self.qFunction2(lamma)
        flux2 = 1e26*Q2*(self.grainSize**2)*self.calculateGrainDistribution2(radius)*self.calculateGrainBlackbody2(radius, lamma)/(self.starDistance**2)
        return flux2

    def calculateFlux(self, lamma):
        """ input lambda wavelength in meters ,return nu*B_nu(lambda) in Jansky*Hz """
        # integrate returns a list of integral value and error, we only want value
        fluxIntegral = integrate.quad(lambda radius: radius*self.calculateGrainDistribution(radius)
                                                 *self.calculateGrainBlackbody(radius, lamma), self.innerRadius, self.outerRadius,
                                                 )[0]
        nu = self.c_const/lamma
        Q = self.qFunction(lamma)
        flux = Q*nu*2*math.pi*fluxIntegral*1e26*(self.grainSize**2)/(self.starDistance**2)
        return flux
    
    def calculateFlux2(self, lamma):
        fluxIntegral2 = integrate.quad(lambda radius: radius*self.calculateGrainDistribution2(radius)
                                                 *self.calculateGrainBlackbody2(radius, lamma), self.ininR, self.ininR+self.ininW,
                                                 )[0]
        nu = self.c_const/lamma # scale by nu
        Q2 = self.qFunction2(lamma)
        flux2 = Q2*nu*2*math.pi*fluxIntegral2*1e26*(self.grainSize**2)/(self.starDistance**2)
        return flux2

    def calculateGrainDistribution(self, radius):
        """ returns the surface density at a given radius """
	if ((radius >= self.innerRadius) and (radius <=self.outerRadius)):
	    surfaceDensity = self.surfaceSigma*((radius/(100*1.496e11))**(-self.powerLaw))
        else:
	    surfaceDensity = 0
	return surfaceDensity/self.grainMass
    
    def calculateGrainDistribution2(self, radius):
	if ((radius >= self.ininR) and (radius <self.innerRadius)):
            surfaceDensity2 = self.surfaceSigma2*((radius/(100*1.496e11))**(-self.powerLaw))
	else:
	    surfaceDensity2 = 0
        return surfaceDensity2/self.grainMass2

    def calculateGrainBlackbody(self, radius, lamma):
        """ returns B_nu(T) in Janskys ( = 10^-26 * Watts / m^2 / Hz ) """
        try:
            nu = self.c_const/lamma
            exponent = self.h_const*nu/(self.k_const*self.calculateGrainTemperature(radius))
            numerator = 2*self.h_const*(nu**3)*math.pi 
            denominator = (self.c_const**2)*(math.e**exponent - 1)
            grainBlackbody = numerator/denominator
        except OverflowError:
            return 0
        return grainBlackbody
   

    def calculateGrainBlackbody2(self, radius, lamma):
        try:
            nu = self.c_const/lamma
            exponent = self.h_const*nu/(self.k_const*self.calculateGrainTemperature2(radius))
            numerator = 2*self.h_const*(nu**3)*math.pi 
            denominator = (self.c_const**2)*(math.e**exponent - 1)
            grainBlackbody2 = numerator/denominator
        except OverflowError:
            return 0
        return grainBlackbody2


    def calculatePlanckFunction(self, temperature, lamma):
        numerator = 2*self.h_const*(self.c_const**2)
        exponent = self.h_const*self.c_const/(lamma*self.k_const*temperature)
        denominator = (lamma**5)*(math.e**exponent-1)
        return numerator/denominator

    def calculateGrainTemperature(self, radius):
        """ Approximates the temperature of a grain using a precalculated table of integrals. """
        return self.temp_function(radius)
    def calculateGrainTemperature2(self, radius):
        return self.temp_function2(radius)
    
    def qFunction(self, lamma):
        """ Returns the emissivity of a grain at a given wavelength. """
        return 1-math.exp(-(2*math.pi*self.grainSize/lamma)**self.grainEfficiency)
    
    def qFunction2(self, lamma):
        """ necessary for a second grain size or beta """
        return 1-math.exp(-(2*math.pi*self.grainSize/lamma)**self.grainEfficiency2)
    
    def calculateStarBlackbody(self, lamma):
        """ returns nu*B_nu(lambda) in Jansky*Hz of the host star """
        nu = self.c_const/lamma
        exponent = self.h_const*nu/(self.k_const*self.starTemperature)
        numerator = 2*self.h_const*(nu**3)
        denominator = (self.c_const**2)*(math.exp(exponent)-1)
        nu = self.c_const/lamma
        # not sure what to scale by (4pi?)
        return numerator/denominator*nu*1e26*(self.starRadius**2)/(self.starDistance**2)
    
    def generateModel(self):
        """ generates lambda and nu*B_nu values in meters and Janksy*Hz, respectively """
        self.generateAsteroids()
        # sample from .1 microns to 10^4 microns
        x = numpy.arange(-7, -2, 0.01)
        x = [10**power for power in x]
        y = [self.calculateFlux(lamma) for lamma in x]
        z = [self.calculateFlux2(lamma) for lamma in x]
        self.disk_lambda = self.convertToMicrons(x)
        self.model_lambda = self.convertToMicrons(x)
        self.disk_flux = y
        self.belt_flux = z
        #z = self.asteroid_flux
        self.model_flux = map(lambda lamma,flux1,flux2: flux1+flux2+self.interpol_funct(lamma*1e6), x, y, z)
   
    def calculateAsteroidBelt(self, lamma):
        ''' for a blackbody belt that is NOT imaged-- not used in this version of the code'''
        try:
            nu = self.c_const/lamma
            exponent = self.h_const*nu/(self.k_const*self.Temp_a)
            numerator = 2*self.h_const*(nu**3)*math.pi*self.n_asteroids #Multiply by number of asteroids
            denominator = (self.c_const**2)*(math.e**exponent - 1)
            asteroidBlackbody = numerator/denominator*nu*1e26*self.asteroid_radius**2/self.starDistance**2
        except OverflowError:
            return 0
        return asteroidBlackbody
    
    def generateAsteroids(self):
        x = numpy.arange(-7, -2, 0.01)
        x = [10**power for power in x]
        y = [self.calculateAsteroidBelt(lamma) for lamma in x]
        self.asteroid_lambda = self.convertToMicrons(x)
        self.asteroid_flux = y
        #print 'Warm Belt Flux at 1.3mm = ', self.calculateAsteroidBelt(1.3e-3)*1.3e-3/self.c_const, 'Jy'

    def generateInterpolation(self):
        # generate interpolated data for the actual data we have
        def calculateInterpol(lam):
            # for less than 10^-6, ignore disk because so faint
            if lam*1e-6 > 1e-6:
                return self.calculateFlux(lam*1e-6)+self.interpol_funct(lam)+self.calculateFlux2(lam*1e-6) #self.calculateAsteroidBelt(lam*1e-6)
            else:
                return self.interpol_funct(lam)
        def calculateInterpolStar(lam):
            # for less than 10^-6, ignore disk because so faint
            if lam*1e-6 > 1e-6:
                return self.interpol_funct(lam)
            else:
                return self.interpol_funct(lam)
        self.interpol_flux = map(calculateInterpol, self.sample_lambda) #for avg IRS points
        self.interpolIRS_flux = map(calculateInterpol, self.IRS_lambda)
        self.interpolIRS_fluxStar = map(calculateInterpolStar, self.IRS_lambda)

    def plotSED(self, diskName, whatbywhat, append):
        """plots the SED of the star and disk """
        self.generateModel()
        plt.plot(self.data_lambda, self.Lsun(self.data_flux), '-.', label = 'Model Photosphere', linewidth=2, color='m')
        plt.loglog(self.model_lambda, self.Lsun(self.model_flux), '-', label = 'Best Fit Model', linewidth=2, color='b')
        plt.loglog(self.IRS_lambda, self.Lsun(self.IRS_flux), label = 'IRS Spectrum', linewidth=2, color='y') # comment out if no IRS spectrum
        plt.loglog(self.disk_lambda, self.Lsun(self.disk_flux), '--', label = 'Disk Model', linewidth=2, color='g')
        plt.loglog(self.disk_lambda, self.Lsun(self.belt_flux), ':', label = 'Inner Belt', linewidth=2, color='k')
        plt.errorbar(self.sample_lambda, self.Lsun(self.sample_flux), yerr=self.Lsun(self.sample_error), fmt='o', label = 'Observed Data', color='k')
        pylab.rcParams['xtick.major.pad']='5'
        pylab.rcParams['ytick.major.pad']='10'
        plt.tick_params(labelsize=16)
        plt.xlabel('$\lambda$ $(\mu m)$', fontsize=26)
        #plt.xlabel('$Wavelength $(\mu m)$', fontsize=20)
        plt.ylabel(r'$L_{\nu}$ $(L_{\odot})$', fontsize=26)
        plt.legend()
        plt.text(0.06, 50, diskName, fontsize=24)
        plt.xlim(5e-2,3e3)
        plt.ylim(1e-7,2e2)
        pylab.gcf().subplots_adjust(left=0.12, bottom=0.15)
        plt.savefig(diskName+'_'+whatbywhat+'_'+append+'_SED.png',dpi=400)
        #plt.show()

    
    def computeChiSquared(self):
        """ computes the chi-squared value of the disk and model data """
        self.generateInterpolation()
        lamma = self.sample_lambda
        model_flux = self.interpol_flux
        actual_flux = self.sample_flux
        error = self.sample_error
        chi_squared = 0
        numToFit = 0
        for i in range(len(lamma)):
            # the others are very close and we don't want them to mess up the chi-squared
            #if lamma[i] > 1e1 and lamma[i] < 1e3:
            if lamma[i] > 5.5 and lamma[i] < 1e4:
                numToFit += 1
        #print 'Wavelength: ',lamma[i]
                #print 'model:', model_flux[i]
        #print 'actual:', actual_flux[i]
                chi_squared += ((model_flux[i]-actual_flux[i])/error[i])**2
        # degrees of freedom = [number of samples] - [number parameters to fit]
        # TEMPO: parameters 6 normally
        dof = numToFit - 6
        chi_squared = chi_squared#/dof
        return chi_squared
    
    def Lsun(self, fluxes):
        ''' Converts fluxes that are in nu*F_nu into Solar Luminosity units '''
        SolLum = [4*math.pi*self.starDistance**2*i/3.839e26/1e26 for i in fluxes]
        return SolLum
    
    def diskMassestimate(self):
        ''' Use Beckwith's formula to get a realistic disk mass '''
        nu = 235.5e9
        kappa_nu = 0.1*(nu/1e12)**self.grainEfficiency
        Mass = 0.1*(0.02/kappa_nu) #Solar Masses
        return Mass/3.00266478958e-06 #Conversion to Earth Masses
    
    def disktemp(self):
        ''' Gets the approximate temperature of the disk '''
        x = numpy.arange(-7, -2, .01)
        x = [10**power for power in x]
        y = [self.calculateFlux(lamma) for lamma in x]
        max = [0,0]
        for i in range(0,len(y)):
            if y[i] > max[1]:
                max[0] = x[i]
                max[1] = y[i]
        lambda_peak = max[0]
        T = 2.898e-3/lambda_peak
        return T
                
    def irsSubtract(self):
        ''' subtract the IRS spectrum from the model to see possible silicate features '''
        self.generateInterpolation()
        plt.clf()
        model_flux = [flux * (lamma*1e-6)/self.c_const for flux,lamma in zip(self.interpolIRS_flux,self.IRS_lambda_raw)]
        star_flux = [flux * (lamma*1e-6)/self.c_const for flux,lamma in zip(self.interpolIRS_fluxStar,self.IRS_lambda_raw)]
        actual_flux = self.IRS_flux_raw ## in Jy
        error = self.IRS_error_raw # not plotted, it overwhelms the plot
        IRS_diff = [m_flux - s_flux for m_flux,s_flux in zip(model_flux,star_flux)]
        IRS_diffStar = [a_flux - s_flux for a_flux,s_flux in zip(actual_flux,star_flux)]
        plt.plot(self.IRS_lambda_raw,IRS_diff,color='g')
        plt.scatter(self.IRS_lambda_raw,IRS_diffStar,marker='.',color='b')
        plt.axhline(y=0,xmin=5,xmax=35,linestyle='-.',linewidth=2,color='m')
        plt.xlim(5,15)
        plt.ylim(-.1,.1)
        plt.xlabel('Wavelength $(\mu m)$')
        plt.ylabel('Excess Flux (Jy)')
        #plt.show()
        plt.savefig('test.png')

    def irsWindow(self):
        ''' plot the window around the IRS spectrum including the disk model '''
        self.generateInterpolation()
        plt.clf()
        model_flux = [flux * (lamma*1e-6)/self.c_const for flux,lamma in zip(self.interpolIRS_flux,self.IRS_lambda_raw)]
        actual_flux = self.IRS_flux_raw ## in Jy
        error = self.IRS_error_raw # not plotted, it overwhelms the plot
        axes = plt.loglog(self.IRS_lambda_raw, model_flux, label = 'Best Fit Disk', color='g')
        plt.scatter(self.IRS_lambda_raw, actual_flux,marker='.', label = 'IRS Spectrum', color='b')
        plt.loglog(self.data_lambda_raw, self.data_flux_raw,linestyle='-.', label = 'Model Photosphere', linewidth=2, color='m')
        plt.xlim(4.8,35)
        plt.ylim(1e-2,1e1)    
        plt.xlabel('Wavelength $(\mu m)$')
        plt.ylabel('Flux (Jy)')
        plt.legend()
        plt.savefig('testwindow.png')
        #plt.show()

    def plotGrainDistribution(self):
        radius = numpy.arange(1,self.outerRadius/1.496e11,.1)
	GDarr = numpy.zeros(len(radius))
	for i in range(len(radius)):
	    GDarr[i] = self.calculateGrainDistribution(radius[i]*1.496e11)*self.grainMass+ self.calculateGrainDistribution2(radius[i]*1.496e11)*self.grainMass2
        plt.loglog(radius,GDarr)
	plt.xlabel('Radius (AU)')
	plt.ylabel('Surface Density (g/m^2)')
	plt.savefig('test.png')
	plt.show()

    def testplot(self):
        '''For debugging'''
        xaxis = [1.5e11*x for x in numpy.arange(40,80,1)]
        yaxis = [self.temp_function(x) + self.temp_function2(x) for x in xaxis]
        plt.plot(xaxis, yaxis)
        #plt.plot(self.sorted_temp,self.integral_list)
        #plt.plot(self.sorted_temp,self.integral_list2,'g')
        plt.savefig('testplot.png')
        #plt.show()
