import math
import numpy
from astropy.io import fits
import os
from disk_spl import Disk

"""
For the generator to work, have folders named "fits" and "images" in
the directory from which you are calling this script. There should
also be a folder named "vis" in which the .vis and .uvf files 
for each spectral window for the disk being modeled are stored. 
"""
class VisibilityGenerator:
    
    """
    creates a visibility generator with width pixels,
    angle inclination and rotation rotation (in degrees)
    """
    def __init__(self, width, angle, rotation, fits_file):
        self.width = width
        # convert to radians
        self.incline = angle*math.pi/180.0
        self.rotation = (rotation+90)*math.pi/180.0
        
        # name files we need
	self.fits_name = fits_file
        self.fits_model_image = 'images/'+fits_file+'.fits'
        self.mir_model_image = 'images/'+fits_file+'.mp'
        self.mir_model_vis = 'images/'+fits_file+'.vis'
        self.fits_model_vis = 'images/'+fits_file+'.uvf'
        self.mir_model_vis2 = 'images/'+fits_file+'2.vis'
        self.fits_model_vis2 = 'images/'+fits_file+'2.uvf'
        self.mir_model_vis3 = 'images/'+fits_file+'3.vis'
        self.fits_model_vis3 = 'images/'+fits_file+'3.uvf'
        # get the data files we need
        self.mir_data_vis = 'vis/HIP79516spw1.vis'
        self.fits_data_vis = 'vis/HIP79516spw1.uvf'
        self.mir_data_vis2 = 'vis/HIP79516spw2.vis'
        self.fits_data_vis2 = 'vis/HIP79516spw2.uvf'
        self.mir_data_vis3 = 'vis/HIP79516spw3.vis'
        self.fits_data_vis3 = 'vis/HIP79516spw3.uvf'
    ## time averaged, weight corrected, channel averaged data in XXYY^^

        # get the data from the observed visibility
        if not os.path.exists(self.fits_data_vis):
            os.system('fits in=' + self.mir_data_vis + ' op=uvout out=' + self.fits_data_vis)
        # store data for computing the chi squared later on
        observed = fits.getdata(self.fits_data_vis, 'PRIMARY').data[:,0,0,0,0]
        observed2 = fits.getdata(self.fits_data_vis2, 'PRIMARY').data[:,0,0,0,0]
        observed3 = fits.getdata(self.fits_data_vis3, 'PRIMARY').data[:,0,0,0,0]
        
	#STOKES I CALCULATED HERE
        numbvis1 = len(observed[:,0,0])
        for i in range(0,numbvis1):
            loc = numpy.array(numpy.where(observed[i,:,2] !=0))
            realI = ((numpy.mean(observed[i,0,0])) + (numpy.mean(observed[i,1,0])))/2
            imagI = ((numpy.mean(observed[i,0,1])) + (numpy.mean(observed[i,1,1])))/2
            realI = (observed[i,0,0] + observed[i,1,0])/2.0
            imagI = (observed[i,0,1] + observed[i,1,1])/2.0
            if (observed[i,0,2] != 0. and observed[i,1,2] != 0.):
                w_I = (4.0 / (1./observed[i,0,2] + 1./observed[i,1,2]))
            else:
                w_I = 0.0
            observed[i,0,0] = realI
            observed[i,0,1] = imagI
            observed[i,0,2] = w_I

        self.data_vis_real = observed[:,0,0]
        self.data_vis_imag = observed[:,0,1]
        self.data_vis_weight = observed[:,0,2]
        self.data_vis_real[numpy.isnan(self.data_vis_real)] = 0
        self.data_vis_imag[numpy.isnan(self.data_vis_imag)] = 0
        self.data_vis_weight[numpy.isnan(self.data_vis_weight)] = 0
    
        numbvis2 = len(observed2[:,0,0])
        for i in range(0,numbvis1):
            loc = numpy.array(numpy.where(observed2[i,:,2] !=0))
            realI = ((numpy.mean(observed2[i,0,0])) + (numpy.mean(observed2[i,1,0])))/2
            imagI = ((numpy.mean(observed2[i,0,1])) + (numpy.mean(observed2[i,1,1])))/2
            realI = (observed2[i,0,0] + observed2[i,1,0])/2.0
            imagI = (observed2[i,0,1] + observed2[i,1,1])/2.0
            if (observed2[i,0,2] != 0. and observed2[i,1,2] != 0.):
                w_I = (4.0 / (1./observed2[i,0,2] + 1./observed2[i,1,2]))
            else:
                w_I = 0.0
            observed2[i,0,0] = realI
            observed2[i,0,1] = imagI
            observed2[i,0,2] = w_I

        self.data_vis_real2 = observed2[:,0,0]
        self.data_vis_imag2 = observed2[:,0,1]
        self.data_vis_weight2 = observed2[:,0,2]
        self.data_vis_real2[numpy.isnan(self.data_vis_real2)] = 0
        self.data_vis_imag2[numpy.isnan(self.data_vis_imag2)] = 0
        self.data_vis_weight2[numpy.isnan(self.data_vis_weight2)] = 0
    
        numbvis3 = len(observed3[:,0,0])
        for i in range(0,numbvis3):
            loc = numpy.array(numpy.where(observed3[i,:,2] !=0))
            realI = ((numpy.mean(observed3[i,0,0])) + (numpy.mean(observed3[i,1,0])))/2
            imagI = ((numpy.mean(observed3[i,0,1])) + (numpy.mean(observed3[i,1,1])))/2
            realI = (observed3[i,0,0] + observed3[i,1,0])/2.0
            imagI = (observed3[i,0,1] + observed3[i,1,1])/2.0
            if (observed3[i,0,2] != 0. and observed3[i,1,2] != 0.):
                w_I = (4.0 / (1./observed3[i,0,2] + 1./observed3[i,1,2]))
            else:
                w_I = 0.0
            observed3[i,0,0] = realI
            observed3[i,0,1] = imagI
            observed3[i,0,2] = w_I

        self.data_vis_real3 = observed3[:,0,0]
        self.data_vis_imag3 = observed3[:,0,1]
        self.data_vis_weight3 = observed3[:,0,2]
        self.data_vis_real3[numpy.isnan(self.data_vis_real)] = 0
        self.data_vis_imag3[numpy.isnan(self.data_vis_imag)] = 0
        self.data_vis_weight3[numpy.isnan(self.data_vis_weight)] = 0

    """
    generate a FITS file with name fits_file
    """
    def generateFITS(self, disk):
        # remove all current files which we will be working with
        os.system('rm -rf images/*')
        #remove if the fits file already exists
        if os.path.exists('fits/' + self.fits_name):
            os.system('rm -rf fits/' + self.fits_name)
            #print 'fits/' + self.fits_name, 'was overwritten.'
        
        # square images
        center = self.width/2.0
        
        # in units of AU
        outerRadius = disk.getOuterRadius()
        
        # now find the scaling factor
        pxwid = 3.0*outerRadius/self.width
        total_flux = 0
        
        # x and y are now transformed into new coordinate plane
        # center is at (0,0) with points sampled at center of each pixel
        
        # initialize an empty 2-D "image" of 0's
        image = [[0]*self.width for x in range(self.width)]
        # exploit the diagonal symmetry for a sub-2x speedup
        for j in xrange(self.width/2):
            y_face = (j-center+0.5)*pxwid

            for i in xrange(self.width):
                x_face = (i-center+0.5)*pxwid
                
                # first rotate, then transform for inclination angle
                x = x_face*math.cos(self.rotation) + y_face*math.sin(self.rotation)
                y = (-x_face*math.sin(self.rotation) + y_face*math.cos(self.rotation)) / math.cos(self.incline)
                
                radius = math.sqrt(x**2+y**2)
                # have flux in Jy/m^2, now need in Jy/pixel
                flux = (disk.calculatePointFlux(radius, 230.538)+disk.calculatePointFlux2(radius, 230.538))/math.cos(self.incline)*((pxwid*1.496e11)**2)
                if flux:
                    # only assign if the flux is NOT 0
                    image[j][i] = flux
                    image[-j-1][-i-1] = flux
                    total_flux += flux*2
        # convert the Python image array into a numpy 4-D array for MIRIAD compatibility
        #print "Total Flux =", total_flux
        image = numpy.array([[image]])
        
        # calculate data necessary for header
        stardist = disk.getStarDistance()
        pix_parsec = pxwid/stardist/3600.0 # AU to degrees
        
        # create hdu object to encapsulate data with the header
        hdu = fits.PrimaryHDU(image)
        head = hdu.header
        
        head.update('OBJECT', 'HD_114082')
        head.update('SIMPLE', 'T')
        head.update('NAXIS', 4)
        head.update('NAXIS1', self.width)
        head.update('NAXIS2', self.width)
        head.update('NAXIS3', 1)
        head.update('NAXIS4', 1)
        head.update('CDELT1', -1.0*pix_parsec)
        head.update('CRPIX1', center+0.5)
        head.update('CRVAL1', (16.0 + 13/60.0 + 34.306189/3600.0)*15)
        head.update('CTYPE1', 'RA---SIN')
        head.update('CDELT2', pix_parsec)
        head.update('CRPIX2', center+0.5)
        head.update('CRVAL2', -1.0*(45.0 + 49/60.0 + 04.06928/3600.0))
        head.update('CTYPE2', 'DEC--SIN')
        head.update('EPOCH', 2000)
        head.update('BSCALE', 1.0)
        head.update('BZERO', 0.0)
        head.update('CDELT3', 1.0)
        head.update('CRVAL3', -1.0)
        head.update('CRPIX3', 1.0)
        head.update('CTYPE3', 'STOKES')
        head.update('CDELT4', 5.625e9)
        head.update('CRPIX4', 1.0)
        head.update('CRVAL4', 240.538e9)
        head.update('CTYPE4', 'FREQ-LSR')
        head.update('BUNIT', 'JY/PIXEL')
        head.update('BTYPE', 'INTENSITY')

        # write to file
        hdu.writeto(self.fits_model_image)
        
        """
    print 'AU^2/pixel:', str(round(pxwid**2,4))
        print 'total flux (Jy):', total_flux
        print 'total flux (Jy*MHz):', str(round(total_flux*225.5e9/1e6,2))
        print 'actual SED flux at 225.5GHz:', str(round(disk.calculateFlux(3.0e8/225.5e9)/1e6,2))
        """
        
    def generateVisibility(self):
        """ generate the visibilities of the model using MIRIAD routines """
        # make sure we have the files we need
        if not os.path.exists(self.mir_data_vis):
            print self.mir_data_vis, "visibility file does not exist."
            exit()
       
        # convert FITS model image into a Miriad image
        os.system('fits in=' + self.fits_model_image + ' op=xyin out=' + self.mir_model_image)
        os.system('uvmodel model=' + self.mir_model_image + ' vis=' + self.mir_data_vis
                  + ' options=replace out=' + self.mir_model_vis)
        # convert the Miriad visibilities into FITS visibilities
        os.system('fits in=' + self.mir_model_vis + ' op=uvout out=' + self.fits_model_vis)
        os.system('uvmodel model=' + self.mir_model_image + ' vis=' + self.mir_data_vis2
                  + ' options=replace out=' + self.mir_model_vis2)
        os.system('fits in=' + self.mir_model_vis2 + ' op=uvout out=' + self.fits_model_vis2)
        os.system('uvmodel model=' + self.mir_model_image + ' vis=' + self.mir_data_vis3
                  + ' options=replace out=' + self.mir_model_vis3)
        os.system('fits in=' + self.mir_model_vis3 + ' op=uvout out=' + self.fits_model_vis3)
    
    def computeChiSquared(self, disk):
        # do all the prereqs
        self.generateFITS(disk)
        self.generateVisibility()
        # get the data from created files
        if not os.path.exists(self.fits_model_vis):
            return numpy.inf # addresses an occasional bug with tempfile on MCMC runs
        if not os.path.exists(self.fits_model_vis2):
            return numpy.inf
        if not os.path.exists(self.fits_model_vis3):
            return numpy.inf
        model_hdu = fits.open(self.fits_model_vis)
        model = model_hdu[0].data.field(5)
        model_hdu2 = fits.open(self.fits_model_vis2)
        model2 = model_hdu2[0].data.field(5)
        model_hdu3 = fits.open(self.fits_model_vis3)
        model3 = model_hdu3[0].data.field(5)
        # compute the chi squared
        # = sum of weight*[(model-data)**2] for real and imaginary
        #numVis = len(model)
        # use numpy vector operations
        real = model[::2,0,0,0,0,0]
        imag = model[::2,0,0,0,0,1]
        real2 = model2[::2,0,0,0,0,0]
        imag2 = model2[::2,0,0,0,0,1]
        real3 = model3[::2,0,0,0,0,0]
        imag3 = model3[::2,0,0,0,0,1]
        chi = float((((real-self.data_vis_real)**2)*self.data_vis_weight).sum())
        chi += float((((imag-self.data_vis_imag)**2)*self.data_vis_weight).sum())
        chi += float((((real2-self.data_vis_real2)**2)*self.data_vis_weight2).sum())
        chi += float((((imag2-self.data_vis_imag2)**2)*self.data_vis_weight2).sum())
        chi += float((((real3-self.data_vis_real3)**2)*self.data_vis_weight3).sum())
        chi += float((((imag3-self.data_vis_imag3)**2)*self.data_vis_weight3).sum())
        print chi
        # degrees of freedom = 2*[number of visibilities] - [number parameters to fit]
        # TEMPO: 4 parameters instead of 6 to fit
        #dof = 2*numVis-4
        chi = chi #/dof -- we need raw chi for MCMC
        model_hdu.close(closed=True)
        return chi
