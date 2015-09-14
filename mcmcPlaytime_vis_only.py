import matplotlib # necessary to save plots remotely; comment out if local
matplotlib.use('Agg') # comment out if local
import triangle
from scipy.stats import mstats
import numpy as np
import emcee
from visgen_master import VisibilityGenerator
from disk_spl import Disk
from astropy.io import fits
import os
import time 
from emcee.utils import MPIPool
import sys
go = ''
import matplotlib.pyplot as plt

append = 'vis_only' # something to add to filenames
whatbywhat= raw_input('MCMC what by what? ') # walkers x steps
disk_name= raw_input('Name of disk? ') #'HIP79516' 
noise = raw_input('RMS noise (in microJy)? ') #'HIP79516'  ##raw_input('What disk? ')
chainDat = fits.getdata('MCMCRUNS/vis_only/'+whatbywhat+'/'+whatbywhat+'.chain.fits')
chiDat = fits.getdata('MCMCRUNS/vis_only/'+whatbywhat+'/'+whatbywhat+'.chi.fits')

ndim = 5

blowout_size = (3 * 2.89  * 3.8e26)/ (16 * 3.14159 * 2.7 * 6.67e-11 * 1. *1.988e33 * 2.9979e8)# <<blowout size

def find_best_fit():
    print '\nMin chi: '+str(np.min(chiDat))
    bestFit = np.ravel(np.array(np.where(chiDat == np.min(chiDat))))
    global best_disk, vis_bestfit, diskBestFit
    print np.array(np.where(chiDat == np.min(chiDat))) #helpful to see where in the chain the min chi was
    best_disk = chainDat[bestFit[0],bestFit[len(bestFit)/2]] # picks a step where walker was at min chi
    vis_bestfit = np.ravel([best_disk[3],best_disk[4]])
    diskBestFit = Disk(0.0001,best_disk[0],best_disk[1],10**best_disk[2],0, 0.8, blowout_size * 1e6,1)

def vis_gen(visVal,disk,append):
    print 'Generating visibilities...'
    visGen = VisibilityGenerator(1024, visVal[0], visVal[1], whatbywhat+'_'+append)
    visGen.generateFITS(disk)
    visGen.generateVisibility()

def miriad_residual(whatbywhat,append,noise,disk_name):
    print 'Making residual image...'
    os.system('cp -rf images/' + whatbywhat +'* MCMCRUNS/vis_only/' +whatbywhat+'/')
    os.system('rm -rf images/' + whatbywhat + '*')
    os.system('uvmodel vis=vis/'+disk_name+'spw123.vis model=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat+'_'+append + '.mp options=subtract out=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.vis')
    os.system('invert vis=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.vis map=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.mp beam=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.bm imsize=1024 cell=.00625 options=mfs robust=.5')
    os.system('clean map=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.mp beam=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.bm out=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.cl niters=100000 cutoff=.0001')
    os.system('restor map=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.mp beam=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.bm model=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.cl out=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.cm')
    os.system('cgdisp in=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.cm,MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_Residual.cm type=pix,con slev=a,'+noise+'e-6 levs1=-3,3,6,9,12,15,18,21 labtyp=arcsec device=/xs beamtyp=b,l,4')

def model_maker(whatbywhat,append,noise,disk_name):
    print 'Making model...'
    os.system('invert vis=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'.vis,MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'2.vis,MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'3.vis map=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_image.mp beam=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_image.bm imsize=1024 cell=.00625 options=mfs robust=.5')
    os.system('clean map=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_image.mp beam=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_image.bm out=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_image.cl niters=100000 cutoff=.0001')
    os.system('restor map=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_image.mp beam=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_image.bm model=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_image.cl out=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_image.cm')
    os.system('cgdisp in=MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_image.cm,MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat + '_'+append+'_image.cm type=pix,con slev=a,'+noise+'e-6 levs1=-3,3,6,9,12,15,18 labtyp=arcsec device=/xs beamtyp=b,l,4')

def image_gen():
    find_best_fit()
    vis_gen(vis_bestfit,diskBestFit,'BestFit')
    prompt= raw_input('Type y to generate best residual image ')
    if prompt == 'y':
        miriad_residual(whatbywhat,'BestFit',noise,disk_name)
    prompt= raw_input('Type y to generate best model image ')
    if prompt == 'y':
        model_maker(whatbywhat,'BestFit',noise,disk_name)

def sed_gen():
    plt.clf()
    prompt= raw_input('Type y to generate best SED ')
    if prompt == 'y':
        diskBestFit.plotSED(disk_name,whatbywhat,append)

def plot_chi():
    x2 = range(len(chiDat[0,:]))
    for i in chiDat: 
        plt.plot(x2, i[:])
    plt.xlabel('Step #')
    plt.ylabel('Raw $\chi^2$')
    plt.ylim(np.min(chiDat),np.min(chiDat)+1000)
    plt.savefig('MCMCRUNS/vis_only/' + whatbywhat + '/' + whatbywhat+ '_'+append+'_burn_in')
    plt.show()

def triangle_gen():
    plt.clf()
    plot_chi()
    burnIn = input('Burn off how many steps? ')
    samples = chainDat[:,burnIn:,:].reshape((-1,ndim))
    triangle.corner(samples, labels=['$R_{In}$ [AU]','$\Delta$R [AU]','log($M_{D}$ [$M_{Earth}$])','i [degrees]','PA [degrees]'], truths=[best_disk[0],best_disk[1],best_disk[2],best_disk[3],best_disk[4]], truth_color='g', quantiles=[0.16,.5,.84], show_titles=True, plot_contours=True)
    plt.savefig('MCMCRUNS/'+append+'/'+whatbywhat+'/'+disk_name+'_'+whatbywhat+'_'+append+'_Triangle.png')
    plt.show()
    
    values = np.zeros((ndim,3))
    for i in range(ndim):
            chainSqueeze=chainDat[:,burnIn:,i].reshape(-1)
        quantiles = mstats.mquantiles(chainSqueeze,prob=[.16,.5,.84],axis=None)
        values[i,0] = quantiles[1] #median
        values[i,1] = quantiles[2]-quantiles[1] #plus one sigma
        values[i,2] = quantiles[0]-quantiles[1] #minus one sigma
    
    values_68 = np.zeros((ndim,2))
    for i in range(ndim):
        chainSqueeze=chainDat[:,burnIn:,i].reshape(-1)
        quantiles = mstats.mquantiles(chainSqueeze,prob=[.32,.68],axis=None)
        values_68[i,0] = quantiles[0]
        values_68[i,1] = quantiles[1]

    #now gen table
    a = r'$R_{In}$  [AU] &'+' {0:.2f}$^{{{1:+.2f}}}_{{{2:+.2f}}}$ & {3:.2f} \\'.format(values[0,0],values[0,1],values[0,2],best_disk[0])
    b = r'$\Delta R$ [AU] &'+' {0:.2f}$^{{{1:+.2f}}}_{{{2:+.2f}}}$ & {3:.2f} \\'.format(values[1,0],values[1,1],values[1,2],best_disk[1])
    c = r'log($M_{Disk}$ [$M_{Earth}$]) &'+' {0:.2f}$^{{{1:+.2f}}}_{{{2:+.2f}}}$ & {3:.2f} \\'.format(values[2,0],values[2,1],values[2,2],best_disk[2])
    d = r'$i$ [$^\circ$] &'+' {0:.2f}$^{{{1:+.2f}}}_{{{2:+.2f}}}$ & {3:.2f} \\'.format(values[3,0],values[3,1],values[3,2],best_disk[3])
    e = r'$PA$ [$^\circ$] &'+' {0:.2f}$^{{{1:+.2f}}}_{{{2:+.2f}}}$ & {3:.2f} \\'.format(values[4,0],values[4,1],values[4,2],best_disk[4])
    
    table_text = open('MCMCRUNS/'+append+'/'+whatbywhat+'/'+disk_name+'_'+whatbywhat+'_'+append+'_Table.txt','w')    
    table_text.write('Median plus/minus 1 sigma & Best fit value.\nBelow is easily copy-able into a latex file. Make sure you edit sig figs on your own!\n')
    table_text.write(a+'\n')
    table_text.write(b+'\n')
    table_text.write(c+'\n')
    table_text.write(d+'\n')
    table_text.write(e+'\n\n')
    table_text.write('One sigma values for single tailed distributions with only an upper or lower constraint. See triangle plot. \n 0.32 & 0.68 \n')
    for num in range(len(values_68)):
        table_text.write(str(values_68[num,0])+' & '+str(values_68[num,1])+'\n')
    table_text.close()

image_gen()
sed_gen()
prompt= raw_input('Type y to generate triangle ')
if prompt == 'y':
    triangle_gen()
