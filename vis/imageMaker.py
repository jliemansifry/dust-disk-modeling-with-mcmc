import os

diskName = raw_input("What's the name of this disk? ")
noiseBox = "region=box'(100,100,300,900)'"
def makeImage(disk):
	os.system('fits in='+disk+'spw1.uvf op=uvin out='+disk+'spw1.vis')
	os.system('fits in='+disk+'spw2.uvf op=uvin out='+disk+'spw2.vis')
	os.system('fits in='+disk+'spw3.uvf op=uvin out='+disk+'spw3.vis')
	os.system('uvcat vis='+disk+'spw1.vis,'+disk+'spw2.vis,'+disk+'spw3.vis out='+disk+'spw123.vis')
	os.system('fits in='+disk+'spw123.vis op=uvout out='+disk+'spw123.uvf')
	os.system('invert vis='+disk+'spw1.vis,'+disk+'spw2.vis,'+disk+'spw3.vis map='+disk+'spw123.mp beam='+disk+'spw123.bm imsize=1024 cell=.00625 options=mfs robust=.5')
	os.system('clean map='+disk+'spw123.mp beam='+disk+'spw123.bm out='+disk+'spw123.cl niters=100000 cutoff=.0001')
	os.system('restor map='+disk+'spw123.mp beam='+disk+'spw123.bm model='+disk+'spw123.cl out='+disk+'spw123.cm')
	os.system('cgdisp in='+disk+'spw123.cm,'+disk+'spw123.cm type=pix,con slev=a,4.5e-5 levs1=-3,3,6,9,12,15,18 labtyp=arcsec device=/xs beamtyp=b,l,4')
	os.system('imstat in='+disk+'spw123.cm region=box'+r"'(100,100,300,900)'")
	#print 'imstat in='+disk+'spw123.cm region=box(100,100,300,900)'
	#print noiseBox

makeImage(diskName)
