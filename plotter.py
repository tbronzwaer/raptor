#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import griddata
from scipy.interpolate import interp2d

filename='output/img_data_1_2.300000e+11_60.00.dat'
outputfile='output/img.png'

print "RAPTOR imaging visualization script\n"

print "Reading header...\n"

f = open(filename)
line = f.readline()
header = line.split()
print "IMAGE INFORMATION"
print "Image size", int(header[0]), "x", int(header[1])
print "Flux of", float(header[2]), "at frequency", float(header[3])

print "\nMODEL PARAMETERS"
print "Used grmhd file:", header[4]
print "Munit:", header[5]
print "Viewing angle:", header[6]
print "Rlow:", header[7]
print "Rhigh:", header[8]

print "\nReading data..."

img_size=int(header[0])
img_size2=int(header[1])

data=np.loadtxt(filename,skiprows=1)+1e-20

print "\nProcessing data..."

if(img_size!=img_size2):
	xi = np.linspace(0,img_size,img_size2)
	yi = np.linspace(0,img_size,img_size2)
	xi = np.tile(xi,img_size2)
	yi = np.repeat(yi,img_size2)
	x  = np.linspace(0,img_size,img_size)
	y  = np.linspace(0,img_size,img_size)
	x = np.tile(x,img_size)
	y = np.repeat(y,img_size)
	img = griddata((x,y),data,(xi,yi),method='cubic',rescale=True)
else:
	img=data

img = np.reshape(img,(-1,img_size2))
img = np.transpose(img)
img = np.flipud(img)

print "\nPlotting data..."

plt.clf()
plt.close('all')

fig = plt.figure(figsize=(img_size,img_size),dpi=1)
halfrange = 20
ax = fig.add_subplot(111)

im = plt.imshow(np.sqrt(img/np.max(img)), interpolation = 'Nearest', cmap='magma', extent=[-halfrange,halfrange,-halfrange,halfrange])
im.axes.get_yaxis().set_visible(False)
ax.set_frame_on(False)
ax.set_xticks([]); ax.set_yticks([])
plt.axis('off')
plt.axis('tight')

fig.set_size_inches(img_size,img_size,forward=True)
fig.tight_layout()
#plt.show()

plt.savefig(outputfile, bbox_inches='tight', transparent=False,pad_inches=0,dpi=1)

print "\nDone!"
