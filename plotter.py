#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy import signal
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
from cython.parallel import parallel, prange
from distutils.core import setup
from distutils.extension import Extension
import cmocean

'''READING THE DATA'''

img_size  = 128
img_size2 = 128
# 2D array to contain image

for i in range (0,1):
	print i
	flux = 0.
	Nimg= i #*25 + 1000
	data=np.loadtxt('test.dat')+1e-20

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
	img = transpose(img)
	img = flipud(img)
	print np.max(img)

	'''PLOT RESULTS'''
# Clear the plot
	plt.clf()
	plt.close('all')
	# Set plot size
	plotsize = 127
	fig = plt.figure(figsize=(img_size/10,img_size/(10)),dpi=10)
	halfrange = 20
	ax = fig.add_subplot(111)


	im = plt.imshow(np.sqrt(img), interpolation = 'Nearest', cmap='magma', extent=[-halfrange,halfrange,-halfrange,halfrange])
	im.axes.get_yaxis().set_visible(False)
	ax.set_frame_on(False)
	ax.set_xticks([]); ax.set_yticks([])
	plt.axis('off')
	plt.axis('tight')

#	plt.clim(1e-3,1.)
	fig.set_size_inches(img_size2/10,img_size2/(10),forward=True)
	fig.tight_layout()
	plt.show()
	#savefig('figures/img_%d.png'%Nimg, bbox_inches='tight', transparent=False,pad_inches=0,dpi=99.28)
