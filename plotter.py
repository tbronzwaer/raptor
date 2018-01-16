#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import griddata
from scipy.interpolate import interp2d

'''READING THE DATA'''


# 2D array to contain image

for i in range (0,1):
	flux = 0.
	Nimg= i #*25 + 1000
	#header=np.loadtxt('output/img_data_0_2.300000e+11_60.00.dat',userows=1)

	f = open('output/img_data_0_2.300000e+11_60.00.dat')
	line = f.readline()
	header = line.split()

	print "Flux of", float(header[2]), "at frequency ", float(header[3])
	print "Munit of ", header[4], "viewing angle of ", header[5], "Rlow ", header[6], "Rhigh", header[7]

	img_size=int(header[0])
	img_size2=int(header[1])

	data=np.loadtxt('output/img_data_0_2.300000e+11_60.00.dat',skiprows=1)+1e-20

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
