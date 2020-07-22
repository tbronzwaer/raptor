#! /usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

'''READING THE DATA'''

img_size = 160

# 2D array to contain image
imgI = [[0 for x in xrange(img_size)] for x in xrange(img_size)]
imgQ = [[0 for x in xrange(img_size)] for x in xrange(img_size)]
imgU = [[0 for x in xrange(img_size)] for x in xrange(img_size)]
imgV = [[0 for x in xrange(img_size)] for x in xrange(img_size)]

fluxI = 0.
fluxQ = 0.
fluxU = 0.
fluxV = 0.

f = open('output/img_data_2.300000e+11_IQUV.dat','r');

for line in f:
    columns=line.split()
    # Pixel coordinates
    x = float(columns[0])
    y = float(columns[1])
    
    # Set image pixel
    imgI[int(x)][int(y)] = float(columns[2])
    imgQ[int(x)][int(y)] = float(columns[3])
    imgU[int(x)][int(y)] = float(columns[4])
    imgV[int(x)][int(y)] = float(columns[5])
    
    fluxI += float(columns[2])
    fluxQ += float(columns[3])
    fluxU += float(columns[4])
    fluxV += float(columns[5])

f.close()

imgI = pl.transpose(imgI)
imgQ = pl.transpose(imgQ)
imgU = pl.transpose(imgU)
imgV = pl.transpose(imgV)


'''PLOT RESULTS'''
# Clear the plot
plt.clf()
plt.close('all')


halfrange = 15



fig, axs = plt.subplots(1,4,figsize=(30,10))

iii1 = axs[0].imshow(imgI, cmap = 'afmhot', origin='lower', interpolation = 'Nearest', extent=[-halfrange,halfrange,-halfrange,halfrange]) #bicubic interp
axs[0].set_title('Stokes I')
axs[0].set_xlabel(r'$\alpha$ ($GM/c^2$)')
axs[0].set_ylabel(r'$\beta$ ($GM/c^2$)')
divider1 = make_axes_locatable(axs[0])
cax1 = divider1.append_axes("right", size="5%", pad=0.1)
cb1=fig.colorbar(iii1, cax=cax1, format='%.1e')

iii2 = axs[1].imshow(imgQ, cmap = 'RdBu', origin='lower', interpolation = 'Nearest', extent=[-halfrange,halfrange,-halfrange,halfrange],vmin=-np.max(imgQ),vmax=np.max(imgQ)) #bicubic interp
axs[1].set_title('Stokes Q')
axs[1].set_xlabel(r'$\alpha$ ($GM/c^2$)')
divider2 = make_axes_locatable(axs[1])
cax2 = divider2.append_axes("right", size="5%", pad=0.1)
cb2=fig.colorbar(iii2, cax=cax2, format='%.1e')

iii3 = axs[2].imshow(imgU, cmap = 'RdBu', origin='lower', interpolation = 'Nearest', extent=[-halfrange,halfrange,-halfrange,halfrange],vmin=-np.max(imgU),vmax=np.max(imgU)) #bicubic interp
axs[2].set_title('Stokes U')
axs[2].set_xlabel(r'$\alpha$ ($GM/c^2$)')
divider3 = make_axes_locatable(axs[2])
cax3 = divider3.append_axes("right", size="5%", pad=0.1)
cb3=fig.colorbar(iii3, cax=cax3, format='%.1e')

iii4 = axs[3].imshow(imgV, cmap = 'RdBu', origin='lower', interpolation = 'Nearest', extent=[-halfrange,halfrange,-halfrange,halfrange],vmin=-np.max(imgV),vmax=np.max(imgV)) #bicubic interp
axs[3].set_title('Stokes V')
axs[3].set_xlabel(r'$\alpha$ ($GM/c^2$)')
divider4 = make_axes_locatable(axs[3])
cax4 = divider4.append_axes("right", size="5%", pad=0.1)
cb4=fig.colorbar(iii4, cax=cax4, format='%.1e')




plt.tight_layout()

print('flux I')
print(fluxI)
print('flux Q')
print(fluxQ)
print('flux U')
print(fluxU)
print('flux V')
print(fluxV)

pl.savefig('figures/IQUV_plot.pdf', transparent=False)

