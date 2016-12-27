#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
from moviepy.video.io.bindings import mplfig_to_npimage

import sys
import numpy as nm
from pylab import *
from matplotlib import pyplot as plt
import pprint

import AthenaModel as ath

import os,sys

# from matplotlib.axes import subplot_class_factory

#from mpl_toolkits.axes_grid1.inset_locator import inset_axes, 

from mpl_toolkits.axes_grid1 import make_axes_locatable


# from read_para import *
# from hdfRadHydroGrid import *

def set_fonts_etc():
    fontsize_bar=14
    fontsize_x = 16
   
    ax.set_ylabel('z', fontsize = 22)
    for ylabel in ax.get_yticklabels():
        ylabel.set_fontsize(fontsize_x)
    
    ax.set_xlabel('R', fontsize = 22)
    for xlabel in ax.get_xticklabels():
        xlabel.set_fontsize(fontsize_x)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)                        

    cb =plt.colorbar(im, cax=cax)        
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(fontsize_bar)
    
Nx=0 
Nz=0

what2do = '2pan'
what2do = '1pan'


# put0='/Volumes/Apps_and_Docs/avdorodn/WORK/DATA_tor8/'
# runs = ['run1', 'n1e10_G01']


pathToFileToReadDir = '/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind/bin/'
print(pathToFileToReadDir)
os.chdir(pathToFileToReadDir)

files = []
dat = []

dat = ath.athDataModel()

f = plt.figure()    
ax = f.add_subplot(111)
ii=0
stp=1
phToSHow = 1
for fileInDir in os.listdir("."):
    if fileInDir.startswith("mhdXwind") & fileInDir.endswith("bin")  :
        Z=[]
        X=[]
        fileToRead = pathToFileToReadDir + fileInDir
        print(fileToRead)
        dat.loadDataFromBinFiles(fileToRead, dat, printDetail = True)
        mx = dat.Mx[dat.i_s:dat.ie:stp, phToSHow, dat.js:dat.je:stp]
        mt = dat.Mt[dat.i_s]    
        mz = dat.Mz[dat.i_s:dat.ie:stp, phToSHow, dat.js:dat.je:stp]
        (X, Z) = meshgrid(dat.x,dat.z)
        x1 = X[dat.i_s:dat.ie:stp, dat.js:dat.je:stp]
        x3 = Z[dat.i_s:dat.ie:stp, dat.js:dat.je:stp]

        plt.streamplot(x1, x3, mx, mz, color='k',linewidth=2)
#         ax.set_title(r'$\log(\rho)$', fontsize = 22)    
#         set_fonts_etc()
        fname = '_tmp%03d.png' %ii
        print( 'Saving frame', fname)
        
        f.savefig(fname)
        files.append(fname)

        ii+=1
        if ii==10:
            break        
  
# show() 
# exit()       





#--
os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc \
        -lavcopts vcodec=wmv2 -oac copy -o animation_torus9.mpg")
                

