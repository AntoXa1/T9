#!/local/data/atorus1/dora/Compilers/epd-7.3-1-rh5-x86_64(1)/bin/python


##!/Library/Frameworks/Python.framework/Versions/Current/bin/python


import socket

from docutils.parsers.rst.directives.misc import Replace
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from numpy import ndarray, zeros, array, size, meshgrid, flipud, floor, where, amin, argmin, int
from pylab import *
from scipy import *

import AthenaModel
import convertTimeStamp
import numpy as nm
from physics1 import *


def round_to_n(x, n):
    " Round x to n significant figures "
    return round(x, -int(nm.floor(nm.sign(x) * nm.log10(abs(x)))) + n)

def roundThenStringToLatexFormat(x, n=2):
    " Format x into nice Latex rounding to n"
    power = int(nm.log10(round_to_n(x, 0)))
    f_SF = round_to_n(x, n) * pow(10, -power)
    return r"${}\cdot 10^{}$".format(f_SF, power)


def set_fonts_etc(ax):
    fontsize_bar=14
    fontsize_x = 16
   
    ax.set_ylabel('z(pc)', fontsize = 22)
    for ylabel in ax.get_yticklabels():
        ylabel.set_fontsize(fontsize_x)
    
    ax.set_xlabel('R(pc)', fontsize = 22)
    for xlabel in ax.get_xticklabels():
        xlabel.set_fontsize(fontsize_x)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)                        

    cb =plt.colorbar(im, cax=cax)        
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(fontsize_bar)
        
def add_inner_title(ax, title, loc, size=None, **kwargs):
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke        
    if size is None:        
        size = dict(size=plt.rcParams['legend.fontsize'])    
    size = dict(size=20)    
    at = AnchoredText(title, loc=loc, prop=size,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)    
    ax.add_artist(at)    
    at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=3)])    
    return at




Nx=0 
Nz=0

what2Do = 'calcTorusMass'
what2Do = None
strmPlot = False

if socket.gethostname()=='atorus':
     putToDataDirs= '/local/data/atorus1/dora/PROJECTS/'
     locdirList = [ 'AthenaWind/bin/', 'AthenaWind_cln2/bin/']
     paramFile =  [  putToDataDirs + x.replace('/bin',"") +'/tst/cylindrical' for  x in locdirList ]
     print paramFile;
     put_out= '/local/data/atorus1/dora/PROJECTS/SCRIPTS/T9/'
     put_FIG = '/local/data/atorus1/dora/PROJECTS/SCRIPTS/T9/'
     dataFileList = [['mhdXwind.0050.bin', 'mhdXwind.0150.bin', 'mhdXwind.0444.bin'], \
                                ['mhdXwind.0050.bin', 'mhdXwind.0150.bin', 'mhdXwind.0444.bin']]
else:
     putToDataDirs= '/Users/dora/WORK/ECLIPSE_SPACE/torus9/DATA/DAT_for_figures/' 
     locdirList = [ 'SL/', 'HW/']
     paramFile =  [  putToDataDirs + x  for  x in locdirList ]                             
     put_out= '/Users/dora/WORK/ECLIPSE_SPACE/torus9'
     put_FIG = '/Users/dora/Documents/TEX/torus9/'
    
     dataFileList = [['mhdXwind.0050.bin', 'mhdXwind.0150.bin', 'mhdXwind.0340.bin'], \
                                ['mhdXwind.0050.bin', 'mhdXwind.0150.bin', 'mhdXwind.0340.bin']]

numRaw=  len(dataFileList[0][:])


  
fig = plt.figure(1, (10, 7))  
grid = AxesGrid(fig, 111, # similar to subplot(132)
                    nrows_ncols = (numRaw, 2),
                    axes_pad = 0.0,
#                     axes_pad=(0.45, 0.35),
                    share_all=True,
                    label_mode = "L",
                    add_all=True,       
                    cbar_pad = 0.,             
                    cbar_location = "right",
                    cbar_mode="single", #"each",
                    direction = "column",
                    )
i_grid =0
timeTeX=[]

for i_dirs in range(len(locdirList)):    
     dirFileToReadBase = putToDataDirs + locdirList[i_dirs]
     dat = AthenaModel.athDataModel()
     dat.loadSimulationParam(paramFile[i_dirs] + '/athinput.torus9_hydro_2D', print_res=True)
     
     
     for fileToOpen in dataFileList[i_dirs] [:]:
#         print(fileToOpen,timeNumeric,arrayTimeTeX)          
        dat.loadDataFromBinFiles(dirFileToReadBase +fileToOpen, dat, printDetail = True)
        print("MBH=", dat.Mbh,   'R0=', dat.Rsc,  'n0=', dat.n0,
               "t_0=", dat.tsc/YR, 'YRs' )
        
        MBH = dat.Mbh
        (X, Z) = meshgrid(dat.x,dat.z)
        var1 = zeros(dat.nx, dat.nz)                    
        offset= 4
        offset_i= offset+50
        ist = offset_i;  ie = dat.nz-offset_i;
        jst = offset;      je =dat.nx-offset
        stp=1
        xmin = dat.x[jst]; xmx = dat.x[je-1]
        zmin = dat.z[ist]; zmx = dat.z[ie-1]
        phToSHow = 1        
        x1 = X[ist:ie:stp, jst:je:stp]
        x3 = Z[ist:ie:stp, jst:je:stp] 
        
        mx = dat.Mx[ist:ie:stp, phToSHow, jst:je:stp]
        mt = dat.Mt[ist:ie:stp, phToSHow, jst:je:stp]    
        mz = dat.Mz[ist:ie:stp, phToSHow, jst:je:stp]
        
        vx = dat.Mx[ist:ie:stp, phToSHow, jst:je:stp]/dat.ro[ist:ie:stp, phToSHow, jst:je:stp]
        vz = dat.Mz[ist:ie:stp, phToSHow, jst:je:stp]/dat.ro[ist:ie:stp, phToSHow, jst:je:stp]
        
        if dat.method == 'MHD':
            Bx = dat.Bx[ist:ie:stp, phToSHow, jst:je:stp]    
            Bz = dat.Bz[ist:ie:stp, phToSHow, jst:je:stp]
            print size(Bx), size(Bz)
            
        shape = (dat.nz,  dat.nt,  dat.nx)
        var1 =  log10(dat.ro[ist:ie, phToSHow, jst:je]) 
        Em = 0.5* (dat.Bt[ist:ie, phToSHow, jst:je]**2+dat.Bx[ist:ie, phToSHow, jst:je]**2+\
                  dat.Bz[ist:ie, phToSHow, jst:je]**2)
                
        Empol = 0.5* (dat.Bx[ist:ie, phToSHow, jst:je]**2+dat.Bz[ist:ie, phToSHow, jst:je]**2)
        
        Emtor = 0.5* (dat.Bt[ist:ie, phToSHow, jst:je]**2)
        
        
        var1 =  log10(Emtor/Empol)         
#         var1 =  log10(Em)
#         grid[i_grid].streamplot(x1, x3, mx, mz, color='r', linewidth=2)                   
        
        vatToShow2D=var1
        stp=17

        # qp1 = grid[i_grid].quiver(X[ist:ie:stp, jst:je:stp], Z[ist:ie:stp, jst:je:stp], (vx[ist:ie:stp, jst:je:stp]), 
        #                 (vz[ist:ie:stp, jst:je:stp]), width=0.008, scale=2,                            
        # pivot='mid', color='black', 
        # units='x' , headwidth =5, headlength =7,
        # linewidths=(0.5,), edgecolors=('black'))
#         cmap=cm.jet,
        im = grid[i_grid].imshow(vatToShow2D , interpolation='bilinear', cmap=plt.get_cmap('nipy_spectral'),    \
                        extent=[xmin, xmx, zmin, zmx] )    

        grid[i_grid].set_xlabel ("R(pc)", fontsize=22)
        grid[i_grid].set_ylabel ("z(pc)", fontsize=22)                 
        
        
        i_grid+=1
        timeNumeric = float(fileToOpen.split('.')[1])*dat.dt_bin*dat.tsc/YR              
        timeTeX.append(roundThenStringToLatexFormat(timeNumeric))         


fig.suptitle(r'Log (Emag$_{\rm tor}$ /Emag$_{\rm pol}$)', y=0.95, fontsize=16)        
grid.cbar_axes[0].colorbar(im)
for cax in grid.cbar_axes: 
    cax.toggle_label(True)

for ax, im_title in zip(grid, timeTeX):
        t = add_inner_title(ax, im_title, loc=2)
        t.patch.set_ec("none")
        t.patch.set_alpha(0.5)


fileNameToSave = 'emagVsTimeSL_HW_G0_5_6panel'
# fig.savefig(put_FIG +fileNameToSave + ".pdf", format='pdf')
        
show()
exit()

