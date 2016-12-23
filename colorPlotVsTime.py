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
#from enstaller.vendor.requests.packages.urllib3 import filepost


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

def formatFileName(timeNum):
    time = str(timeNum)
    strLen= len(time)

    if strLen==1:
        return('000'+time)
    elif strLen==2:
        return('00'+time)
    elif strLen==3:
        return('000'+time)
    else:
        return(time)
    
def timeToStrTeX(timeNumeric, t_0):
    lenOfFilePostfix = 4
    arrayTimeTeX = []
    id0 = []
    for time, i in zip(timeNumeric,xrange(len(timeNumeric))):    
        timePhys =   timeNumeric[i]*t_0/YR
        arrayTimeTeX.append(roundThenStringToLatexFormat(timePhys))              
        print('time in yrs=',  timePhys, arrayTimeTeX)
        timeStr= str.zfill(timeNumeric[i], lenOfFilePostfix)
        id0.append(timeStr)
        print timeStr
    return(arrayTimeTeX)

Nx=0 
Nz=0

what2Do = 'calcTorusMass'
what2Do = None
strmPlot = False

if socket.gethostname()=='atorus':
    locdirList = ['bin/']    
    
    putToDataDirs= '/local/data/atorus1/dora/PROJECTS/AthenaWind_cln3/'
    
    putToParamFile = putToDataDirs + 'tst/cylindrical/'

    put_out= '/local/data/atorus1/dora/PROJECTS/SCRIPTS/T9/'
    put_FIG = '/local/data/atorus1/dora/PROJECTS/SCRIPTS/T9/'

else:
     putToDataDirs= '/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind/'
     putToParamFile = '/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind/tst/cylindrical/'
     locdirList = ['bin/' ]
                             
     put_out= '/Users/dora/WORK/ECLIPSE_SPACE/torus9'
     put_FIG = '/Users/dora/Documents/TEX/torus9/'
    

filePostFix = ['40', '90', '368', '492']

#filePostFix = ['60', '90', '368', '736']

#filePostFix = ['60', '90', '368', '625']



filePostFix = [str.zfill(x,4) for x in filePostFix]    
timeNumeric =[int(x) for x in filePostFix]


dirFileToReadBase = putToDataDirs #+ locdirList[0]

print(dirFileToReadBase)         
 
dat = AthenaModel.athDataModel()

dat.loadSimulationParam(putToParamFile + 'athinput.torus9_hydro_2D', print_res=True)

print("MBH=", dat.Mbh,   'R0=', dat.Rsc,  'n0=', dat.n0,      
       "t_0=", dat.tsc/YR, 'YRs',  dat.tsc ) 

MBH = dat.Mbh

timeTeX= []
phToSHow = 1        
 
fig = plt.figure(1, (10, 10))  
numRaw=  len(timeNumeric)
numRaw=  2

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


fileNamePrefix ="mhdXwind."
for postFix,i  in zip(filePostFix, range(len(filePostFix))):
    
    fileToOpen = putToDataDirs +locdirList[0] +fileNamePrefix +postFix+'.bin'
    print fileToOpen
    dat.loadDataFromBinFiles(fileToOpen, dat, printDetail = True)
    
    
    var1 = zeros(dat.nx, dat.nz)                    
    (X, Z) = meshgrid(dat.x,dat.z)
    offset= 4
    offset_i= offset+0
    ist = offset_i;  ie = dat.nz-offset_i;
    jst = offset;      je =dat.nx-offset
    stp=1
    xmin = dat.x[jst]; xmx = dat.x[je-1]
    zmin = dat.z[ist]; zmx = dat.z[ie-1]    
    x1 = X[ist:ie:stp, jst:je:stp]
    x3 = Z[ist:ie:stp, jst:je:stp] 

    mx = dat.Mx[ist:ie:stp, phToSHow, jst:je:stp]
    mt = dat.Mt[ist:ie:stp, phToSHow, jst:je:stp ]    
    mz = dat.Mz[ist:ie:stp, phToSHow, jst:je:stp]
    
    Bx = dat.Bx[ist:ie:stp, phToSHow, jst:je:stp]    
    Bz = dat.Bz[ist:ie:stp, phToSHow, jst:je:stp]
    print size(Bx), size(Bz)
    
    var1 =  (dat.ro[ist:ie, phToSHow, jst:je])    
    vatToShow2D=log10(var1)
    
    im = grid[i].imshow(vatToShow2D , interpolation='bilinear',cmap=cm.jet, 
                        extent=[xmin, xmx, zmin, zmx] )
 
    grid[i].set_xlabel ("R(pc)", fontsize=22)
    grid[i].set_ylabel ("z(pc)", fontsize=22)
    

    # timeNumeric = float(fileToOpen.split('.')[1])*dat.dt_bin*dat.tsc/YR              
    timeNumeric = float(fileToOpen.split('.')[1])*dat.dt_bin              

    timeTeX.append("t="+roundThenStringToLatexFormat(timeNumeric)+"$t_0$")         


grid.cbar_axes[0].colorbar(im)
for cax in grid.cbar_axes: 
    cax.toggle_label(True)

for ax, im_title in zip(grid, timeTeX):
        t = add_inner_title(ax, im_title, loc=2)
        t.patch.set_ec("none")
        t.patch.set_alpha(0.5)

       
fileNameToSave = 'rhoSOL_G0.0_4panel'
fig.savefig(put_FIG +fileNameToSave + ".pdf", format='pdf')
        
show()
exit()

       
