#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

from scipy import *
import numpy as nm
from numpy import ndarray, zeros,array,size,meshgrid,flipud,floor,where,amin,argmin,int
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import AxesGrid,ImageGrid
from physics1 import *
import AthenaModel
import convertTimeStamp
from docutils.parsers.rst.directives.misc import Replace

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
    


Nx=0 
Nz=0

what2Do = 'calcTorusMass'
what2Do = None
strmPlot = False


put0= '/Users/dora/WORK/ECLIPSE_SPACE/torus9/DATA/'
put_out= '/Users/dora/Documents/TEX/torus9/'

locdir = 'SolovievSep201615_256x8x256_L0.n10e10'
put0 +=locdir



dat = AthenaModel.athDataModel()
dat.loadSimulationParam(put0 + '/athinput.torus9_hydro_2D', print_res=True)

 

timeNumeric =[120]

MBH = dat.Mbh


id0 = []
arrayTimeTeX = []

for time, i in zip(timeNumeric,xrange(len(timeNumeric))):   
          
#     arrayTimeTeX.append(roundThenStringToLatexFormat(timePhys))
    arrayTimeTeX.append(roundThenStringToLatexFormat(time))
                  
#     print('time in yrs=',  timePhys, arrayTimeTeX)
    print('time in yrs=',  arrayTimeTeX)
    
    timeStr= formatFileName(timeNumeric[i])
    
    id0.append(timeStr)
    
    print timeStr

fileNamePrefix ="mhdXwind."
fileToOpen = put0 +'/'+fileNamePrefix+id0[0]+'.bin'
# print(i, ' ', fileToOpen)          
fileToOpen='/Users/dora/WORK/ECLIPSE_SPACE/torus9/DATA/SolovievSep201615_256x8x256_L0.n10e10/mhdXwind.0120.bin'

dat.loadDataFromBinFiles(fileToOpen, dat, printDetail = True)


print("MBH=", dat.Mbh,   'R0=', dat.Rsc,  'n0=', dat.n0,
               "t_0=", dat.tsc/YR, 'YRs' )

    
timePhys=   float( timeNumeric[0] )*dat.tsc*dat.dt_bin/YR
print("timePhys = ", timePhys)                  
# exit()
#**************************



if what2Do == 'calcTorusMass':
    mTorus= dat.torusMass()
    print('torusMass=', mTorus/MSUN); exit()
    
    Nz = dat.dd[:,10].size
    plt.plot(dat.x, dat.Usc/1e5*dat.u1[Nz-2,:])
    show() 
    exit()


fig = plt.figure(1, (10, 10))  
grid = AxesGrid(fig, 111, # similar to subplot(132)
                    nrows_ncols = (2, 2),
#                     axes_pad = 0.0,
                    axes_pad=(0.45, 0.35),
                    share_all=True,
                    label_mode = "L",
                    add_all=True,       
                    cbar_pad = 0.,             
                    cbar_location = "right",
                    cbar_mode="each",
                    )
    

(X, Z) = meshgrid(dat.x,dat.z)
var1 = zeros(dat.nx, dat.nz)    


offset= 4
offset_i= offset+50
ist = offset_i;  ie = dat.nz-offset_i;
jst = offset;      je =dat.nx-offset
stp=1
xmin = dat.x[jst]; xmx = dat.x[je-1]
zmin = dat.z[ist]; zmx = dat.z[ie-1]

(X1, Z1) = meshgrid(dat.x[jst:je], dat.z[ist:ie])

Lk = zeros(dat.nx, dat.nz)
Lk = 1./sqrt(Z1**2 + X1**2)

# eps = 1.e-2;
# dat.ro[dat.ro<eps] = eps 

for i in range(4):
            
#     xmin = dat.x[jst]; xmx = dat.x[je-1]
#     zmin = -xmx/2.; zmx = xmx/2
     
    phToSHow = 1
    
    x1 = X[ist:ie:stp, jst:je:stp]
    x3 = Z[ist:ie:stp, jst:je:stp] 
    
    mx = dat.Mx[ist:ie:stp, phToSHow, jst:je:stp]
    mt = dat.Mt[ist:ie:stp, phToSHow, jst:je:stp]    
    mz = dat.Mz[ist:ie:stp, phToSHow, jst:je:stp]
    if dat.method == 'MHD':
        Bx = dat.Bx[ist:ie:stp, phToSHow, jst:je:stp]    
        Bz = dat.Bz[ist:ie:stp, phToSHow, jst:je:stp]
        print size(Bx), size(Bz)
        
    shape = (dat.nz,  dat.nt,  dat.nx)
        
    if (i==0):  
        var1 =  (dat.ro[ist:ie, phToSHow, jst:je]) 
        grid[0].streamplot(x1, x3, mx, mz, color='r', linewidth=2)
        
        grid[0].set_title(r'$\log(\rho)$', fontsize = 22)

    if (i==1):  
        var1 = log(dat.ro[ist:ie, phToSHow, jst:je]**(GAM-1.))
        var1 = log(dat.etot[ist:ie, phToSHow, jst:je]) #/ dat.ro[ist:ie, phToSHow, jst:je]**GAM )
        var1 =  log10(fabs(var1))
        grid[1].set_title(r'$\log(E_{\rm tot})$', fontsize = 22)
        
    if (i==2):  
                
        var1 = dat.Mt[ist:ie, phToSHow, jst:je]/Lk/(dat.ro[ist:ie, phToSHow, jst:je])                                
                
        var1 =  log10(fabs(var1))                
        
        grid[2].set_title(r'$\log( l/l_{\rm k})$', fontsize = 22)
        
    if (i==3):  
        varEm =  (dat.Bt[ist:ie, phToSHow, jst:je]**2+dat.Bx[ist:ie, phToSHow, jst:je]**2+\
                  dat.Bz[ist:ie, phToSHow, jst:je]**2)
        
        var1 =  log10(varEm) 
        grid[3].set_title(r'$\log( E_{\rm m})$', fontsize = 22)
        
    
    vatToShow2D=var1
    im = grid[i].imshow(vatToShow2D , interpolation='bilinear',cmap=cm.jet, 
                    extent=[xmin, xmx, zmin, zmx] )
    
    grid[i].set_xlabel ("R(pc)", fontsize=22)
    grid[i].set_ylabel ("z(pc)", fontsize=22)
     
       
        
    grid.cbar_axes[i].colorbar(im)


fileNameToSave = put_out+'fourPanelResultsFig2'
# fig.savefig(fileNameToSave + ".pdf", format='pdf')

show()
exit()
    
