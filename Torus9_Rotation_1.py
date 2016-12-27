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

put0= '/Users/dora/WORK/ECLIPSE_SPACE/torus9/DATA/'
put_out= '/Users/dora/Documents/TEX/torus9/'

locdir = 'SolovievSep201615_256x8x256_L0.n10e10/'
put0 +=locdir



dat = AthenaModel.athDataModel()
dat.loadSimulationParam(put0 + '/athinput.torus9_hydro_2D', print_res=True)

 
timeStr_list=['0063', '0120']


MBH = dat.Mbh
# n0=dat.n0
F2Fedd=dat.F2Fedd
Ledd = 4*pi*CL*G*dat.Mbh/KPE


tsc=sqrt(dat.Rsc**3/G/dat.Mbh)
Lx = Ledd* dat.F2Fedd
print("MBH=", dat.Mbh,   'R0=', dat.Rsc,  'n0=', dat.n0, 'F2Fedd=', dat.F2Fedd,
               "t_0=", tsc/YR, 'YRs' ,
               "L= ", Lx)


fileNamePrefix ="mhdXwind."
strmPlot = False




arrayTimeTeX = []
lenFileList = len(timeStr_list)
timeNumeric = zeros(lenFileList) 
fileToOpenList = []                   

for timeStr, i in zip(timeNumeric,xrange(len(timeNumeric))):    
    timeNumeric[i] = float(timeStr_list[i])
    timePhys =   timeNumeric[i]*tsc*dat.dt_bin/YR    
    arrayTimeTeX.append(roundThenStringToLatexFormat(timePhys))              
    print('time in yrs=',  timePhys, arrayTimeTeX)
    timeStr= formatFileName(timeNumeric[i])
    
    fileToOpen = put0 +fileNamePrefix+timeStr_list[i]+'.bin'
    fileToOpenList.append(fileToOpen)
              
    print timeStr

#**************************


# for timeStr in timeStr_list:
if what2Do == 'calcTorusMass':
    mTorus= dat.torusMass()
    print('torusMass=', mTorus/MSUN); exit()
    
    Nz = dat.dd[:,10].size
    plt.plot(dat.x, dat.Usc/1e5*dat.u1[Nz-2,:])
    show() 
    exit()

fig = plt.figure(1, (10, 4))

grid = AxesGrid(fig, 111, # similar to subplot(132)
                    nrows_ncols = (1, 2),
                    axes_pad=(0.3, 0.),
                    share_all=True,
                    label_mode = "L",
                    add_all=True #,       
#                     cbar_pad = 0.,             
#                     cbar_location = "right",
#                     cbar_mode="each",
                    )
print fileToOpenList
for fileToOpen, i in zip(fileToOpenList, range(lenFileList)):
    dat.loadDataFromBinFiles(fileToOpen, dat, printDetail = True)
    print(fileToOpen)
    var1 = zeros(dat.nx, dat.nz)    
    (X, Z) = meshgrid(dat.x,dat.z)
    
    
    offset= 4
    offset_i= offset+50
    ist = offset_i;  ie = dat.nz-offset_i;
    jst = offset;      je =dat.nx-offset
    stp=1
    xmin = dat.x[jst]; xmx = dat.x[je-1]
    zmin = dat.z[ist]; zmx = dat.z[ie-1]

# (X1, Z1) = meshgrid(dat.x[jst:je], dat.z[ist:ie])
# Lk = zeros(dat.nx, dat.nz)
# Lk = Z1**2 + X1**2

    phToSHow = 1
    
    zToShow = [-1. , -0.5,  0.,  0.5,  1.]
#     f = plt.figure()
#     ax = f.add_subplot(111)

    for z_i in zToShow:
        
        i_z = nm.argmin( fabs ( dat.z - z_i))    
        print("i_z=",  i_z, dat.z[i_z])
        
    # print ("arg= ",  i_z , dat.z[i_z]); exit()
    # exit()
    
        i_show = i_z
    
        var = dat.Mt[i_show, phToSHow, jst:je]/dat.ro[i_show, phToSHow, jst:je]  #/Lk
        grid[i].plot(dat.x[jst:je],   var )
    
    
    var = 1./sqrt(dat.x[jst:je]) 
    grid[i].plot(dat.x[jst:je], var,  '--', color ='k')
    grid[i].legend(('z=-1', 'z=-0.5', 'z=0',   'z=0.5',  'z=1', '$V_{\\rm K} $'),fontsize=10)
    grid[i].set_title('$V_\phi$ (t='+ arrayTimeTeX[i]+' yr)'  , fontsize=18)
    grid[i].set_xlabel("R(pc)", fontsize=18)
    grid[i].set_ylabel("z(pc)", fontsize=18)
    
    
    fileNameToSave = put_out+'rotationVelocityFig_1'
    plt.savefig(fileNameToSave + ".pdf", format='pdf')
plt.show()
exit()


