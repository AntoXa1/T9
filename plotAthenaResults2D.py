#!/local/data/atorus1/dora/Compilers/epd-7.3-1-rh5-x86_64(1)/bin/python


##!/Users/dora/anaconda3/bin/python3

##!#!/Users/dora/Library/Enthought/Canopy_32bit/User/bin/python

# import sys

import numpy as nm
from pylab import *
from matplotlib import pyplot as plt
import pprint

import AthenaModel as ath
    
method = 'GAS'   
method = 'MHD' 
STREAM = True
# STREAM = False

pathToBin = '/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind/bin/RES1/mhdwind1/'
pathToBin = "/Users/dora/WORK/ECLIPSE_SPACE/torus9/DATA/SolovievSep201615_256x8x256_L0.n10e10/"


pathToBin = '/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind/bin/'
pathToBin = "/Users/dora/WORK/ECLIPSE_SPACE/torus9/DATA/runSolDec201612_256x8x256_L0.5n10e10/"



pathToBin ='/local/data/atorus1/dora/PROJECTS/AthenaWind/bin/'


# pathToBin = '/Users/dora/WORK/ECLIPSE_SPACE/torus9/DATA/runDec201608_256x8x256_L0.5n10e8/'


# pathToBin = '/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind/bin/lowresDat/'
# fileToOpen=pathToBin + 'HKDisk.00228.bin'

fileToOpen=pathToBin + 'mhdXwind.0228.bin'  



try:
    file = open(fileToOpen, 'rb')
except IOError:
    print('cannot open', fileToOpen)
    exit()


dat = ath.athDataModel()

#  READ COORDINATE SYSTEM INFORMATION
coordsys = nm.fromfile(file, dtype=nm.int32,  count=1)[0]

dat.nx, dat.nt, dat.nz, dat.nvar,dat.nscalars, dat.ifselfg, dat.ifpart =  nm.fromfile(file,dtype=nm.int32,count=7)
# ifselfg, ifpart = particles and selfgravity
gamma1,  iso_csound = nm.fromfile(file,dtype=nm.float_,count=2)


t,dt = nm.fromfile(file,dtype=nm.float_,count=2)
dat.x = nm.fromfile(file,dtype=nm.float,count=dat.nx)
dat.ph = nm.fromfile(file,dtype=nm.float,count=dat.nt)
dat.z = nm.fromfile(file,dtype=nm.float,count=dat.nz)

shape = (dat.nz,  dat.nt,  dat.nx)
print(shape)


count = nm.prod(shape)

dat.ro =  nm.fromfile(file,dtype=nm.float,  count=count).reshape( shape )
dat.Mx = nm.fromfile(file,dtype=nm.float,  count=count).reshape( shape )
dat.Mt = nm.fromfile(file,dtype=nm.float,  count=count).reshape( shape )
dat.Mz = nm.fromfile(file,dtype=nm.float,  count=count).reshape( shape )
dat.etot = nm.fromfile(file,dtype=nm.float,  count=count).reshape( shape )
if method == 'MHD':
    dat.Bx = nm.fromfile(file,dtype=nm.float,  count=count).reshape( shape )
    dat.Bt = nm.fromfile(file,dtype=nm.float,  count=count).reshape( shape )
    dat.Bz = nm.fromfile(file,dtype=nm.float,  count=count).reshape( shape )

(X, Z) = meshgrid(dat.x,dat.z)

offset= 2
ist = offset;  ie = dat.nz-offset;
jst = offset;      je =dat.nx-offset

# print(je)

stp=1

xmin = dat.x[jst]; xmx = dat.x[je-1]
zmin = dat.z[ist]; zmx = dat.z[ie-1]


phToSHow = 1

# plt. matshow(datToPlot  )

f = plt.figure()    
ax = f.add_subplot(111)

x1 = X[ist:ie:stp, jst:je:stp]
x3 = Z[ist:ie:stp, jst:je:stp]

mx = dat.Mx[ist:ie:stp, phToSHow, jst:je:stp]
mt = dat.Mt[ist:ie:stp, phToSHow, jst:je:stp]    
mz = dat.Mz[ist:ie:stp, phToSHow, jst:je:stp]
if method == 'MHD':
    Bx = dat.Bx[ist:ie:stp, phToSHow, jst:je:stp]    
    Bz = dat.Bz[ist:ie:stp, phToSHow, jst:je:stp]
    print (size(Bx), size(Bz))

# print(size(Bx,0),size(Bx,1)); exit()
shape = (dat.nz,  dat.nt,  dat.nx)
var1 = zeros(dat.nx, dat.nz)

Vk= sqrt(1./ (dat.x[ist:ie]**2 +dat.z[jst:je] **2)    )
var1 = dat.Mt[ist:ie, phToSHow, jst:je]/dat.ro[ist:ie, phToSHow, jst:je]/Vk


# var1 =  (dat.Bt[ist:ie, phToSHow, jst:je])
var1 =  (dat.ro[ist:ie, phToSHow, jst:je])
vatToShow2D =  log10(var1)
# var1 =  (dat.Bt[ist:ie, phToSHow, jst:je]**2+dat.Bx[ist:ie, phToSHow, jst:je]**2+dat.Bz[ist:ie, phToSHow, jst:je]**2)

# dat.Mt[ist:ie, phToSHow, jst:je]


# vatToShow2D=log10(var1/dat.etot[ist:ie, phToSHow, jst:je])
# vatToShow2D=dat.etot[ist:ie, phToSHow, jst:je]
# vatToShow2D=mt
  
im = plt.imshow(vatToShow2D , interpolation='bilinear',cmap=cm.jet, 
                extent=[xmin, xmx, zmin, zmx] )

plt.colorbar()

# im = plt.imshow( dat.Bz[ist:ie, phToSHow, jst:je]  , interpolation='bilinear',cmap=cm.jet, 
#                 extent=[xmin, xmx, zmin, zmx] )    


if method == 'MHD' and STREAM:
#     plt.streamplot(x1, x3, Bx, Bz, color='k',linewidth=2)
#    plt.streamplot(x1, x3, mx, mz, color='k',linewidth=2)

         qp1 = plt.quiver(x1, x3, mx/dat.ro[ist:ie, phToSHow, jst:je], mz/dat.ro[ist:ie, phToSHow, jst:je], width=0.008, scale=10,                            
         pivot='mid', color='black',
         units='x' , headwidth =5, headlength =7,
         linewidths=(0.5,), edgecolors=('black'))

# qp1 = plt.quiver(x1, x3, Bx, Bz, width=0.008, scale=0.8,                            
#         pivot='mid', color='black',
#         units='x' , headwidth =5, headlength =7,
#         linewidths=(0.5,), edgecolors=('black'))


show()
#exit()

                     
                     
#show()

# print ro.shape
#exit()


#if file.tell() != eof: print ('Error: Too few bytes read.')

#file.close()
#print("done")
