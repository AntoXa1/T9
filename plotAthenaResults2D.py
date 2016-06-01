import sys
import numpy as nm
from pylab import *
from matplotlib import pyplot as plt
import pprint

class athDataModel:
  
  """ reads data from hdf file
        fields: r, rth, dd, ee, u1,u2,u3"""
        
  def __init__(self):   
    self.nx=None
    self.nt=None
    self.nz=None
           
    self.x=None
    self.z=None
    self.ph=None
    
    self.dd=None
    self.etot=None    
    
    self.Mx=None
    self.Mt=None
    self.Mz=None
    
    self.Bx=None
    self.Bt=None
    self.Bz=None
    self.nvar =0    
    self.nscalars=0
    self.ifselfg=0
    self.ifpart=0
    
method = 'GAS'   
method = 'MHD' 
STREAM = True
# STREAM = False

pathToBin = '/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind/bin/RES1/mhdwind1/'
pathToBin = '/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind/bin/'

# fileToOpen=pathToBin + 'HKDisk.0022.bin'
fileToOpen=pathToBin + 'mhdXwind.0092.bin'  

try:
    file = open(fileToOpen, 'rb')
except IOError:
    print 'cannot open', fileToOpen
    exit()

# file.seek(0,2)
# eof = file.tell()
# file.seek(0,0)

dat = athDataModel()

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
    print size(Bx), size(Bz)

# print(size(Bx,0),size(Bx,1)); exit()
shape = (dat.nz,  dat.nt,  dat.nx)
var1 = zeros(dat.nx, dat.nz)

var1 = dat.Mt[ist:ie, phToSHow, jst:je]/dat.ro[ist:ie, phToSHow, jst:je]

var1 =  (dat.ro[ist:ie, phToSHow, jst:je])
var1 =  log10(var1)


# var1 =  dat.Bt[ist:ie, phToSHow, jst:je]

# dat.Mt[ist:ie, phToSHow, jst:je]


vatToShow2D=var1
# vatToShow2D=dat.etot[ist:ie, phToSHow, jst:je]
# vatToShow2D=mt
  
im = plt.imshow(vatToShow2D , interpolation='bilinear',cmap=cm.jet, 
                extent=[xmin, xmx, zmin, zmx] )

plt.colorbar()

# im = plt.imshow( dat.Bz[ist:ie, phToSHow, jst:je]  , interpolation='bilinear',cmap=cm.jet, 
#                 extent=[xmin, xmx, zmin, zmx] )    
if method == 'MHD' and STREAM:
    plt.streamplot(x1, x3, Bx, Bz, color='k',linewidth=2)

qp1 = plt.quiver(x1, x3, mx/dat.ro[ist:ie, phToSHow, jst:je], mz/dat.ro[ist:ie, phToSHow, jst:je], width=0.008, scale=2,                            
        pivot='mid', color='black',
        units='x' , headwidth =5, headlength =7,
        linewidths=(0.5,), edgecolors=('black'))

# qp1 = plt.quiver(x1, x3, Bx, Bz, width=0.008, scale=0.8,                            
#         pivot='mid', color='black',
#         units='x' , headwidth =5, headlength =7,
#         linewidths=(0.5,), edgecolors=('black'))


show()
exit()

                     
                     
show()

# print ro.shape
exit()







if file.tell() != eof: print 'Error: Too few bytes read.'



file.close()
print("done")