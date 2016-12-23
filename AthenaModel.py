import sys
import numpy as nm
# from pylab import *
# from matplotlib import pyplot as plt
from physics1 import *
import pprint, os
import re
from numpy.distutils.misc_util import fortran_ext_match

class athDataModel:
  
  """ reads data from hdf file
        fields: r, rth, dd, ee, u1,u2,u3"""
        
  def __init__(self):
    self.method='MHD'   
    self.debug = False
    self.nx=None
    self.nt=None
    self.nz=None
    self.i_s=None
    self.ie=None
    self.js=None
    self.je=None
    self.dt_bin=None       
    self.x=None
    self.z=None
    self.ph=None
    self.ro=None
    self.etot=None    
    
    self.Mx=None
    self.Mt=None
    self.Mz=None
    
    self.Bx=None
    self.Bt=None
    self.Bz=None
    
    self.Mbh=0
    self.beta=0
    self.n0=0
    self.Rsc=0
    self.Usc=0.
    self.Dsc=0.
    self.tsc=None
    
    self.nvar =0    
    self.nscalars=0
    self.ifselfg=0
    self.ifpart=0
    self.par=None
  
  def loadSimulationParam(self, fileNameFullPath, print_res):
      
    locdir=os.getcwd()
    
#    print ("file=", fileNameFullPath)
    
#     
    try:

#         for dataLine in open( fileNameFullPath, 'r').read().split('\n'):
        for dataLine in open( fileNameFullPath, 'r').read().splitlines():
            
            print(dataLine)

            
            if 'Nx1' in dataLine:
                self.nx = int(re.findall('Nx1\s*=\s*(\d*)', dataLine)[0])            
                self.i_s=0
                self.ie=self.nx-1

    
            if 'Nx2' in dataLine:
                self.nt = int(re.findall('Nx2\s*=\s*(\d*)', dataLine)[0])
    #             print nx2
    
            if 'Nx3' in dataLine:
                self.nz = int(re.findall('Nx3\s*=\s*(\d*)', dataLine)[0])            
                self.js=0
                self.je=self.nz-1

            
            if 'nc0' in dataLine:            
                                 
                self.n0 = float(re.findall('nc0\s*=\s*(\d*\.\d*e\d*)', dataLine)[0])
                
#                 print re.sub(r'\.' , "", re.findall('nc0\s*=\s*(\d*\.\d*e\d*)', dataLine)[0])
#                 print re.findall('nc0\s*=\s*(\d*\.\d*e\d*)', dataLine)[0], "la"; 
                self.n0  =1.e8
                                
                self.Dsc = MP* self.n0            
                print ("self.Dsc=", self.Dsc)
                
            if 'F2Fedd' in dataLine:         
                F2Fedd = re.findall('F2Fedd\s*=\s*(\d*\.\d*)', dataLine)[0]                
                F2Fedd_str=re.sub(r'\.', "", F2Fedd)
                self.F2Fedd = float(F2Fedd)

            if 'dt' in dataLine:         
                dt_bin = re.findall('dt\s*=\s*(\d*\.\d*)', dataLine)[0]                
                dt_bin_str=re.sub(r'\.', "", dt_bin)
                self.dt_bin = float(dt_bin)               
#                print(" inforcing dt =0.1  this needs to be corrected \n" )
                self.dt_bin = 0.1

        
            if 'r0' in dataLine:                                 
                     self.Rsc =  float( re.findall('r0\s*=\s*(\d*.*\d)', dataLine)[0] )        
                    
                     
                     
                     print ("Rc=", self.Rsc)
                     rg = 2.*G*self.Mbh/CL**2
                     
                     self.Rsc *= rg
                     
                     self.Rsc = PC
                     
                     print ("self.Rsc is scaled ro 1pc")

            if 'M2Msun' in dataLine:                                 
                     self.Mbh =  float(  re.findall('M2Msun\s*=\s*(\d*.*\d)', dataLine)[0] )
#                      print ("M=", self.Mbh)        
                     self.Mbh *= MSUN  
                                       
            if 'beta' in dataLine:         
                beta = re.findall('beta\s*=\s*(\d*.\d*)', dataLine)[0]
                beta_str=re.sub(r'\.', "", beta)
                self.beta = float(beta)
#                 print("beta=", beta)
#             self.tsc=sqrt( self.Rsc**3/G/self.Mbh/MSOL)
#             
#            
#             
#             mass = self.Mbh*MSOL
#             
#             self.Esc=self.Dsc*G*mass/self.Rsc
#             self.Usc = sqrt( G*mass / self.Rsc )
   
        if print_res:
            print('Mbh =', self.Mbh);
            print( 'nc0 = %10.3e' % (self.n0))
            print("r0 = %10.3e cm "% (self.Rsc))      
            print('F2Fedd = %1.1f ' % self.F2Fedd)
            print( "beta = %10.2f" % (self.beta)  )    
            print( "dt_bin = %10.2f" % (self.dt_bin)  )    
            
            
                
    
    except Exception as e:
        print(str(e))
    
    self.Usc = nm.sqrt(G*self.Mbh/self.Rsc)            
    self.tsc=self.Rsc/self.Usc
#     print(self.tsc/YR)
    self.par= {'Dsc': self.Dsc, 'tsc':self.tsc,  'Usc':self.Usc, 'Mbh':self.Mbh,  'Usc' :self.Usc,
                    'nc0':self.n0,' Rsc': self.Rsc
                    }

# try:
#     for dataLine in open( locdir+'/torus1/zmp_inp', 'r').read().split('\n'):
#         if '&ggen1' in dataLine:
#             nbl1 = re.findall('nbl=(\d*),', dataLine)[0]            
#         if '&ggen2' in dataLine:
#             nbl2 = re.findall('nbl=(\d*),', dataLine)[0]            
# except Exception, e:
#     print(str(e))
     
    
      
       
  
  def loadDataFromBinFiles(self, fileToOpen, dat, printDetail = True):
    try:
        file = open(fileToOpen, 'rb')
    
    except IOError:
            print ('cannot open', fileToOpen)
            exit()

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
    if self.method == 'MHD':
        dat.Bx = nm.fromfile(file,dtype=nm.float,  count=count).reshape( shape )
        dat.Bt = nm.fromfile(file,dtype=nm.float,  count=count).reshape( shape )
        dat.Bz = nm.fromfile(file,dtype=nm.float,  count=count).reshape( shape )

    

  def massAccrRate(self, dat, mdot, phToSHow = 1):
         js =2                        
         
         mdot[:] = 2.*nm.pi* self.Rsc*self.x[js]* self.Mx[:, phToSHow,js]*self.Usc*self.Dsc
                                
         
         mdotTot =-nm.trapz(mdot, self.z*self.Rsc)         
          
         return(mdotTot)
  
  def massLossRate(self, dat, phToSHow = 1):

     #inner boundary
     mdot_a = nm.zeros(dat.nz)
     mdot_a[:] = 2.*nm.pi* dat.Rsc*dat.x[dat.js]* dat.Mx[:, phToSHow,dat.js]*dat.Usc*dat.Dsc                                         
     mdotTot_a =-nm.trapz(mdot_a, dat.z*dat.Rsc)         

     #  outer x boundary
     mout1 = nm.zeros(dat.nz)           
     mout1[:] =  2.*nm.pi* dat.Rsc*dat.x[dat.je]* dat.Mx[:, phToSHow,dat.js]*dat.Usc*self.Dsc                                       
     mout1[:]=[0. if x<0. else x for x in mout1 ]  
     
     mdotTot1 = nm.trapz(mout1, dat.z*dat.Rsc)
     mout2 = nm.zeros(dat.nx)           
     
     #  upper z boundary
     mout2[:] = 2.*nm.pi* dat.Rsc*dat.x[:]*self.Mz[dat.ie, phToSHow,:]*self.Usc*self.Dsc                                    
     mdotTotUp = nm.trapz(mout2, dat.x*dat.Rsc)
     
     mout2[:] = 0.
     #  lower z boundary
     mout2[:] = -2.*nm.pi* dat.Rsc*dat.x[:]*self.Mz[dat.i_s, phToSHow,:]*self.Usc*self.Dsc                           
     mdotTotBot = nm.trapz(mout2, dat.x*dat.Rsc)
          
     return(mdotTot1, mdotTotUp, mdotTotBot, mdotTot_a )
     
