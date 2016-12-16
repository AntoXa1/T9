#!/local/data/atorus1/dora/Compilers/epd-7.3-1-rh5-x86_64(1)/bin/python

import scipy
from  numpy import ndarray, zeros, array, size, meshgrid, flipud, floor, where, amin, argmin, int
import numpy as nm
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import os
import time
import pprint

from scipy import trapz
from matplotlib import pyplot as plt
from physics1 import *
import glob
import socket

import AthenaModel as ath
import time
import cPickle as pickle




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




def movingAverage(interval, window_size):
    window= nm.ones(int(window_size))/float(window_size)
    return ( nm.convolve(interval, window, 'same') )
    
    
    
        
def torusMass(dat):
    from math import pi
    
    m=0.
    m2 = 0
    Dsc2 = 1e-4*dat.Dsc #background
    
    for i in range(1, dat.Nz):        
        for j in range(1, dat.Nx):
        
            r =  dat.x[i]            
            dr = (dat.x[i] - dat.x[i-1])            
            dh = (dat.z[i] - dat.z[i-1])             
            dv = 2*pi*r*dr*dh                        
            dm2 = dat.Rsc**3* Dsc2 *dv* dat.dd[i,j]            
            dm = dat.Rsc**3* dat.Dsc *dv* dat.dd[i,j]
            
#             print(dat.Nz,dat.Nx)
#             print (r, dr, dh, dv, dm, dat.Dsc)
            
            m+=dm
            m2+=dm2
    print('torusMass=', m)
    return (m)
        
       
        
    
def getFilePostfixFromNum(num, needLen):
    s1 = str(num)
    return(s1.zfill(needLen))
       


def kinEnergyLoss(dat):    
    
    i0  = dat.Nz - 3
    dEkn =  pi* dat.x[:]*dat.Rsc * dat.dd[:,js]*dat.Dsc  * (dat.u1[i0, :]*dat.Usc)**3
    ekinTot =trapz(  2.*dEkn, dat.x*dat.Rsc) 
    return(ekinTot)


     
def velocAvr(dat, mtor):
        
    i_eq = nm.argmin(abs(dat.z-1.))
#     print(i_eq, '  ', mtor);exit()    
    vAvr=0.
        
    for i in range(i_eq, dat.Nz):        
        for j in range(1, dat.Nx):
            
            r =  dat.x[i]            
            dr = (dat.x[i] - dat.x[i-1])            
            dh = (dat.z[i] - dat.z[i-1])             
            dvol = 2*pi*r*dr*dh                        
            
            dvzAvr1 = dat.u1[i, j]*dat.dd[i,j] *dvol            
            dvzAvr1 *= dat.Rsc**3* dat.Dsc*dat.Usc
            
            #print(dvzAvr, dvol, dat.u1[i, j])
            
            vAvr += dvzAvr1
     
    res = vAvr/mtor    
    print('v1Avr=', res)    
    return(res)
                
    
#     v1av = 4*pi*dvav1 * dat.Usc**2 * dat.Dsc* (dat.Rsc**2) / dmdt1
            
whatToDo = 'processTauAndColDensFromFile'

dirFileToReadBase = ''
dataDir = ''



dat =ath.athDataModel()

if socket.gethostname()=='atorus':
    locdirList = ['/local/data/atorus1/dora/PROJECTS/AthenaWind/']

    locdir = locdirList[0]   
    dirToRead=locdir+'tst/cylindrical/'
    print(dirToRead)
    dat.loadSimulationParam(dirToRead + 'athinput.torus9_hydro_2D', print_res=True)
   
    dirToRead=locdir+'bin/'
    print(dirToRead)
    
    put_out= '/local/data/atorus1/dora/PROJECTS/SCRIPTS/T9'

else:
     locdirList = [ 'SolovievSep201615_256x8x256_L0.n10e10/']
     locdirList = [ 'runDec201608_256x8x256_L0.5n10e8/']
     put_out= '/Users/dora/WORK/ECLIPSE_SPACE/torus9'
     put_FIG = '/Users/dora/Documents/TEX/torus9/'
     locdir = locdirList[0]
     dirFileToReadBase = os.getcwd()
     dataDir = '/DATA/'

     dirToRead = dirFileToReadBase + dataDir+locdirList[0]
     
     dat.loadSimulationParam(dirToRead + 'athinput.torus9_hydro_2D', print_res=True)



#--------------------------------------
files = []
nFile=1.
maxNumFile=10000

       
mdotZi =[]
mdotAccrOverTime=[]
mdotOverTimeUD=[]
mdotOverTimeR=[]
ekinOverTime=[]
avrVelOverTime=[]
mv1=[]
simTime=[]
#--------------------------------------

mdotFromPickledFile = False 

plotMdotFromTxtFile = True


calcFromDataFiles = False
writeToTxtFile    =calcFromDataFiles


writeToPickledFile=False

# calcFromDataFiles = True
calcTorusMass = False
calcMdot = calcFromDataFiles

plotMdotNow = False
plotEkin = False
plotAvrVel =False

firstTime = 0
fileToReadPrefix="mhdXwind"

if plotMdotFromTxtFile:
    
    filesToRead = ['torus9_mdot_HW_tot.dat', 'torus9_mdot_SL_tot.dat']        
#     f = plt.figure()
    f = plt.figure(1, (15,7))  
    
    for filename,i in zip(filesToRead, range(len(filesToRead))):
                
        filename = put_out +'/'+ filename
    
        try:
            timeMdot=nm.loadtxt(filename)
            print('loading data from',filename)
        except IOError:
            print('cannot load data from',filename)
            exit()
        
        time = timeMdot[:,0]            
        mdotAccr = timeMdot[:,1]
        mdotWin = timeMdot[:,2]
        
        mdotAvr=mdotAccr
        ax = f.add_subplot(1,2,i+1)
#         if i==0:
        
#         if i==1:
#             ax = f.add_subplot(122)
        y1 = mdotAvr
        y2 = mdotWin        
        time*= dat.tsc/31536e3
        ax.plot(time,  log10(fabs(y1)))
        ax.plot(time,  log10(fabs(y2)), '--',linewidth=2)        
        ax.legend(('${\dot M_{a}}$', '${\dot M}_{w}$'), loc='lower right')
#         plt.xlim((0,113000))

        ax.set_xlabel ("$t [yr]$",  fontsize=16 )
        ax.set_ylabel ("log ${\dot M}[M_\odot/yr]$", fontsize=16)
        ax.ticklabel_format(style = 'sci', useOffset=True)
#     ax.ticklabel_format(style = 'sci', useOffset=False)
    f.suptitle('Accretion rate and wind mass-loss rate', fontsize=16)
    
    fileNameToSave = put_FIG+'mdotAccrAndWindTwoPanel'
#     f.savefig(fileNameToSave + ".pdf", format='pdf')
    
    show(); 
    exit()

    
if calcFromDataFiles :
    
    filelist =  glob.glob(os.path.join(dirToRead, 'mhdXwind*.bin') )
#     print(filelist)
  
#     for fileInDir in os.listdir(dirToRead):
    for fileInDir in sorted(filelist):
        print(fileInDir)        
#         if fileInDir.startswith(fileToReadPrefix) and  fileInDir.endswith("bin") :            
#         fileToOpen = dirToRead + fileInDir              
        fileToOpen =    fileInDir
        print("file to open:", fileToOpen)  
        dat.loadDataFromBinFiles(fileToOpen, dat, printDetail=False ) 
                    
                    
        if (calcTorusMass): 
            torMass= torusMass(dat)
            print("mass/Msol= ", torMass/MSOL)
            exit()

#         print(dat.Nz, dat.Nx); exit()        
        #ieq = int( round(dat.Nz/2., 0))
        ieq = nm.argmin(abs(dat.z-1.))
        
        if (calcMdot):                                 
#             torMass= torusMass(dat)
                
            [mdotTot1, mdotTotUp, mdotTotBot,mDotAccr] = dat.massLossRate(dat)
            print('mdotAccr=', mDotAccr/6.e25,
                     'out=', mdotTot1, 
                    'Up=', mdotTotUp, 
                    'Bot=', mdotTotBot)    
            
            mdotAccrOverTime.append(mDotAccr/6.e25)       
            mdotOverTimeUD.append((mdotTotUp +mdotTotBot)/6.e25)
            mdotOverTimeR.append((mdotTot1)/6.e25)            
            
        if (plotEkin): 
            ekinTot =kinEnergyLoss(dat)
            print('ekinTot=', ekinTot)
            ekinOverTime.append(ekinTot)
        
        if (plotAvrVel):
            if (firstTime ==0):
                torMass= torusMass(dat)
                firstTime =1
            vAvr=velocAvr(dat,torMass)
            avrVelOverTime.append(vAvr)
                              
            #vin1 = 0.2*( dat.u2[ieq, js] + dat.u2[ieq+1, 1] +dat.u2[ieq+2, 1] +dat.u2[ieq-1, 1]+dat.u2[ieq-2, 1] )  *dat.Usc        
            # mv1.append( vin1 )        
                        
        simTime.append(nFile)                                      
        
        nFile += 1
        if (nFile > maxNumFile):
            break
        mdotZi = []
    

    timeArr=nm.array(simTime)*dat.dt_bin
    totData = zip(timeArr,mdotAccrOverTime,mdotOverTimeUD, mdotOverTimeR)
    
    if writeToPickledFile:
        filename = 'torus9_mdot_tot'+'.p'
        print("mdot_accr pickle file",  filename)
        pickle.dump( totData, open(filename, "wb")) 
        print("mdot_accr saved to  ", filename)
        
    if writeToTxtFile:
        # filename =  put_out +'/'+ 'torus9_mdot_HW_tot'+'.dat'        
        filename =  put_out +'/'+ 'torus9_mdot_SL_tot'+'.dat'        
        try:
            nm.savetxt(filename,  list(totData))
            print("mdot_accr txt file",  filename)
        except IOError:
            print('cannot save to', filename)
            
            
    # filename = 'torus9_mdotaccr_'+'.p'
    # print("mdot_accr pickle file",  filename)
    # pickle.dump(zip(timeArr,mdotAccrOverTime), open(filename, "wb")) 
    # print("mdot_accr saved to  ", filename)
            
            
    # try:
    #     mdotAccrAvr=pickle.load( open (filename, "rb"))
    #     print("cannot open pickled file: mdot_accr")
    # 
    # except IOError:
    #     mdotAccrAvr = movingAverage(mdotOverTime, 20)
    #     pickle.dump(mdotAccrAvr, open(filename, "wb"))      
    #     print("mdot_accr saved to  ", filename)
        
        
    
    
    
    
    #ax.plot(timeArr,  avrVelOverTime)
    #ax.plot(timeArr,  log10(ekinOverTime))


    
if mdotFromPickledFile:   
    
#     f = plt.figure()
#     ax = f.add_subplot(111)

    filesToRead = ['torus9_mdot_tot.p']        
    iiter = range(size(filesToRead))
    lineType = ['-']
        
    t0 = 0    
    
    
    for filename, i in zip(filesToRead, iiter):
        try:            
            filename = dirFileToReadBase +'/'+ filename
            print(filename)                        
            timeMdot = nm.asarray( pickle.load( open( filename, "rb" ) ) )            
          
            time = timeMdot[:,0]            
            mdotAccr = timeMdot[:,1]
            mdotWin = timeMdot[:,2]
#             print mdotWin

            print("loaded mdot_accr form a pickled file",  filename)
   
        except IOError:
            print("cannot open pickled file", filename)
            
        convol_window = 10
        mdotAvr=mdotAccr
#         mdotAvr = movingAverage(mdotAccr, convol_window)
        ax = f.add_subplot(111)        
#         ax.plot(time,  log10(mdot))        
        y1 = mdotAvr
        y2 = mdotWin
        
#         y = mdot
#         t0 = max(t0, time[0])
        
        time*= dat.tsc/31536e3
        ax.plot(time,  log10(fabs(y1)))
        ax.plot(time,  log10(fabs(y2)), '--',linewidth=2)        
        ax.legend(('${\dot M_{a}}$', '${\dot M}_{w}$'), loc='lower right')
        plt.xlim((0,113000))
#         ax.legend(('A', 'B'), loc='lower right')
        
        
#         t0 = abs(time - 3).argmin()      
#         time = time[0 : t0]
# #         
#         ax.plot( time ,  log10(  abs(y[0:t0])  ) , lineType[i])
    
#     ax.annotate('$\Gamma=0.05$', xy=(6000,0),  xycoords='data',
#             xytext=(6000, 0), textcoords='data',            
#             horizontalalignment='right', verticalalignment='top',
#             )
#     
#     ax.annotate('$\Gamma=0.1$', xy=(6000,0),  xycoords='data',
#             xytext=(10000, -2.5), textcoords='data',            
#             horizontalalignment='right', verticalalignment='top',
#             )
     
     
    ax.set_xlabel ("$t [yr]$",  fontsize=16 )
    ax.set_ylabel ("log ${\dot M}[M_\odot/yr]$", fontsize=16)
    ax.ticklabel_format(style = 'sci', useOffset=True)
#     ax.ticklabel_format(style = 'sci', useOffset=False)
    ax.set_title("Accretion rate and outer disk mass-loss rate")
        


    fileNameToSave = put_out+'mdotAccrAndWindOnePanel'
#     f.savefig(fileNameToSave + ".pdf", format='pdf')
    show(); 
    exit()

if plotMdotNow:     
    f = plt.figure()    
    ax = f.add_subplot(111)

#     ax.plot(timeArr,  mdotOverTime)         
    y=log10(abs( nm.array(mdotAccrOverTime) ) )                     
#     y=mdotOverTime
    
#     ax.plot(timeArr,  y)
#     ax.plot(timeArr,  log10(mdotOverTimeR))
#     ax.plot(timeArr,  (mdotOverTimeUD))    
#     set_fonts_etc(ax)
    ax.plot(timeArr,  y)
    show()    
    exit()
    
    


# if (plotMassLoss):     
#     
#     mLost = trapz(mdotOverTime, timeArr* 4.8*10**3)
#     print("mLost=", mLost)
#     
#     ax.plot(timeArr,  (mdotOverTime))
#     
# #ax.plot(simTime, mv1)

plt.show()
