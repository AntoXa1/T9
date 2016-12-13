import scipy
from  numpy import ndarray, zeros, array, size, sqrt, meshgrid, flipud, floor, where, amin, argmin, int
import numpy as nm
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import pprint
import os
import glob
import time
import os
import time
import cPickle as pickle
from matplotlib import colors



from scipy import trapz
from matplotlib import pyplot as plt
from physics1 import *

import physics1 as ph

phToSHow = 1
import AthenaModel as ath


# def inclinationToPositionInGrid(dat, th):

def plotFormat(ax1,ax2=None,ax3=None,im=None):
#     from mpl_toolkits.axes_grid1 import make_axes_locatable

    fontsize_bar=14
    fz = 16
    
#     ax1.set_ylabel('$\tau N(\tau(\theta)>1)$', fontsize = 22)
#     ax1.set_title('Obscuring models', fontsize = 19)
#         
#     ax1.set_ylabel('log (Number of obscuring models)', fontsize = 19)
#     for ylabel in ax.get_yticklabels():
#         ylabel.set_fontsize(fontsize_x)
#     
    ax1.set_xlabel('Inclination angle', fontsize=fz)
    ax1.set_ylabel('$\log(N_{col})$', fontsize=fz)
#     for xlabel in ax1.get_xticklabels():
#         xlabel.set_fontsize(fontsize_x)
# 
    ax1.set_title('Column density $(cm^{-2})$', fontsize =fz)

    fig.subplots_adjust(wspace=0.2)       
#     plt.setp([a.get_yticklabels() for a in fig.axes[1:] ], visible=False)

#     ax1.set_ylabel('$log (N_{col})$', fontsize = 19)
    
#     for ylabel in ax.get_yticklabels():
#         ylabel.set_fontsize(fontsize_x)
#     
#     ax2.set_xlabel('Inclination angle', fontsize = 19)
#     for xlabel in ax.get_xticklabels():
#         xlabel.set_fontsize(fontsize_x)
    
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", size="5%", pad=0.05)                        
# 
#     if (im):
#         cb =plt.colorbar(im, cax=cax)        
#         for t in cb.ax.get_yticklabels():
#             t.set_fontsize(fontsize_bar)

    
def cellTrace(dat, th, Xsrc_x_z, phToSHow):
        tiny = 1e-18
        i_s = dat.i_s
        ie = dat.ie
        js = dat.js
        je = dat.je
        
        
                                
        zout = dat.x[je] / nm.tan(th) + Xsrc_x_z[1]
        
        ic = abs(dat.z - zout).argmin()          
        jc = je
        
#         xcur = nm.array([ dat.x[jc], dat.z[ic] ])
        
        xcur = nm.array([ dat.x[dat.je], zout ])                
        
        xsrc0= Xsrc_x_z                 
        
        a = xcur - xsrc0
        
        
        norm = (xcur - xsrc0)/nm.dot( xcur - xsrc0, xcur - xsrc0)
        #crossing x_{js} plane
        xz1 =nm.array([ dat.x[js+1], norm[1]/max(norm[0], tiny)*(dat.x[js+1] - xsrc0[0]) + xsrc0[1] ])            
                
#         i =  minloc((dat. z -xz1[1]) **2)        
        i = abs(dat.z - xz1[1]).argmin()        
        if i>=dat.ie or i==dat.i_s: return(0.,0.)
        
        j = dat.js
        xz0 = nm.array([ dat.x[j], dat.z[i] ])
        
        
        stp = 0
        dstot = 0.
        tau = 0.
        dtau = 0.
        dcol_dns = 0.
        col_dns = 0.

        dl = 0.
      
        while True:        
                
                xz1 =nm.array([ dat.x[j+1], #crossing x_{j+1} plane   
                    norm[1]/max(norm[0], tiny)*(dat.x[j+1] - xsrc0[0]) + xsrc0[1] ])            
                
                Sx = sqrt(nm.dot(xz1 - xz0, xz1-xz0))
                
                itm = i+  nm.int( nm.sign( norm[1] ) )# crossing z_{i+-1} plane 

                if (norm[1] != 0.):
                    xz2 =nm.array([ norm[0]/norm[1] *(dat.z[itm] - xsrc0[1]) + xsrc0[0], 
                                dat.z[itm] ])
                else:
                    xz2 =nm.array([ norm[0]/nm.max(norm[1],tiny)*(dat.z[itm] - xsrc0[1]) + xsrc0[0], 
                                dat.z[itm] ])
                
                Sz = sqrt(nm.dot(xz2 - xz0, xz2-xz0))

                if (Sz>Sx):
                    dl = Sx  #right
                    xz0 = xz1
                    j += 1
                else:
                    dl = Sz #up or down
                    xz0 = xz2
                    i = itm
                
                dstot += dl                                
                stp += 1
                opac = ph.KPE
                            
#                 print "here", dat.Rsc, dat.Dsc,  dat.ro[i,phToSHow, j],  opac, dl,i,j                
                                
                dtau = dat.Rsc*dat.Dsc* dat.ro[i,phToSHow,j]* opac * dl
            
                dcol_dns = nm.fabs( dat.ro[i, phToSHow, j]*dat.n0* dat.Rsc*dl)
                
                tau += dtau                
                col_dns += dcol_dns
#                 print i, j, stp
                
                if  ( i <= dat.i_s or i >= dat.ie or j <= dat.js or j >= dat.je-1):
                    break
                if  i==ic and j ==jc :
                    break
        
    
        return(tau, col_dns)
        
def transmisionFunctionsAngleGrid(dat, printResults=False):
    Nth = 100         
    
    Xsrc = array([0, 0])  
    
#     res = dat.x[dat.js]/dat.z[dat.ie]
#     thmax = nm.pi / 2
    thmin = nm.arctan2(dat.x[dat.js], dat.z[dat.ie]-Xsrc[0] )
    
    thmax =nm.arctan2(dat.x[dat.js], dat.z[dat.i_s]-Xsrc[0])
        
    angle = linspace(thmin, thmax, Nth) 
    
#     print angle, range(1,Nth)
#     time.sleep(3)
    
    tauTheta = zeros(Nth)
    colDens = zeros(Nth)
                                                 
    for k,th in zip(range(1,Nth), angle):        
#         th = nm.pi/2              
#         print k,th
#         time.sleep(3)
        tauTheta[k], colDens[k] = cellTrace(dat, th, Xsrc_x_z = [Xsrc[1], Xsrc[0]], phToSHow =phToSHow)        
#         print("tau=", tauTheta[k], "colDens=", colDens[k], th, Nth)
        
#         exit(); time.sleep(3)
        
        
    if printResults:
        fig  = plt.figure()
        ax = fig.add_subplot(111)
        funcToPlot =  log10(colDens)
        funcToPlot =  tauTheta
        ax.plot(angle*180./nm.pi, funcToPlot)
        show()
    
    return(tauTheta,colDens, angle)
    


def iteratorOverDataDirectoriesOverHDfFiles(dirFileToReadBase,
                                        locdirList2, funcToCalculate, fileToSavePrefix):
                                            
    mod={}
    maxNumFile = 100

    dat =ath.athDataModel()   
    
    for dirName, i in zip(locdirList2, range(size(locdirList2))):
        
        mod.update({ locdirList2[i]:{'ang':[], 'tau':[], 'cdens':[] }})
        
        simTime = []
        
        dirToRead = dirFileToReadBase + dirName
        
#         print "dirToRead=",  dirToRead
        
        
        
        nFile = 0.
        dat.loadSimulationParam(dirToRead + 'athinput.torus9_hydro_2D', print_res=True)   
        
        
        for fileInDir in os.listdir(dirToRead):
            
#             print("1  ",  fileInDir)
            if fileInDir.startswith("mhdXwind"):                         

                dat.loadDataFromBinFiles(dirToRead + fileInDir, dat, 
                                                printDetail=False)        
                
                print("file to open:", dirToRead + fileInDir)
                               
                tau, cdens, ang  = funcToCalculate(dat, printResults=False)
                
#                 print tau, cdens, ang
                                              
                
                mod[locdirList2[i]]['tau'].append(tau.tolist())
                mod[locdirList2[i]]['cdens'].append(cdens.tolist())
                                               
                if not mod[ locdirList2[i]]['ang']  :
                    mod[locdirList2[i]]['ang'].append(ang.tolist())
            
                simTime.append(nFile)                                                         
                nFile += 1
                if (nFile > maxNumFile):
                    print "maxNumFile reached"
                    break
                
        
        mod[ locdirList2[i] ].update( {'par': dat.par  })  
              
        mod[ locdirList2[i] ]['par'].update( {'dt_bin': dat.dt_bin })
       
    return(mod)
        
        
# ----------------------------------------------------------------------------                
#                                     MAIN            
# ----------------------------------------------------------------------------
            
# whatToDo = 'calculTauAndColDens'
whatToDo = 'processTauAndColDensFromFile'

dirFileToReadBase = os.getcwd()
dataDir = '/DATA'




locdirList2 = [ '/SolovievSep201615_256x8x256_L0.n10e10/']

# locdirList2 = ['/SolovievSep201614_128x8x128_L0.n10e8/']
# dirFileToReadBase ='/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind' 
# dataDir = ''
# locdirList2 = ['/bin/']

if whatToDo == 'calculTauAndColDens':

        multiDat=iteratorOverDataDirectoriesOverHDfFiles(dirFileToReadBase + dataDir, locdirList2,
                                            transmisionFunctionsAngleGrid, fileToSavePrefix=False)

#         for i_type_mod in xrange(len(locdirList2)):
            
        fileToSavePrefix ='multiDat_TauColDensVsAngle.p'        
        filename = dirFileToReadBase+dataDir +'/'+ fileToSavePrefix                    
        pickle.dump(multiDat, open(filename, "wb"))                     
        print("saved to  ", filename)        
        
        multiDat=[]
        
        print("calculTauAndColDens ..done"); exit()

if whatToDo =='processTauAndColDensFromFile':
    filename = dirFileToReadBase+dataDir+'/'+ 'multiDat_TauColDensVsAngle.p'
    
    mdat = pickle.load( open( filename, "rb" ) )                             
    Nd = len(locdirList2)    
        
    Na = len(mdat[locdirList2[0]] ['ang'][0])                
    Ny = len(mdat[locdirList2[0]] ['tau'])
    Nm=0
    distribFun1 = zeros(Na) 
    distribFun2 = zeros(Na)
    col_max=0
    col_min = 1e30
    
    fig  = plt.figure()
    lineType = ['-', '--', '-o', '-*']
    color = ['k', 'b', 'g', 'r']
    markerWidth = [2,2,2,2]
    
    for i_types in xrange(Nd):
                        
        for j_y in  xrange(Ny):            
            for k_ang in  xrange(Na):
                angl = mdat[locdirList2[i_types]] ['ang'][0][k_ang]
                tau= mdat[locdirList2[i_types]] ['tau'][j_y][k_ang]                
                col= mdat[locdirList2[i_types]] ['cdens'][j_y][k_ang]                                             

                if col > col_max:
                    col_max = col                
                if col <= col_min:
                    col_min = col                                   
                Nm +=1
                if tau>1.:
                    distribFun1[k_ang] +=1

        
        NBINS = 10
        colDens, dColdens = nm.linspace(col_min, col_max, NBINS, retstep=True)
        
        NANG = 50
        angle = nm.linspace(0., nm.pi, NANG, retstep=True)
        
        ang_scat = []
        Ncol_scat = []
                
        simTime = []
        
        told  = 0
        dt_bin = mdat[locdirList2[0]] ['par']['dt_bin']
         
        for j_y in  xrange(Ny):            
            
            tnew = told + dt_bin 
            simTime.append(tnew)
            told=tnew
            
            for k_ang in  xrange(Na):
                
                angl = mdat[locdirList2[i_types]] ['ang'][0][k_ang]                                                    
                col= mdat[locdirList2[i_types]] ['cdens'][j_y][k_ang]
                
                Ncol_scat.append(col)
                ang_scat.append(angl)
                                      
        x = nm.array(mdat[locdirList2[0]]['ang'][0])    
        imax=  distribFun1.argmax()    
        thmax = x[imax] - nm.pi/2.    
        x=x-thmax        
        
#         print(thmax)
        
        
        loc1='lower center'                        
        simTimeArr =nm.array(simTime)
                     
        print len(Ncol_scat)
#         eps = 1.e-2
#         Ncol_scat[Ncol_scat<eps]=eps
        
        
        if i_types==0:
            ax1 = fig.add_subplot(111)
            
#             color = [str(y/100. ) for y in simTimeArr]            
#             color = matplotlib.cm.rainbow(color)                        
#             print color        
#             color = [simTime,simTime]
#             x = np.random.rand(100)
#             y = np.random.rand(100)
#             t=x            
#             print t
#             plt.scatter(x, y, c=t)
#             plt.show()               
            
            print nm.array(ang_scat)

            x1 =(nm.array(ang_scat))*180./nm.pi
            
            y1= log10( Ncol_scat )
            
#             print  len(x),len(y), len(simTime);exit()
            
            cm = plt.cm.get_cmap('RdYlBu')
            ax1.scatter( x1, y1, cmap=cm)
            plt.xlim((0,185))
            
#             fig.colorbar(sc)
#             ax1.scatter( x,y, c=simTime)
 
                    
# #             
#         ax1 = fig.add_subplot(131)      
#         ax1.plot(x*180./nm.pi, log10(distribFun1),  lineType[i_types], color=color[i_types], linewidth=markerWidth[i_types]        

#         ax1.legend(('$\Gamma=0.01$', '$\Gamma=0.05$', '$\Gamma=0.1$', '$\Gamma=0.3$'), 
#                         shadow = False, loc = loc1)
        if i_types==1:
            ax2 = fig.add_subplot(142,sharey=ax1)
            ax2.scatter((nm.array(ang_scat)-thmax)*180./nm.pi, log10(Ncol_scat), color=color[i_types])
        
#         ax2.legend(('$\Gamma=0.01$', '$\Gamma=0.05$', '$\Gamma=0.1$', '$\Gamma=0.3$'), 
#                         shadow = False, loc = loc1)
        if i_types==2:
            ax3 = fig.add_subplot(143,sharey=ax1)
            ax3.scatter((nm.array(ang_scat)-thmax)*180./nm.pi, log10(Ncol_scat), color=color[i_types])
#         ax2.legend(('$\Gamma=0.01$', '$\Gamma=0.05$', '$\Gamma=0.1$', '$\Gamma=0.3$'), 
#                         shadow = False, loc = loc1)
        if i_types==3:
            ax3 = fig.add_subplot(144,sharey=ax1)
            ax3.scatter((nm.array(ang_scat)-thmax)*180./nm.pi, log10(Ncol_scat), color=color[i_types])

    
    plotFormat(ax1)
    
    put_out= '/Users/dora/Documents/TEX/torus9/'    
    fileNameToSave = put_out+'colDensInclAngleAllModels'
    fig.savefig(fileNameToSave + ".pdf", format='pdf')

    show()
    

print('END OF PROGRAM')

#                 idxa = (np.abs(angle-angl)).argmin()
#                 idxc = (np.abs(colDens-col)).argmin()
       
#     Ncol = nm.sort(distribFun2)
#     indx = nm.argsort(distribFun2)    
#     for i in indx:                
#     cmin = distribFun2[0]    
#     cmax = distribFun2[len(distribFun2)-1]    

#     for i_type_mod in xrange(len(locdirList2)):
#         print("model", i_type_mod)            
#         for yin_dir, j in zip(mdat[locdirList2[i_type_mod]]['tau'], xrange(N)):                            
#             ax.plot(x[0], log10(yin_dir))
                                   