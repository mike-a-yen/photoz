# Run aperture photometry on images on all filters
# Detect object using SEP
# Generate RA, DEC, Flux, Flux Error for each object
# Match objects across different filters assign an ID
# Make catalog ready to be passed to BPZ and EAZY
import os
import sys
import time
import re
from subprocess import call
import logging
import shutil
from seechange.controller import Controller
from seechange.config import config
from seechange.model import *
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import quad
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import sep
import cPickle
import progressbar as pb
from colors import bcolors
import threedhst.eazyPy as eazy

widgets = [pb.ETA(),pb.Percentage()]

class ClusterData(object):
    def __init__(self,clusterName,detectionFilter='F105W'):
        ''' Create a catalog from data inside see_change_db containing Mag 
        and MagErr for cluster which is then used to run BPZ and EAZY to
        generate files containing P(z) for each object. BPZ and EAZY require 
        multiple 'helper' files in order to set parameters, etc. Then run these
        files through matplotlib to make P(z) plots for each object. This will 
        create all necessary directories, required to run.
        '''
        self.detectionFilter = detectionFilter
        #load latest images from database in all filters
        self.clusterName = clusterName
        self.IDnum = self.ClusterID(verbose=True)
        self.CheckRequiredDirectories()
        config.load()
        controller = Controller()
        self.field = controller.get_field(name=self.clusterName)
        deepsets = [self.field.get_deepset(filter_name=str(filter['filter'])) for filter in self.field.filters]
        self.filters = [str(filter['filter']) for filter in self.field.filters]
        self.CheckDetectionFilter()
        self.images = {str(image.filter):image.open() for image in deepsets}
        hduList = {str(image.filter):image.hdulist for image in self.images.values()}
        self.dataList = {filter:hdu[1].data for filter,hdu in hduList.items()}
        self.dataList = {filter:data.byteswap(True).newbyteorder() for filter,data in self.dataList.items()}
        for hdu in hduList.values(): hdu.close()        
        self.zeroPoints = {filter:self.ZeroPoint(filter) for filter in self.filters}
        #set zero point for and pixel size
        # perform background subtraction on all filters
        self.bkgRMS = {filter:self.Background(filter) for filter in self.filters}
        #perform sep obj detection on detection filter image
        self.objs = self.ImageObjDetection(self.detectionFilter)
        self.detectPixelSize = self.images[self.detectionFilter].pixel_scale
        
        genTypes = ['ID','x','y','RA','DEC']
        magTypes = list(np.array([(filter+'mag',filter+'magerr') for filter in self.filters]).flatten())
        dataTypes = {'names':genTypes+magTypes,'formats':['>i4']+['float' for i in xrange(len(genTypes+magTypes)-1)]}
        nObjs = len(self.objs)
        self.objInfo = np.zeros((nObjs),dtype=dataTypes)
        for filter in self.filters:
            if filter ==  self.detectionFilter:
                self.objInfo['ID'] = self.IDnum + np.arange(0,int(nObjs))
                self.objInfo['x'] = self.objs['x']
                self.objInfo['y'] = self.objs['y']
                self.objInfo['RA'] = self.objs['RA']
                self.objInfo['DEC'] = self.objs['DEC']
                self.objInfo[filter+'mag'],self.objInfo[filter+'magerr'] = self.AperaturePhoto(filter,self.objs)
            else:
                scale = self.detectPixelSize/self.images[filter].pixel_scale
                objects = self.objs[:]
                objects['x'], objects['y'] = self.images[filter].rd_to_xy(self.objs['RA'],
                                                                          self.objs['DEC'])
                objects['a'] = scale*self.objs['a']
                objects['b'] = scale*self.objs['b']
                self.objInfo[filter+'mag'],self.objInfo[filter+'magerr'] = self.AperaturePhoto(filter,objects)

    def FromScratch(self):
        ''' Run all checks and create all necessary directories
            and files to make P(z) plots.
            ***Best to run this if you are not sure what to do***
            More detailed instructions in HOWITWORKS
        '''
        self.CheckRequiredDirectories()
        self.CheckDetectionFilter()
        self.WriteCatalog()
        self.WriteColumns()
        self.ParamBPZ()
        self.CheckBPZRequirements()
        self.ParamEAZY(ASCII=False)
        self.TranslateEAZY()
        self.ZeroPointEAZY()
        self.CheckEAZYRequirements()
        self.RunBPZ()
        self.RunEAZY()
        # ASCII = False for Binary (speed)
        # ASCII = True for ASCII (human readable)
        self.ProcessBPZOutput(ASCII=False)
        self.ProcessEAZYOutput(ASCII=False)
        self.CombinePZ(ASCII=False)
        self.PlotPicklePZ('BPZ',ASCII=False)
        self.PlotPicklePZ('EAZY',ASCII=False)
        self.PlotPicklePZ('COMBINED',ASCII=False)
        self.PlotAllMethods()

    def RunBPZ(self):
        goodToGo = self.CheckBPZRequirements()
        if goodToGo == False:
            print bcolors.WARNING+'Can not Run BPZ yet, missing required files'+bcolors.ENDC
        else:
            print 'Calling BPZ......'
            #callStr = 'python /home/scpdata05/clustersn/photoz/bpz-1.99.3/bpz.py CATALOGS/%s.cat -P CATALOGS/%s.pars'%(self.clusterName,self.clusterName)
            callStr = 'python BPZSRC CATALOGS/%s.cat -P CATALOGS/%s.pars'%(self.clusterName,self.clusterName)
            
            print callStr
            call(callStr,shell=True)

    def RunEAZY(self):
        goodToGo = self.CheckEAZYRequirements()
        if goodToGo == False:
            print bcolors.WARNING+'Can not run EAZY yet, missing required files'+bcolors.ENDC
        else:
            print 'Calling EAZY......'
            # more details in ParamEAZY function or .param file
            #callStr = '/home/scpdata05/clustersn/photoz/eazy-photoz/src/eazy -p CATALOGS/%s.param -t CATALOGS/%s.translate -z CATALOGS/%s.zeropoint'%(self.clusterName,self.clusterName,self.clusterName)
            callStr = './EAZYSRC -p CATALOGS/%s.param -t CATALOGS/%s.translate -z CATALOGS/%s.zeropoint'%(self.clusterName,self.clusterName,self.clusterName)
            print callStr
            call(callStr,shell=True)
            
    def ProcessBPZOutput(self,ASCII=False):
        ''' Rewrite the standard BPZ output to a more useable format.
            Writes a .pkl file with the obj P(z)
            pz = {ID:np.array(P(z)),z:np.array(z)}
            Saved as binary(fast) or ASCII(human readable)
        '''
        if ASCII == False: ending = '.pklb'
        elif ASCII == True: ending = '.pkl'

        fileName = 'BPZ/'+self.clusterName+'.probslite'
        if os.path.exists(fileName) == False:
            print bcolors.WARNING+fileName+' does not exist\n Need to run RunBPZ()'+bcolors.ENDC
        else:
            f = open(fileName,'r')
            lines = f.readlines()
            # lines are formated as
            # ID P(z) DataPoints
            # xrange for P(z) is x=range(0.01,10.01,0.01)
            pz = {}
            z = np.arange(0.01,10.01,0.01)
            pz['z'] = z[:]
            for line in lines:
                if '#' in line: continue
                id = line.split()[0]
                data = np.array([float(num) for num in line.split()[1:]])
                pz[id] = self.NormalizePZ(data,z)

            WritePickle(pz,'BPZ/pzpickles/'+self.clusterName+ending,ASCII=ASCII)
            print 'Saved BPZ Pickle as '+'BPZ/pzpickles/'+self.clusterName+ending

    def ProcessEAZYOutput(self,ASCII=False):
        ''' Rewrite the standard EAZY output to a more useable format.
            Must have the EAZY binary output. Set in EAZY .param file.
            Writes a .pkl files with the obj P(z).
            pz = {ID:np.array(P(z)),z:np.array(z)}
            Saved in binary or ASCII
        '''
        if ASCII == False: ending='.pklb'
        elif ASCII == True: ending='.pkl'
        if os.path.exists('EAZY/%s/%s.pz'%(self.clusterName,self.clusterName)) == False:
            print bcolors.FAIL+'Can not Process EAZY'+bcolors.ENDC
            print bcolors.FAIL+'Do not have standard EAZY binary output'+bcolors.ENDC
            raise RuntimeError('Must create binary ouput from EAZY, set in .param file, and run cl.RunEAZY()')
        if os.path.exists('CATALOGS/%s.cat'%self.clusterName) == False:
            raise RuntimeError('Must write catalog, try running cl.WriteCatalog()')
        fi = np.genfromtxt('CATALOGS/%s.cat'%self.clusterName,names=True)
        pz = {}
        for idx,obj in enumerate(fi):
            zrange, data = eazy.getEazyPz(idx,
                                          MAIN_OUTPUT_FILE='EAZY/%s/%s'%(self.clusterName,self.clusterName),
                                          OUTPUT_DIRECTORY='./')
            if pz.has_key('z') == False:
                pz['z'] = zrange
            elif np.all(pz['z']==zrange) == False:
                print bcolors.FAIL+'EAZY objects have different z ranges!'+bcolors.ENDC
            id = str(int(obj['ID']))
            pz[id] = data
        WritePickle(pz,'EAZY/pzpickles/'+self.clusterName+ending,ASCII=ASCII)
        print 'Saved EAZY Pickle as '+'EAZY/pzpickles/'+self.clusterName+ending


    def CombinePZ(self,ASCII=False):
        ''' Combine the P(z)'s from BPZ and EAZY.
            Done by interpolating with 1-D spline.
            Output a pickle file containing the 
            renormalized P(z) from BPZ and EAZY.
            pz = {'z':zrange,ID#:P(z)}
        '''
        if ASCII == True: ending = '.pkl'
        if ASCII == False: ending = '.pklb'
        bpz = self.GetPickledOutput('BPZ')
        eazy = self.GetPickledOutput('EAZY')
        goodIDs,badIDs = self.CheckObjects(bpz.keys(),eazy.keys())
        pz = {'z':np.arange(0.01,10.01,0.01)}
        print bcolors.CYAN+'Applying Spline to object P(z)......'+bcolors.ENDC
        for id in goodIDs:
            if id == 'z': continue
            splBPZ = InterpolatedUnivariateSpline(bpz['z'],bpz[id],k=1)
            splEAZY = InterpolatedUnivariateSpline(eazy['z'],eazy[id],k=1)
            splCom = np.multiply(splBPZ(pz['z']), splEAZY(pz['z']))
            normalization = 1.0/np.trapz(splCom,pz['z'])
            pz[id] = normalization*splCom
        WritePickle(pz,'COMBINED/pzpickles/'+self.clusterName+ending,ASCII=ASCII)
        print bcolors.GREEN+'Saved combined P(z) to COMBINED/'+\
            self.clusterName+ending+bcolors.ENDC

    def GetPickledOutput(self,pzMethod,ASCII=False):
        ''' Get P(z) for each object from the .pkl or .pklb file
            created by ProcessEAZYOutput() or ProcessBPZOutput().
            Must specify if you want the EAZY or BPZ .pkl file
            Also specify if you want the ASCII or Binary file, default
            is Binary
        '''
        pzMethod = pzMethod.upper()
        methods = ['BPZ','EAZY','COMBINED']
        if isinstance(pzMethod,str) == False:
            raise TypeError('Must enter a string for argument, must be BPZ, EAZY, or COMBINED')
        if pzMethod not in methods:
            raise RuntimeError('Must enter either \'BPZ\', \'EAZY\', or \'COMBINED\' for argument')
        if ASCII == False:
            pklFile = pzMethod+'/pzpickles/'+self.clusterName+'.pklb'
        elif ASCII == True: 
            pklFile = pzMethod+'/pzpickles/'+self.clusterName+'.pkl'
        if os.path.exists(pklFile) == False:
            print bcolors.WARNING+pklFile+' does not exist\n Making it now......'+bcolors.ENDC
            if pzMethod == 'BPZ':
                self.ProcessBPZOutput(ASCII=ASCII)
            elif pzMethod == 'EAZY':
                self.ProcessEAZYOutput(ASCII=ASCII)
            elif pzMethod == 'COMBINED':
                self.CombinePZ(ASCII=ASCII)
        pz = OpenPickle(pklFile)
        return pz

    def PlotPicklePZ(self,pzMethod,ASCII=False):
        ''' Plot P(z) from EAZY and BPZ output files.
            Gets data from the .pkl or .pklb files created 
            by ProcessBPZOutput() or ProcessEAZYOutput().
            Requires running RunEAZY() or RunBPZ().
            Specify ASCII or Binary file
        '''
        pzMethod = pzMethod.upper()
        pz = self.GetPickledOutput(pzMethod,ASCII=ASCII)
        z = pz['z'] # initialize to some range
        # to speed up matplotlib generate axis
        # then update data in axis each iteration
        fig,ax = plt.subplots()
        ax.set_xlabel('z')
        ax.set_ylabel('PDF(z)')
        line, = ax.plot(z,z) # initialize line
        #plt.show(block=False)
        pbar = pb.ProgressBar(widgets=widgets,maxval=len(pz.keys())-1).start()
        for i,key in enumerate(pz.keys()):
            if key == 'z': continue
            id = key
            data = pz[id]
            ax.set_title('%s PDF(z) %s %s'%(pzMethod,self.clusterName,id))
            line.set_ydata(data)
            lowerLim = min([0.0,1.3*min(data)])
            upperLim = 1.3*max(data)
            ax.set_ylim(lowerLim,upperLim)
            #ax.draw_artist(ax.patch)
            #ax.draw_artist(line)
            #fig.canvas.update()
            #fig.canvas.flush_events()
            saveName = 'pzplots/%s_%s_%s'%(id,self.clusterName,pzMethod)
            fig.savefig(saveName+'.png')
            cPickle.dump(ax,file(saveName+'.pklb','wb'))
            pbar.update(i)
        pbar.finish()
        print bcolors.GREEN+str(len(pz)-1)+' '+pzMethod+' PDF(z) Plots Created!'+bcolors.ENDC
  
    def PlotAllMethods(self):
        ''' Plot P(z) from all Methods on a single figure '''
        bpz = self.GetPickledOutput('BPZ')
        eazy = self.GetPickledOutput('EAZY')
        comb = self.GetPickledOutput('COMBINED')
        goodIDs, badIDs = self.CheckObjects(bpz.keys(),eazy.keys(),comb.keys())
        z = comb['z']
        fig,(ax1,ax2,ax3) = plt.subplots(1, 3, figsize = (12,6),
                                         sharex=False, sharey=False)
        ax1.set_xlabel('z');ax2.set_xlabel('z');ax3.set_xlabel('z')
        ax1.set_ylabel('PDF(z)');ax2.set_ylabel('PDF(z)');ax3.set_ylabel('PDF(z)')
        line1, = ax1.plot(bpz['z'],z)
        line2, = ax2.plot(eazy['z'],z)
        line3, = ax3.plot(comb['z'],z)
        #plt.show(block=False)
        pbar = pb.ProgressBar(widgets=widgets,maxval=len(goodIDs)-1).start()
        for i,id in enumerate(goodIDs):
            if id == 'z': continue
            ax1.set_title('BPZ PDF(z) %s %s'%(self.clusterName,id))
            line1.set_ydata(bpz[id])
            top1, bottom1 = 1.3*max(bpz[id]), min([0.0,1.3*min(bpz[id])])
            ax1.set_ylim(bottom1,top1)
            #ax1.draw_artist(ax1.patch); ax1.draw_artist(line1)
            ax2.set_title('EAZY PDF(z) %s %s'%(self.clusterName,id))
            line2.set_ydata(eazy[id])
            top2, bottom2 = 1.3*max(eazy[id]), min([0.0,1.3*min(eazy[id])])
            ax2.set_ylim(bottom2,top2)
            #ax2.draw_artist(ax2.patch); ax2.draw_artist(line2)
            ax3.set_title('COMBINED PDF(z) %s %s'%(self.clusterName,id))
            line3.set_ydata(comb[id])
            top3, bottom3 = 1.3*max(comb[id]), min([0.0,1.3*min(comb[id])])
            ax3.set_ylim(bottom3,top3)
            #ax3.draw_artist(ax3.patch); ax3.draw_artist(line3)
            #fig.canvas.update()
            #fig.canvas.flush_events()
            saveName = 'pzplots/%s_%s_ALL'%(id,self.clusterName)
            fig.savefig(saveName+'.png')
            cPickle.dump(fig,file(saveName+'.pklb','wb'))
            pbar.update(i)
        pbar.finish()
        print bcolors.GREEN+str(len(comb)-1)+' PDF(z) Triple Plots Created!'+bcolors.ENDC
            
    def Background(self,filter):
        print 'SEP Background subtraction for %s......'%filter
        data = self.dataList[filter]
        mask = (data == 0)
        bkg = sep.Background(data, mask=mask)
        bkg.subfrom(data)
        bkgRMS = bkg.globalrms
        return bkgRMS

    def ImageObjDetection(self,filter):
        print 'SEP object detection on %s......'%filter
        sigma = 3.0
        threshold = sigma * self.bkgRMS[filter]
        objs = sep.extract(self.dataList[filter], threshold, minarea=20)
        # add fields for mag and magerr
        desc = np.dtype([('RA','float64'),('DEC','float64')])
        newObjs = np.zeros(objs.shape, dtype = objs.dtype.descr + desc.descr)
        for name in objs.dtype.names:
            newObjs[name] = objs[name]
        # add RA and DEC to objs
        newObjs['RA'], newObjs['DEC'] = self.images[filter].xy_to_rd(objs['x'],objs['y'])

        return newObjs

    def AperaturePhoto(self,filter,objs):
        print 'Running Aperature Photometry on %s......'%filter
        kronrad, krflag = sep.kron_radius(self.dataList[filter],
                                          objs['x'], objs['y'],
                                          objs['a'], objs['b'],
                                          objs['theta'],
                                          6.0)

        flux, fluxerr, flag = sep.sum_ellipse(self.dataList[filter],
                                              objs['x'], objs['y'],
                                              objs['a'], objs['b'],
                                              objs['theta'],
                                              2.5*kronrad, subpix=1,
                                              err=self.bkgRMS[filter])
        #use circular aperature photometry if the kronradius is too small. see http://sep.readthedocs.org/en/v0.2.x/apertures.html
        r_min = 1.75 # minimum diameter = 3.5
        use_circle = kronrad * np.sqrt(objs['a']*objs['b']) < r_min
        cflux, cfluxerr, cflag = sep.sum_circle(self.dataList[filter],
                                                objs['x'][use_circle],
                                                objs['y'][use_circle],
                                                r_min, subpix=1,err=self.bkgRMS[filter])
        flux[use_circle] = cflux
        fluxerr[use_circle] = cfluxerr
        flag[use_circle] = cflag
        #convert flux to magnitudes using the appropriate zeropoint
        #absolute flux measurement (AB for Z-PEG)
        mag = -2.5*np.log10(flux)+self.zeroPoints[filter]
        #calculate magerr
        fluxdown = flux - fluxerr
        fluxup = flux + fluxerr
        magup = -2.5*np.log10(fluxdown) + self.zeroPoints[filter]
        magdown = -2.5*np.log10(fluxup) + self.zeroPoints[filter]
        magerr = ((magup - mag) + (mag-magdown))/2.

        return mag, magerr

    def NormalizePZ(self,pz,zrange):
        ''' I do not trust BPZ and EAZY to normalize 
            their P(z)'s to 1.0, so I'll do it myself!
            pz is a np.ndarray of the P(z) values
            zrange is the redshift sample steps
        '''
        if type(pz) == list: pz = np.array(pz)
        if type(zrange) == list: zrange = np.array(zrange)
        if isinstance(pz,np.ndarray) == False:
            raise IOError('pz must be a np.array or list, got %s'%str(type(pz)))
        if isinstance(zrange,np.ndarray) == False:
            raise IOError('zrange must be a np.array or list, got %s'%str(type(zrange)))
        splPZ = InterpolatedUnivariateSpline(zrange,pz,k=1)
        normalization = abs(1.0/np.trapz(splPZ(zrange),zrange))
        normPZ = normalization*pz[:]
        return normPZ

    def ZeroPoint(self,filter):
        print 'Finding Zero Point for %s......'%filter
        if filter == 'F105W':
            zeroPoint = 26.235
        elif filter == 'F110W':
            zeroPoint = 26.803
        elif filter == 'F125W':
            zeroPoint = 26.213
        elif filter == 'F140W':
            zeroPoint = 26.437
        elif filter == 'F160W':
            zeroPoint = 25.921
        elif filter == 'F606W':
            zeroPoint = 26.49113
        elif filter == 'F814W':
            zeroPoint = 25.1
        else:
            raise Exception("Do not have information for filter %s" %filter)
        return zeroPoint

    def CheckRequiredDirectories(self):
        '''Check for CATALOGS, BPZ, and EAZY directories
           in pwd, if they do not exist, make them'''
        nothingMade = True
        reqDirs = ['CATALOGS', 'BADCATS','pzplots',
                   'COMBINED', 'COMBINED/pzpickles',
                   'EAZY', 'EAZY/pzpickles','EAZY/'+self.clusterName,
                   'BPZ','BPZ/pzpickles']
        madeDirs = []
        for dir in reqDirs:
            if os.path.isdir(dir) == False:
                print bcolors.WARNING+'Making '+dir+' Directory'+bcolors.ENDC
                os.mkdir(dir)
                nothingMade = False; madeDirs.append(dir)
        if nothingMade == True:
            print bcolors.GREEN+'All Directories Exist\n Good to Go......'+bcolors.ENDC
        else:
            print bcolors.WARNING+'Had to create '+str(len(madeDirs))+' directory(ies)'+bcolors.ENDC
            for dir in madeDirs: print bcolors.FAIL+dir+bcolors.ENDC

    def CheckDetectionFilter(self):
        if self.detectionFilter not in self.filters:
            raise Exception("Do not have %s in images, can not use for detection filter."  %self.detectionFilter)

    def CheckObjects(self,*keyLists):
        ''' Check to make sure that the P(z) methods
            contain exactly the same objs. One
            can not have more objs than the other.
            And they must have the same obj ID#'s.
            bpz is a np.array of obj ID#'s
            eazy is a np.array of obj ID#'s
            *** Does NOT remove the 'z' key ***
        '''
        keyLists = list(keyLists)
        if len(keyLists) < 2:
            raise IOError('Require at least 2 lists of ID numbers, you only entered %d'%len(keyLists))
        for i,kL in enumerate(keyLists):
            if isinstance(kL,(np.ndarray,list)) == False:
                raise IOError('Require list or np.ndarray, not %s'%type(kL).mro()[0])
            if type(kL) == list: keyLists[i] = np.array(kL)
        # find all objects NOT in all keyLists
        diff = np.array([np.setdiff1d(kL1,kL2) for kL1 in keyLists for kL2 in keyLists if kL1 is not kL2])
        diff = np.concatenate(tuple(diffSet for diffSet in diff))
        diff = np.unique(diff)
        # find all objects in all keyLists
        same = np.array([np.intersect1d(kL1,kL2) for kL1 in keyLists for kL2 in keyLists if kL1 is not kL2])
        same = np.concatenate(tuple(sameSet for sameSet in same))
        same = np.unique(same)
        # remove any elements in same that are in diff
        same = np.delete(same,[i for i in xrange(len(same)) for j in xrange(len(diff)) if same[i]==diff[j] ])        
        if len(diff) == 0:
            print bcolors.GREEN+'All objects match!'+bcolors.ENDC
        else:
            print bcolors.WARNING+'PDF(z) methods have different objects'+bcolors.ENDC
            for id in diff: print id,
            print ''
        return same, diff

    def CheckBPZRequirements(self):
        '''check that .cat and .columns files have been created
           for <cluster Name>, so BPZ can be run'''
        cat = os.path.exists('CATALOGS/%s.cat'%self.clusterName)
        col = os.path.exists('CATALOGS/%s.columns'%self.clusterName)
        param = os.path.exists('CATALOGS/%s.pars'%self.clusterName)
        if cat == True and col == True and param == True:
            goodToGo = True
            catTime = time.ctime(os.path.getmtime('CATALOGS/%s.cat'%self.clusterName))
            colTime = time.ctime(os.path.getmtime('CATALOGS/%s.columns'%self.clusterName))
            paramTime = time.ctime(os.path.getmtime('CATALOGS/%s.pars'%self.clusterName))
            print bcolors.GREEN+'***BPZ Requirements Met***'+bcolors.ENDC
            print bcolors.BLUE+'.cat last updated on\n'+catTime+bcolors.ENDC
            print bcolors.BLUE+'.columns last updated on\n'+colTime+bcolors.ENDC
            print bcolors.BLUE+'.pars last updated on\n'+paramTime+bcolors.ENDC
        else:
            goodToGo = False
            if cat == False:
                print bcolors.FAIL+'<<<Missing .cat>>>'+bcolors.ENDC
                print bcolors.FAIL+'<<<Try running cl.WriteCatalog()>>>'+bcolors.ENDC
            if col == False:
                print bcolors.FAIL+'<<<Missing .columns for BPZ>>>'+bcolors.ENDC
                print bcolors.FAIL+'<<<Try running cl.WriteColumns()>>>'+bcolors.ENDC
            if param == False:
                print bcolors.FAIL+'<<<Missing .pars for BPZ>>>'+bcolors.ENDC
                print bcolors.FAIL+'<<<Try running cl.ParamBPZ()>>>'+bcolors.ENDC
        return goodToGo
    
    def CheckEAZYRequirements(self):
        '''check that .cat, .param, .translate, .zeropoint 
           files have been created for <cluster Name>, so 
           EAZY can be run'''
        cat = os.path.exists('CATALOGS/%s.cat'%self.clusterName)
        param = os.path.exists('CATALOGS/%s.param'%self.clusterName)
        trans = os.path.exists('CATALOGS/%s.translate'%self.clusterName)
        zp = os.path.exists('CATALOGS/%s.zeropoint'%self.clusterName)
        checks = np.array([cat,param,trans,zp])
        if np.all(checks):
            print bcolors.GREEN+'***EAZY Requirements Met***'+bcolors.ENDC
            goodToGo = True
            catTime = time.ctime(os.path.getmtime('CATALOGS/%s.cat'%self.clusterName))
            paramTime = time.ctime(os.path.getmtime('CATALOGS/%s.param'%self.clusterName))
            transTime = time.ctime(os.path.getmtime('CATALOGS/%s.translate'%self.clusterName))
            zpTime = time.ctime(os.path.getmtime('CATALOGS/%s.zeropoint'%self.clusterName))
            print bcolors.BLUE+'.cat last updated on \n'+catTime+bcolors.ENDC
            print bcolors.BLUE+'.param last updated on \n'+paramTime+bcolors.ENDC
            print bcolors.BLUE+'.translate last updated on \n'+transTime+bcolors.ENDC
            print bcolors.BLUE+'.zeropoint last updated on \n'+zpTime+bcolors.ENDC
        else:
            goodToGo = False
            if cat == False:
                print bcolors.FAIL+'<<<Missing .cat>>>'+bcolors.ENDC
                print bcolors.FAIL+'<<<Try running cl.WriteCatalog()>>>'+bcolors.ENDC
            if param == False:
                print bcolors.FAIL+'<<<Missing .param>>>'+bcolors.ENDC
                print bcolors.FAIL+'<<<Try running cl.ParamEAZY()>>>'+bcolors.ENDC
            if trans == False:
                print bcolors.FAIL+'<<<Missing .translate>>>'+bcolors.ENDC
                print bcolors.FAIL+'<<<Try running cl.TranslateEAZY()'+bcolors.ENDC
            if zp == False:
                print bcolors.FAIL+'<<<Missing .zeropoint>>>'+bcolors.ENDC
                print bcolors.FAIL+'<<<Try running cl.ZeroPointEAZY()>>>'+bcolors.ENDC
        return goodToGo

    def CheckSpecZ(self):
        if os.path.exists('specz/'+self.clusterName+'_infofile.cat') == False:
            raise Exception('Do not have SpecZ file for %s'%self.clusterName)
        fo = open(self.clusterName+'_infofile.cat','r')
        catalog = np.genfromtxt(fo,names=True)
        for obj in catalog:
            match = self.MatchRADECtoObj(obj['RA'],obj['DEC'])
            if match == None: continue
            else:
                id = match['ID']
                
    def MatchRADECtoOBJ(self,inRA,inDEC,threshold=1.0):
        ''' Given an RA and DEC, return the closest object
            found by sep. 
            threshold is how close an object must be to the
            given RA and DEC to be considered a match
        '''
        targetCoord = SkyCoord(ra=inRA*u.degree,dec=inDEC*u.degree,frame='icrs')
        c1 = SkyCoord(ra=self.objInfo['RA']*u.degree, dec=self.objInfo['DEC']*u.degree,frame='icrs')
        separation = c1.separation(targetCoord).to(u.arcsec)
        if np.sort(separation)[0] > threshold*u.arcsec:
            return None
        else:
            return self.objInfo[np.argsort(separation)[0]],np.sort(separation)[0]

    def ClusterID(self,verbose=True):
        name = self.clusterName
        if name == 'SPT0205':
            number = 10000
        elif name == 'XMM44':
            number = 20000
        elif name == 'SPT-CLJ2106-5844':
            number = 30000
        elif name == 'SPT2040':
            number = 40000
        elif name == 'SPARCS-J1049':
            number = 50000
        elif name == 'MOO-1014':
            number = 60000
        elif name == 'ISCSJ-1432+3253':
            number = 70000
        elif name == 'IDCSJ1426':
            number = 80000
        if verbose == True:
            print 'Setting Cluster ID number to: ', number
        return number

    def WriteColumns(self):
        outName = 'CATALOGS/'+self.clusterName+'.columns'
        badOutputName = 'BADCATS/'+self.clusterName+'.columns'
        output = open(outName,'w')
        badOutput = open(badOutputName,'w')
        columns = np.array(self.objInfo.dtype.names)
        filters = self.filters
        print 'Writing Columns......'
    # write a columns file, for help read $BPZPATH/test/UDFtest.columns
        output.write('# Filter    Columns    AB/Vega    zp_error    zp_offset\n')
        badOutput.write('# Filter    Columns    AB/Vega    zp_error    zp_offset\n')
        for filter in filters:
            cols = ''
            for i,name in enumerate(columns):
                if filter in name:
                    cols += '%d, '%(i+1)
            cols=cols[:-2]# trim off additional ', '
            instrument = str(self.field.filters[filters.index(filter)]['instrument'])
            detector = str(self.field.filters[filters.index(filter)]['detector'])
            fullFilterName = 'HST_'+instrument+'_'+detector+'_'+filter
            zpOff = self.ZeroPoint(filter)
            output.write('%s    %s    AB    0.01     %f\n'%(fullFilterName,cols,zpOff))
            badOutput.write('%s    %s    AB    0.01     %f\n'%(fullFilterName,cols,zpOff))
        output.write('ID        1\n')
        #m_0 should be F814W according to $BPZPATH/prior_xxxx.py
        m_0Index = np.where(columns == 'F814Wmag')[0]
        output.write('M_0        %d\n'%(m_0Index+1))
        badOutput.write('ID        1\n')
        badOutput.write('M_0        %d\n'%(m_0Index+1))
        output.close()
        badOutput.close()
        print 'Saved '+outName
        print 'Saved'+badOutputName

    def WriteCatalog(self):
        outName = 'CATALOGS/'+self.clusterName+'.cat'
        badOutputName = 'BADCATS/'+self.clusterName+'.cat'
        output = open(outName,'w')
        badOutput = open(badOutputName,'w')
        columns = np.array(self.objInfo.dtype.names)
        print 'Writing Header......'
        output.write('# ')
        badOutput.write('# ')
        for col in columns:
            output.write(col+'    ')
            badOutput.write(col+'    ')
        output.write('\n')
        badOutput.write('\n')
        print 'Writing Data......'
        # cluster.objInfo is an np.recarray (ID,x,y,RA,DEC,Filermag,Filtermagerr,...)
        nObjs = len(self.objInfo)
        for i in xrange(nObjs):
            obj = self.objInfo[i]
            mag = np.array([obj[filter+'mag'] for filter in self.filters])
            err = np.array([obj[filter+'magerr'] for filter in self.filters])
            if np.any(np.isnan(list(obj))) == True: #write data to bad catalog if it has nan
                writeTo = badOutput
            else:
                writeTo = output
            for name in columns:
                if name == 'ID':
                    writeTo.write('%d    '%obj['ID'])
                else:
                    writeTo.write('%f    '%obj[name])
            writeTo.write('\n')
        output.close()
        badOutput.close()
        print 'Saved '+outName
        print 'Saved '+badOutputName

    def TranslateEAZY(self):
        ''' Make .translate file, so EAZY can read Catalog '''
        print 'Making .translate file for EAZY......'
        fo = open('CATALOGS/'+self.clusterName+'.translate','w')
        fo.write('ID id\n')
        for filter in self.filters:
            code = self.EAZYFilterCode(filter)
            fo.write('%smag F%s\n'%(filter,code))
            fo.write('%smagerr E%s\n'%(filter,code))
            print '%smag F%s'%(filter,code)
            print '%smagerr E%s'%(filter,code)
        fo.close()
        print 'Saved CATALOGS/%s.translate'%self.clusterName

    def ZeroPointEAZY(self):
        ''' Make .zeropoint for EAZY '''
        print 'Making .zeropoint file for EAZY......'
        fo=open('CATALOGS/%s.zeropoint'%self.clusterName,'w')
        for filter in self.filters:
            code = self.EAZYFilterCode(filter)
            zp = self.ZeroPoint(filter)
            print zp
            fo.write('M%s %f\n'%(code,zp))
            print 'M%s %f'%(code,zp)
        fo.close()
        print 'Saved CATALOGS/%s.zeropoint'%self.clusterName

    def ParamBPZ(self):
        ''' Make .pars file for BPZ'''
        print 'Making .pars file for BPZ......'
        parsFile = 'CATALOGS/'+self.clusterName+'.pars'
        fo = open(parsFile,'w')
        fo.write('COLUMNS CATALOGS/%s.columns\n'%self.clusterName)
        fo.write('OUTPUT BPZ/%s.bpz\n'%self.clusterName)
        fo.write('MAG yes\n')
        fo.write('CHECK BPZ/%s.flux_comparison\n'%self.clusterName)
        fo.write('VERBOSE yes\n')
        fo.write('INTERP 2\n')
        fo.write('PROBS no\n')
        fo.write('PROBS_LITE BPZ/%s.probslite\n'%self.clusterName)
        fo.write('GET_Z yes\n')
        fo.write('PLOTS on\n')
        fo.close()
        print 'Saved %s'%fo.name

    def ParamEAZY(self,ASCII=False):
        if ASCII == False: binary = 'y'
        if ASCII == True: binary = 'n'
        ''' Make .param file for EAZY '''
        print 'Making .param file for EAZY......'
        paramFile = 'CATALOGS/'+self.clusterName+'.param'
        fo = open(paramFile,'w')
        fo.write('#### EAZY %s PARAM\n\n'%self.clusterName)
        
        fo.write('## Filters\n')
        fo.write('FILTERS_RES    FILTER.RES.latest # Filter transmission data\n')
        fo.write('FILTER_FORMAT    1 # 0: energy 1:photon-count\n')
        fo.write('SMOOTH_FILTERS    n # Smooth filter curves\n')
        fo.write('SMOOTH_SIGMA    100. #Gaussian sigma (Angstroms) to smooth filters\n\n')
        
        fo.write('## Templates\n')
        fo.write('TEMPLATES_FILE    templates/eazy_v1.2_dusty.spectra.param # Template definition file\n')
        fo.write('TEMPLATE_COMBOS    a # Template combinations\n')
        fo.write('                     # 1: one template at a time\n')
        fo.write('                     # 2: two templates, read allowed in combinations from TEMPLATES_FILE\n')
        fo.write('                     # a <or> 99: all templates simultaneously\n')
        fo.write('NMF_TOLERANCE    1.e-4 # Tolerance for non-negative combinations\n')
        fo.write('WAVELENGTH_FILE    templates/EAZY_v1.1_lines/lambda_v1.1.def # Wavelength grid definition file\n')
        fo.write('TEMP_ERR_FILE    templates/TEMPLATE_ERROR.eazy_v1.0 # Template error definition file\n')
        fo.write('TEMP_ERR_A2    0.50 # Template error amplitude\n')
        fo.write('SYS_ERR    0.00 # Systematic flux error (% of flux)\n')
        fo.write('APPLY_IGM    y # Apply Madau 1995 IGM absoption\n')
        fo.write('SCALE_2175_BUMP    0.00 # Scaling of 2175A bump. Values 0.13 (0.27) absorb ~10 (20) % at peak\n\n')
        fo.write('DUMP_TEMPLATE_CACHE    n # Write binary template cache\n')
        fo.write('USE_TEMPLATE_CACHE    n # Load in template cache\n')
        fo.write('CACHE_FILE    %s.tempfilt # Template cache file in output directory\n'%self.clusterName)
        
        fo.write('## Input Files\n')
        fo.write('CATALOG_FILE    CATALOGS/%s.cat # Catalog data file\n'%self.clusterName)
        fo.write('MAGNITUDES    y # Catalog photometry in flux or magnitudes\n')
        fo.write('NOT_OBS_THRESHOLD    -90 # Ignore flux point if < NOT_OBS_THRESHOLD\n')
        fo.write('N_MIN_COLORS    %d # Number of available filters\n\n'%(len(self.filters)-2))
        
        fo.write('##Output Files \n')
        fo.write('OUTPUT_DIRECTORY    EAZY/%s/ # where to put output\n'%self.clusterName)
        fo.write('MAIN_OUTPUT_FILE    %s # Main output file, .zout\n'%self.clusterName)
        fo.write('PRINT_ERRORS    y # Print 68, 95, 99% confidence intervals\n')
        fo.write('CHI2_SCALE    1.0 # Scale ML Chi-Squared values to improve confidence intervals\n')
        fo.write('VERBOSE_LOG    y # Dump information form the run into [MAIN_OUTPUT_FILE].param\n')
        fo.write('OBS_SED_FILE    y # Write out observed SED/object, .obs_sed\n')
        fo.write('TEMP_SED_FILE    y # Write out best template fit/object, temp_sed\n')
        fo.write('POFZ_FILE    y # Write out PofZ/object, .pz\n')#need this
        fo.write('BINARY_OUTPUT    %s # SAVE output files in binary format to save space\n\n'%binary)
        
        fo.write('## Redshift / Mag prior\n')
        fo.write('APPLY_PRIOR    y # Apply apparent magnitude prior\n')
        fo.write('PRIOR_FILE    templates/prior_K_extend.dat # File containing prior grid \n')
        fo.write('PRIOR_FILTER    217 # Filter from FILTER_RES corresponding to the columns in PRIOR_FILE\n')#maybe edit to choose the correct prior filter(217=F814W)
        fo.write('PRIOR_ABZP    25.1 # AB zeropoint of fluxes in catalog\n\n')
        
        fo.write('## Redshift Grid\n')
        fo.write('FIX_ZSPEC    n # Fix redshift to catalog zspec\n')
        fo.write('Z_MIN    0.01 # Minimum redshift\n')
        fo.write('Z_MAX    10.0 # Maximum redshift\n')
        fo.write('Z_STEP    0.01 # Redshift step size\n')
        fo.write('Z_STEP_TYPE    0 # 0=ZSTEP, 1=Z_STEP*(1+z)\n\n')
        
        fo.write('## Zeropoint Offsets\n')
        fo.write('GET_ZP_OFFSETS    y # Look for .zeropoint file and compute zeropoint offsets\n')
        fo.write('ZP_OFFSET_TOL    1.e-4 # Tolerance for iterative fit for zeropoint offsets\n\n')
        fo.write('## Rest-frame colors\n')
        fo.write('REST_FILTERS    --- # Comma separated list of rest frame filters to compute\n')
        fo.write('RF_PADDING    1000. # Padding (Ang) for choosing observed filters around specified rest-frame pair.\n')
        fo.write('RF_ERRORS    n # Compute RF color errors from p(z)\n')
        fo.write('Z_COLUMN    z_peak # Redshift to use for rest-frame color calculation (z_a,z_p,z_m1,z_m2,z_peak)\n')
        fo.write('USE_ZSPEC_FOR_REST    y # Use zspec when available for rest-frame colors\n')
        fo.write('READ_ZBIN    n # Get redshifts from OUTPUT_DIRECTORY/MAIN_OUTPUT_FILE.zbin rather than fitting them.\n\n')
        
        fo.write('## Cosmology\n')
        fo.write('H0    70.0 # Hubble Constant (km/s/Mpc)\n')
        fo.write('OMEGA_M    0.3 # Omega_matter\n')
        fo.write('OMEGA_L    0.7 # Omega_lambda\n')
        fo.close()
        print 'Saved %s'%fo.name

    def EAZYFilterCode(self,filter):
        pathToFilterInfo = 'FILTER.RES.latest.info'
        fo = open(pathToFilterInfo,'r')
        lines = fo.readlines()
        for f in self.field.filters:
            detector = str(f['detector']).lower()
            instrument = str(f['instrument']).lower()
            band = str(f['filter']).lower()
            if filter.lower() != band: continue
            for line in lines:
                line = line.lower()
                if 'hst' not in line: continue
                if band not in line: continue
                if instrument not in line: continue
                if detector not in line: continue
                components = line.split()
                code = components[0]
                return code

    def SetDetectionFilter(self,filter):
        print 'Setting detection filter to: ',filter
        print 'Reinitializing ClusterData()......'
        self.__init__(self.clusterName,detectionFilter = filter)

def WritePickle(dict,file,ASCII=False):
    ''' Takes in dictionary and file name(with path)
        and saves it as a .pkl or .pklb file
    '''
    if ASCII == True:
        file = file.rsplit('.',1)[0]+'.pkl'
        f = open(file,'w')
        cPickle.dump(dict, f, protocol = 0)
    elif ASCII == False:
        file = file.rsplit('.',1)[0]+'.pklb'
        f = open(file,'wb')
        cPickle.dump(dict, f, protocol = -1)
    f.close()
    return f.name

def OpenPickle(file):
    ''' Open file(with path), either ASCII[.pkl] 
        or Binary[.pklb], and return the dictionary
    '''
    # check if file is binary or ASCII
    if 'pklb' in file.split('.')[-1]:
        # Binary file
        fo = open(file,'rb')
        dict = cPickle.load(fo)
        print 'Binary unpickled: '+file
    elif 'pkl' in file.split('.')[-1]:
        # ASCII file
        fo = open(file,'r')
        dict = cPickle.load(fo)
        print 'ASCII unpickled: '+file
    else:
        raise IOError('File %s not understood, require a .pkl or .pklb extension'%file)
    fo.close()
    return dict
#if __name__ == '__main__':
#    cl = ClusterData('SPT0205')
#    cl.ParamEAZY()
#    cl.TranslateEAZY()
#    cl.ZeroPointEAZY()
#    cl.RunEAZY()
