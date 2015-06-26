import sep
import numpy as np
from seechange.controller import Controller
from seechange.config import config
from seechange.model import *
from pythonImports import bcolors

class SEPObjects(object):
    def __init__(self,clusterName,detectionFilter='F105W'):
        self.detectionFilter = detectionFilter
        self.clusterName = clusterName
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
        self.zeroPoints = {filter:self.ZeroPoint(filter) for filter in self.filters}
        # perform background subtraction on all filters
        self.bkgRMS = {filter:self.Background(filter) for filter in self.filters}
        # perform sep obj detection on detection filter image
        self.objs = self.ImageObjDetection(self.detectionFilter)

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
        sigma = 2.5
        threshold = sigma * self.bkgRMS[filter]
        objs = sep.extract(self.dataList[filter], threshold, minarea=15)
        # add fields for mag and magerr
        desc = np.dtype([('RA','float64'),('DEC','float64')])
        newObjs = np.zeros(objs.shape, dtype = objs.dtype.descr + desc.descr)
        for name in objs.dtype.names:
            newObjs[name] = objs[name]
        # add RA and DEC to objs
        newObjs['RA'], newObjs['DEC'] = self.images[filter].xy_to_rd(objs['x'],objs['y'])

        return newObjs

    def CheckDetectionFilter(self):
        if self.detectionFilter not in self.filters:
            raise Exception("Do not have %s in images, can not use for detection filter."  %self.detectionFilter)

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

