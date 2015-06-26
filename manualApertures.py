import numpy as np
from pythonImports import bcolors

class Apertures(object):
    ''' An helper class to better handle hand drawn apertures
        Include a dictionary of fits file images in the form
        images = {FILTER:seechange.data.fits.,...}
    '''
    def __init__(self,**args):
        self.wantedArgs = {'RA','DEC','a','b','theta'}
        self.apertures = self.FormatArgs(args)

    def Add(self,apFile = None,**newargs):
        ''' Add an aperture by hand after __init__
            Must add all wanted arguments, to check
            what they are do cl.wantedArgs
            Alternatively, write an aperture file with a header
            # RA DEC a b theta
            where each row is an aperture
        '''
        newargs = self.FormatArgs(newargs)
        if apFile != None:
            fi = np.genfromtxt(apFile,names=True)
            newargs = np.concatenate((newargs,fi))
        if self.CheckArgs(newargs) == True:
            self.apertures = np.concatenate((self.apertures,newargs))
            print bcolors.GREEN+'Added new apertures'+bcolors.ENDC
        else:
            print bcolors.FAIL+'Aperture arguments not added'+bcolors.ENDC
            print bcolors.FAIL+'Check the argument validity by passing to ' \
                  'cl.CheckArgs(args)'+bcolors.ENDC
    
    def Remove(self,RA=[],DEC=[],purge=False):
        ''' Remove any hand drawn apertures by specifying its ra and dec
            turn purge to True to delete all apertures
        '''
        if len(RA) > 0 and len(DEC) > 0:
            idx1 = np.where(self.apertures['RA'] == RA)[0]
            idx2 = np.where(self.apertures['DEC'] == DEC)[0]
            killIdx = np.intersect1d(idx1,idx2)
            if len(killIdx) == 0:
                print bcolors.FAIL+'Could not delete apertures, ' \
                    'specified RA and DEC do no agree on the same aperture'+bcolors.ENDC
            else:
                if len(killIdx) >  1:
                    print bcolors.WARNING+'Deleting %d apertures'%len(idx1)+bcolors.ENDC
                self.apertures = np.delete(self.apertures,killIdx)
                print bcolors.CYAN+'Deleted apertures'+bcolors.ENDC
                print bcolors.CYAN+'Removed apertures from objInfo'+bcolors.ENDC
        if purge == True:
            print bcolors.WARNING+'Deleting all apertures'+bcolors.ENDC
            self.apertures = np.delete(self.apertures,np.arange(0,len(self.apertures)))
            print bcolors.WARNING+'Removed all apertures from objInfo'+bcolors.ENDC
    
    def ApplyToImages(self,images):
        ''' Apply apertures to specified images '''
        pass

    def PrintApertures(self):
        ''' Print a table of manual aperture properties '''
        print '####################################'
        print 'RA        DEC        a    b    theta'
        for ap in self.apertures:
            print ap['RA'], ap['DEC'], ap['a'], ap['b'], ap['theta']
        print '####################################'

    def Save(self,output):
        ''' Save a file with the manual aperture parameters '''
        fo = open(output,'w')
        fo.write('# RA       DEC        a    b    theta\n')
        for ap in self.apertures:
            fo.write('%f    '%ap['RA'])
            fo.write('%f    '%ap['DEC'])
            fo.write('%f    '%ap['a'])
            fo.write('%f    '%ap['b'])
            fo.write('%f    \n'%ap['theta'])
        fo.close()
        print 'Saved Manual Apertures as: %s'%output

    def FormatArgs(self,args):
        ''' Take the optional arguments, delete the unwanted arguments,
            and reformat the wanted arguments to a numpy record array.
              ** This will not add any wanted, but unspecified arguments
              ** Only takes dictionaries and np.arrays as input
        '''
        if isinstance(args,(dict,np.ndarray)) == False:
            print bcolors.FAIL+'Optional input arguments not passed in correctly.'+bcolors.ENDC
            raise IOError('Did not receive a dictionary or np.array, ' \
                              'got a %s'%str(type(args)))
        elif isinstance(args,np.ndarray):
            args = ArrayToDict(args)
        args = self.DeleteUnwantedArgs(args)
        argKeys = {k for k in args.keys()}
        if len(args) == 0:
            print bcolors.CYAN+'No args for manual aperture added by hand'+bcolors.ENDC
            args = {k:np.array([],dtype='<f8') for k in self.wantedArgs}
        # check if the arguments in args are number types 
        # and make them a np.array
        if np.all([isinstance(i,(int,float,long)) for i in args.values()]) == True:
            args = {k:np.array([i],dtype='<f8') for k,i in args.items()}
        # check if the arguments in args are lists and make them a np.array
        elif np.all([isinstance(i,list) for i in args.values()]) == True:
            args = {k:np.array(i,dtype='<f8') for k,i in args.items()}
        # all args must be of the same data type
        elif len(set(map(type,args.values()))) != 1:
            print bcolors.FAIL+'Args are not of the same data type'+bcolors.ENDC
            print bcolors.FAIL+''.join(set(map(str,map(type,args.values()))))+bcolors.ENDC
            raise IOError('Args not of the same data type, '+ \
                          'args must be all number types, lists, or np arrays')
        # if the args are not any of these, they are something that can not be used
        elif np.any(isinstance(i,(int,float,long,list,np.ndarray)) == False \
                        for i in args.values()) == True:
            print bcolors.FAIL+'Did not recieve required data types for manual apertures' \
                +bcolors.ENDC
            raise IOError('Can only take number types, lists, or numpy arrays')
        # make sure all arguments have the same number of entries
        if len(set(map(len,args.values()))) != 1:
            print bcolors.FAIL+'Did not receive equal number of paramters ' \
                'for all arguments'+bcolors.ENDC
            raise IOError('RA DEC a b theta must have equal lengths')
        # turn dict into numpy record array
        args = DictToArray(args)
        return args

    def CheckArgs(self,args):
        ''' Check the args passed to the __init__ are appropriate 
            and can be used. args should be
            RA, DEC - must be inside image
            a,b - major,minor axis of ellipse
            theta - tilt of ellipse
            all args can be entered as ints, floats, lists, or numpy arrays
        '''
        goodToGo = False
        if isinstance(args,(dict,np.ndarray)):
            args = self.FormatArgs(args)
        else:
            print bcolors.FAIL+'args not received as dictionary or np.array'+bcolors.ENDC
            raise IOError('Expected dictionary or np.array,' \
                          'got a %s'%str(type(args)))
        argKeys = set(args.dtype.names)
        # if no arguments are passed in argsEmpty == True
        argsEmpty = bool(len(args)==False)
        if argsEmpty == True:
            goodToGo = True
        elif self.wantedArgs == argKeys and argsEmpty == False:
            print bcolors.GREEN+'All required arguments for Manual Aperture met'+bcolors.ENDC
            goodToGo = True
        elif len(self.wantedArgs.difference(argKeys)) != 0:
            print bcolors.FAIL+'Missing necessary arguments'+bcolors.ENDC
            print bcolors.FAIL+','.join(self.wantedArgs.difference(argKeys))+bcolors.ENDC
            raise IOError('Did not receive all required arguments')
#        # check if ra and dec in all filters
#        if np.all([image.contains_rd(r,d) for image in self.images.values() \
#                       for r,d in zip(args['RA'],args['DEC'])]) == False:
#            raise IOError('Not all RA and DEC are inside of image')
        # check if a > b for all manual apertures
        if np.all([a >= b for a,b in zip(args['a'],args['b'])]) == False:
            print bcolors.FAIL+'Input for ellipse minor axis '+ \
                ' greater than major axis!'+bcolors.ENDC
            raise IOError('Ellipse minor axis greater than major axis, '\
                              'check a,b values in manual apertures')
        return goodToGo

    def DeleteUnwantedArgs(self,args):
        ''' Deletes unwanted keys and values from a dictionary
            Also fix the case sensitivity of args
            i.e. If you enter ra=val, val will be assigned to RA
        '''
        toUpper = {k for k in args.keys() if k.islower() and (k.upper() in self.wantedArgs)}
        toLower = {k for k in args.keys() if k.isupper() and (k.lower() in self.wantedArgs)}
        toDelete = {k for k in args.keys() if k.upper() not in map(str.upper,self.wantedArgs)}
        newargs = {}
        for k,v in args.items():
            if k in toDelete: continue
            elif k in toUpper:
                newargs[k.upper()] = v
            elif k in toLower:
                newargs[k.lower()] = v
            else:
                newargs[k] = v
        if len(toDelete) != 0:
            print bcolors.WARNING+'Found unwated arguments'+bcolors.ENDC
            print bcolors.WARNING+' '.join(toDelete)+bcolors.ENDC
            print bcolors.WARNING+'Deleted all unwanted arguments'+bcolors.ENDC
        return newargs


def ArrayToDict(array):
    ''' Take a numpy record array and convert it to
        a dictionary of np arrays
    '''
    names = np.array(array.dtype.names)
    dic = {name:array[name] for name in names}
    return dic
def DictToArray(dic):
    ''' Take a dictionary and convert it to a numpy
        record array. Values in dictionary must be 
        of the same type.
    '''
    names = dic.keys()
    # check to see if the dictionary contains a list or np array
    hasArrays = np.all([isinstance(dic[name],np.ndarray) \
                            for name in names])
    hasLists = np.all([isinstance(dic[name],list) \
                           for name in names])
    hasNums = np.all([isinstance(dic[name],(int,float,long)) \
                          for name in names])
    if hasArrays == True:
        formats = [dic[name].dtype.name for name in names]
    elif hasLists == True:
        formats = [type(dic[name][0]) for name in names]
    elif hasNums == True:
        formats = [type(i) for i in dic.values()]
        dic = {k:np.array([v]) for k,v in dic.items()}
    else:
        print bcolor.FAIL+'Did not receive numbers, lists, or np.arrays '\
            'as values in dictionary'
        raise IOError('Did not recieve expected types in dictionary values')
    dtype = dict(names=names,formats=formats)
    array = np.array(zip(*dic.values()),dtype=dtype)
    return array

def AddColumnToRecArray(array,dt):
    ''' Given a numpy record array, add new
        columns with given names and data types.
        Must give a np.dtype([('name','type'),...])
        New columns will be filled with zeros.
    '''
    if not isinstance(dt,np.dtype):
        raise IOError('Did not receive a np.dtype, check doc string')
    # convert array to dictionary, add fields to 
    # dictionary, then convert back to rec array
    dic = ArrayToDict(array)
    for name,dtype in dt.descr:
        dic[name] = np.zeros(len(array),)
    ar = DictToArray(dic)
    return ar

def DeleteColumnFromRecArray(array,col):
    ''' Given a numpy record array, remove
        a specified column by name.
    '''
    # convert array to dictionary, remove fields
    # then convert back to array
    if isinstance(col,str):
        col = np.array([col])
    elif not isinstance(col,(list,np.ndarray)):
        raise IOError('Did not receive str, list, or np.array')
    dic = ArrayToDict(array)
    dic = {k:v for k,v in dic.items() if k not in col}
    ar = DictToArray(dic)
    return ar
