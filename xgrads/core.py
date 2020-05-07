# -*- coding: utf-8 -*-
"""
Created on 2020.03.02

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
import os, sys, re
import numpy as np
import cartopy.crs as ccrs
from datetime import datetime
from numpy import datetime64, timedelta64


"""
Core classes are defined below
"""
class CtlDescriptor(object):
    """
    This class represents a descriptor file like the .ctl file for GrADS.
    It generally includes a multi-dimensional spherical dataset.
    
    Attributes:
        dsetPath: dataset file
        descPath: descriptor file
        indxPath: index file (for GRIB file)
        stnmPath: station map file (for station file)
        
        pdef: projection-definition
        tdef: t-definition
        zdef: z-definition
        ydef: y-definition
        xdef: x-definition
        vdef: variable-definition
        
        zrev: z-dimension reverse (i.e., from north to south)
        yrev: y-dimension reverse (i.e., from upper to lower levels)
        
        hasData: whether the corresponding binary data file exist
        vcount : variable count
        
        dtype: data type
        periodicX: whether xdef is periodic
        
        cal365Days: whether the calendar is always 365 days (no leap year)
        template  : whether it is a template for multiple binary files
        sequential: whether it is a sequential file (Fortran style)
        byteOrder : byte order, little-endian or big-endian
        storage   : storage type, '99' or '-1,20' or others
        
        totalZCount: total number of horizontal slice
        zRecLength : record length of a single horizontal slice
        tRecLength : record length of a single time (including all variables)
    """
    def __init__(self, **kwargs):
        """
        Constructor.
        
        Parameters
        ----------
        file = fileName : str
            The ctl path/file name.
        content = content: str
            A string representation for the ctl file contents.

        Returns
        ----------
        ctl : CtlDescriptor
            An object represents the ctl file
        """
        self.vcount = 0
        
        self.pdef = None
        self.tdef = None
        self.zdef = None
        self.ydef = None
        self.xdef = None
        self.vdef = None
        
        self.dsetPath = ''
        self.descPath = ''
        self.indxPath = ''
        self.stnmPath = ''
        self.storage  = ''
        self.dtype    = ''
        self.title    = ''
        self.incre    = ''
        
        self.zrev     = False
        self.yrev     = False
        self.hasData  = False
        
        self.periodicX   = False
        self.cal365Days  = False
        self.template    = False
        self.sequential  = False
        
        self.totalZCount = 0
        self.zRecLength  = 0
        self.tRecLength  = 0
        self.byteOrder   = sys.byteorder
        
        if kwargs.get('file'):
            abspath = kwargs['file']
            
            if os.path.getsize(abspath) / (1024.0*1024.0) > 2:
                raise Exception('ctl file is too large (> 2 MB)')
            
            with open(abspath, 'r') as f:
                fileContent = f.readlines()
                
                fileContent = [line.lower() for line in fileContent]
                
                for i, line in enumerate(fileContent):
                    if (line.startswith('dset') or
                        line.startswith('index') or
                        line.startswith('stnmap')) and '^' in line:
                        dirname = os.path.dirname(abspath)
                        if dirname[-1] != '/':
                            fileContent[i] = line.replace('^',
                                               os.path.dirname(abspath) + '/')
                        else:
                            fileContent[i] = line.replace('^',
                                               os.path.dirname(abspath))
            self.descPath=abspath
        
        elif kwargs.get('content'):
            fileContent = kwargs['content'].lower().splitlines()
        else:
            raise Exception('invalid key word '+
                            '(file='' or content='' is allowed)')
        
        self.parse(fileContent)
        
    def get_data_projection(self):
        """
        Return the data projection indicated in PDEF for plot using cartopy.
        """
        if self.pdef is None:
            return ccrs.PlateCarree()
        else:
            return self.pdef.get_projection()
    
    def parse(self, fileContent):
        for oneline in fileContent:
            if   oneline.startswith('dset'):
                dpath_str = oneline.split()[1]
            elif oneline.startswith('index'):
                self._processIndex(oneline)
            elif oneline.startswith('stnmap'):
                self._processStnmap(oneline)
            elif oneline.startswith('dtype'):
                self.dtype = oneline[5:].strip()
            elif oneline.startswith('pdef'):
                self._processPDEF(oneline)
            elif oneline.startswith('title'):
                self.title = oneline.split()[1].strip()
            elif oneline.startswith('undef'):
                self.undef = float(oneline.split()[1].strip())
            elif oneline.startswith('options'):
                self._processOptions(oneline)
            elif oneline.startswith('byteswapped'):
                self.byteOrder  = 'big' \
                    if sys.byteorder == 'little' else 'little'
            elif oneline.startswith('xdef'):
                self._processXDef(oneline, fileContent)
            elif oneline.startswith('ydef'):
                self._processYDef(oneline,fileContent)
            elif oneline.startswith('zdef'):
                self._processZDef(oneline, fileContent)
            elif oneline.startswith('tdef'):
                self._processTDef(oneline)
            elif oneline.startswith('vars'):
                self._processVars(oneline, fileContent)
        
        if self.template:
            self._processDSets(dpath_str)
        else:
            self._processDSet(dpath_str)
        
        if self.yrev:
            self.ydef.samples = np.flip(self.ydef.samples)
        
        if self.zrev:
            self.zdef = np.flip(self.zdef)
    
    def _processDSets(self, dpath_str):
        times = self.tdef.samples
        
        strPos = dpath_str.index('%')
        endPos = len(dpath_str) - dpath_str[::-1].index('%') + 2
        
        template = dpath_str[strPos:endPos]
        
        tokens = []
        for token in self._split_by_len(template, 3):
            tokens.append(self._get_template_format(token))
        
        fmt = ''.join(tokens)
        
        fileList = []
        
        for l in range(len(times)):
            part = times[l].item().strftime(fmt)
            
            fname = dpath_str[:strPos] + part + dpath_str[endPos:]
            
            # remove duplicated file
            if fname not in fileList:
                fileList.append(fname)
        
        self.dsetPath = np.array(fileList)
        
        has = True
        for file in fileList:
            if not os.path.exists(file):
                has = False
                break

        self.hasData = has
        
    def _processDSet(self, dpath_str):
        self.dsetPath = dpath_str
        self.hasData  = os.path.exists(self.dsetPath)
    
    def _processIndex(self, oneline):
        self.indxPath = oneline.split()[1]
    
    def _processStnmap(self, oneline):
        self.stnmPath = oneline.split()[1]
    
    def _processPDEF(self, oneline):
        self.pdef = PDEF(oneline)
    
    def _processOptions(self, oneline):
        if 'yrev'             in oneline: self.yrev       = True
        if 'zrev'             in oneline: self.zrev       = True
        if 'template'         in oneline: self.template   = True
        if 'sequential'       in oneline: self.sequential = True
        if '365_day_calendar' in oneline: self.cal365Days = True
        if 'big_endian'       in oneline: self.byteOrder  = 'big'
        if 'byteswapped'      in oneline: self.byteOrder  = \
            'big' if sys.byteorder == 'little' else 'little'
    
    def _processXDef(self, oneline, fileContent):
        tokens = oneline.split()
        xnum   = int(tokens[1])
        
        if   tokens[2] == 'linear': xlnr = True
        elif tokens[2] == 'levels': xlnr = False
        else: raise Exception('invalid type for xdef (linear or levels): ' +
                              tokens[2])

        if xlnr:
            start, intv = float(tokens[3]), float(tokens[4])
            self.xdef = Coordinate('xdef', np.arange(start,
                                                     start + intv * xnum,
                                                     intv, dtype=np.float32))
        else:
            values = [float(i) for i in tokens[3:]]
            count  = len(values)
            index  = fileContent.index(oneline) + 1
            
            while count < xnum:
                split = fileContent[index].split()
                values += [float(v) for v in split]
                count  += len(split)
                index  += 1
            
            if count != xnum:
                raise Exception('not parse xdef correctly')
            
            self.xdef = Coordinate('xdef', np.array(values))

        self.periodicX = self.xdef.isPeriodic(360)

    def _processYDef(self, oneline, fileContent):
        tokens = oneline.split()
        ynum   = int(tokens[1])
        
        if   tokens[2] == 'linear': ylnr = True
        elif tokens[2] == 'levels': ylnr = False
        else: raise Exception('invalid type for ydef (linear or levels): ' +
                              tokens[2])
        
        if ylnr:
            start, intv = float(tokens[3]), float(tokens[4])
            self.ydef = Coordinate('ydef', np.arange(start,
                                                     start + intv * ynum,
                                                     intv, dtype=np.float32))
        else:
            values = [float(i) for i in tokens[3:]]
            count  = len(values)
            index  = fileContent.index(oneline) + 1
            
            while count < ynum:
                split = fileContent[index].split()
                values += [float(v) for v in split]
                count  += len(split)
                index  += 1
            
            if count != ynum:
                raise Exception(('not parse ydef correctly, count={0} '+
                                 'while ynum={1}').format(count, ynum))
            
            self.ydef = Coordinate('ydef', np.array(values))
    
    def _processZDef(self, oneline, fileContent):
        tokens = oneline.split()
        znum   = int(tokens[1])
        
        if   tokens[2] == 'linear': zlnr = True
        elif tokens[2] == 'levels': zlnr = False
        else: raise Exception('invalid type for zdef (linear or levels): ' +
                              tokens[2])
        
        if zlnr:
            start, intv = float(tokens[3]), float(tokens[4])
            self.zdef = Coordinate('zdef', np.arange(start,
                                                     start + intv * znum,
                                                     intv, dtype=np.float32))
            
        else:
            values = [float(i) for i in tokens[3:]]
            count  = len(values)
            index  = fileContent.index(oneline) + 1
            
            while count < znum:
                split = fileContent[index].split()
                values += [float(v) for v in split]
                count  += len(split)
                index  += 1
            
            if count != znum:
                raise Exception('not parse zdef correctly')
            
            self.zdef = Coordinate('zdef', np.array(values))

    def _processTDef(self, oneline):
        tokens = oneline.split()
        tnum   = int(tokens[1])

        if tokens[2]!='linear':
            raise Exception('nonlinear tdef is not supported')

        times = self._times_to_array(tokens[3], tokens[4], tnum)
        
        self.incre = GrADS_increment_to_timedelta64(tokens[4])
        self.tdef  = Coordinate('tdef', times)
    
    def _processVars(self, oneline, fileContent):
        if (self.dtype != 'station' and 
            not all([self.tdef, self.zdef, self.ydef, self.xdef])):
            raise Exception('vdef should be after x, y, z and t definitions')
        
        t = self.tdef.length()
        
        if self.dtype != 'station':
            y = self.ydef.length() if self.pdef is None else self.pdef.jsize
            x = self.xdef.length() if self.pdef is None else self.pdef.isize
        
            self.zRecLength = x * y * 4
        
        tokens = oneline.split()
        vnum   = int(tokens[1])
        start  = fileContent.index(oneline) + 1
        
        if vnum < 1:
            raise Exception('there should be at least one CtlVar')
        
        self.vcount = vnum
        self.vdef   = [None] * vnum
        
        v = CtlVar(fileContent[start])
        v.index =0
        v.strPos=0
        v.tcount=t
        
        if self.dtype != 'station':
            v.ycount=y
            v.xcount=x
        
        self.totalZCount += v.zcount
        self.storage      = v.storage
        self.vdef[0]      = v
        
        type1 = type2 = type3 = False
        
        for i in range(1, vnum):
            v  = CtlVar(fileContent[start+i])
            vs = self.vdef[i-1]
            
            v.index  = i
            v.undef  = self.undef
            v.tcount = t
            
            if self.dtype != 'station':
                v.ycount = y
                v.xcount = x
            
            if v.storage in ['99', '0']:
                v.strPos = vs.zcount * self.zRecLength + vs.strPos
                type1 = True
            elif v.storage == '-1,20':
                v.strPos = vs.zcount * self.zRecLength * t + vs.strPos
                type2 = True
            else:
                type3 = True
            
            if not type3 and type1 == type2:
                raise Exception('storage type should be the same')
            
            self.totalZCount += v.zcount
            self.vdef[i] = v
        
        self.tRecLength = self.zRecLength * self.totalZCount
        
        if fileContent[start + vnum].strip() != 'endvars':
            raise Exception('endvars is expected')

    def _get_template_format(self, part):
        """
        Get time format string.  See the following URL for reference:
        http://cola.gmu.edu/grads/gadoc/templates.html
        
        %x1   1 digit decade
        %x3   3 digit decade
        %y2   2 digit year
        %y4   4 digit year
        %m1   1 or 2 digit month
        %m2   2 digit month (leading zero if needed)
        %mc   3 character month abbreviation
        %d1   1 or 2 digit day
        %d2   2 digit day (leading zero if needed)
        %h1   1 or 2 digit hour
        %h2   2 digit hour
        %h3   3 digit hour (e.g., 120 or 012)
        %n2   2 digit minute; leading zero if needed
        %f2   2 digit forecast hour; leading zero if needed; more digits added
              for hours >99; hour values increase indefinitely
        %f3   3 digit forecast hour; leading zeros if needed; more digits added
              for hours >999; hour values increase indefinitely
        %fn2  2 digit forecast minute; leading zero if needed; more digits added
              for minutes > 99; minute values increase indefinitely (2.0.a9+)
        %fhn  forecast time expressed in hours and minutes (hhnn) where minute
              value (nn) is always <=59 and hour value (hh) increases indefinitely.
              If hh or nn are <=9, they are padded with a 0, so they are always at
              least 2 digits; more digits added for hours >99. (2.0.a9+)
        %fdhn forecast time expressed in days, hours, and minutes (ddhhnn) where
              minute value (nn) is always <=59, hour value (hh) is always <=23 and
              day value (dd) increases indefinitely. If dd, hh, or nn are <=9, they
              are padded with a 0 so they are always at least 2 digits; more digits
              added for days >99. (2.0.a9+)
        %j3   3 digit julian day (day of year) (2.0.a7+)
        %t1   1 or 2 digit time index (file names contain number sequences that
              begin with 1 or 01) (2.0.a7+)
        %t2   2 digit time index (file names contain number sequences that begin
              with 01) (2.0.a7+)
        %t3   3 digit time index (file names contain number sequences that begin
              with 001) (2.0.a7+)
        %t4   4 digit time index (file names contain number sequences that begin
              with 0001) (2.0.a8+)
        %t5   5 digit time index (file names contain number sequences that begin
              with 00001) (2.0.a8+)
        %t6   6 digit time index (file names contain number sequences that begin
              with 000001) (2.0.a8+)
        %tm1  1 or 2 digit time index (file names contain number sequences that
              begin with 0 or 00) (2.0.a7+)
        %tm2  2 digit time index (file names contain number sequences that begin
              with 00) (2.0.a7+)
        %tm3  3 digit time index (file names contain number sequences that begin
              with 000) (2.0.a7+)
        %tm4  4 digit time index (file names contain number sequences that begin
              with 0000) (2.0.a8+)
        %tm5  5 digit time index (file names contain number sequences that begin
              with 00000) (2.0.a8+)
        %tm6  6 digit time index (file names contain number sequences that begin
              with 000000) (2.0.a8+)
    
        Parameters
        ----------
        part : str
            A string in the above format started with %.
    
        Returns
        -------
        re : str
            A string represents the format in python datetime
        """
        if   part == '%y2':
            return '%y'
        elif part == '%y4':
            return '%Y'
        elif part == '%m1':
            return '%m'
        elif part == '%m2':
            return '%m'
        elif part == '%mc':
            return '%b'
        elif part == '%d1':
            return '%d'
        elif part == '%d2':
            return '%d'
        elif part == '%h1':
            return '%H'
        elif part == '%h2':
            return '%H'
        elif part == '%n2':
            return '%M'
        else:
            raise Exception('unsupported format: ' + part)

    def _split_by_len(self, s, size):
        """
        Split a string by a given size.
    
        Parameters
        ----------
        s : str
            A given string.
        size : int
            A given size.
    
        Returns
        -------
        re : list
            A list contains the splitted strings.
        """
        chunks = len(s)
        
        return [s[i:i + size] for i in range(0, chunks, size)]


    def _times_to_array(self, strTime, incre, tnum):
        """
        Convert GrADS time string of strart time and increment
        to an array of numpy.datetime64.
        
        Parameters
        ----------
        strTime : str
            Grads start time e.g., 00:00z01Jan2000.
        incre : str
            Grads time increment in str format e.g., 1dy.
        tnum : int
            Grads time increment in str format e.g., 1dy.
        
        Returns
        ----------
            re : numpy array of datetime64
        """
        if 'mo' in incre:
            start = GrADStime_to_datetime(strTime)
            
            lst = []
            for l in range(tnum):
                y, m = start.year, start.month
                y, m = y+int((m+l-1)/12), int((m+l-1)%12)+1
                lst.append(start.replace(year=y, month=m))
            
            return np.asarray(lst, dtype='datetime64[s]')
            
        elif 'yr' in incre:
            start = GrADStime_to_datetime(strTime)
            
            lst = []
            for l in range(tnum):
                y = start.year + l
                lst.append(start.replace(year=y))
            
            return np.asarray(lst, dtype='datetime64[s]')
        
        else:
            start = GrADStime_to_datetime64(strTime)
            intv  = GrADS_increment_to_timedelta64(incre)
            
            return np.arange(start, start + intv * tnum, intv)

    def __str__(self):
        """
        Print this class as a string.
        """
        vdef = np.array(self.vdef)
        pdef = self.pdef.proj if self.pdef is not None else ''
        
        return \
            '   dsetPath: ' + str(self.dsetPath)  + '\n'\
            '   descPath: ' + str(self.descPath)  + '\n'\
            '   indxPath: ' + str(self.indxPath)  + '\n'\
            '   stnmPath: ' + str(self.stnmPath)  + '\n'\
            '      title: ' + str(self.title)     + '\n'\
            '      undef: ' + str(self.undef)     + '\n'\
            '       zrev: ' + str(self.zrev)      + '\n'\
            '       yrev: ' + str(self.yrev)      + '\n'\
            '      dtype: ' + str(self.dtype)     + '\n'\
            '   template: ' + str(self.template)  + '\n'\
            '  periodicX: ' + str(self.periodicX) + '\n'\
            ' cal365Days: ' + str(self.cal365Days)+ '\n'\
            ' sequential: ' + str(self.sequential)+ '\n'\
            '  byteOrder: ' + str(self.byteOrder) + '\n'\
            '       xdef: ' + str(self.xdef)      + '\n'\
            '       ydef: ' + str(self.ydef)      + '\n'\
            '       zdef: ' + str(self.zdef)      + '\n'\
            '       tdef: ' + str(self.tdef)      + '\n'\
            '       pdef: ' + str(pdef)           + '\n'\
            '       vdef: ' + str(vdef)



class PDEF(object):
    """
    PDEF class.  Parse necessary info in PDEF.
    
    Reference: http://cola.gmu.edu/grads/gadoc/pdef.html
    """
    def __init__(self, oneline):
        """
        Constructor.
        
        Parameters
        ----------
        oneline : str
            The ASCII line of PDEF in ctl file.
        """
        if 'nps' in oneline or 'sps' in oneline:
            token = oneline.split()
            
            if len(token) != 8:
                raise Exception('not enough tokens for PDEF, ' +
                                'expected 8 but found ' + str(len(token)))
            
            self.isize   = int  (token[1]) # size of native grid in x direction
            self.jsize   = int  (token[2]) # size of native grid in y direction
            self.proj    =      (token[3]) # type of projection
            self.ipole   = int  (token[4]) # i-coord of pole ref to ll corner
            self.jpole   = int  (token[5]) # j-coord of pole ref to ll corner
            self.lonref  = float(token[6]) # reference longitude
            self.gridinc = float(token[7]) # distance between gripoints in km
            
        elif 'lccr' in oneline or 'lcc' in oneline:
            token = oneline.split()
            
            if len(token) != 13:
                raise Exception('not enough tokens for PDEF, ' +
                                'expected 13 but found ' + str(len(token)))
            
            self.isize   = int  (token[1]) # size of native grid in x direction
            self.jsize   = int  (token[2]) # size of native grid in y direction
            self.proj    =      (token[3]) # type of projection
            self.latref  = float(token[4]) # ref latitude
            self.lonref  = float(token[5]) # ref longitude (E>0, W<0)
            self.iref    = float(token[6]) # i of ref point
            self.jref    = float(token[7]) # j of ref point
            self.Struelat= float(token[8]) # S true lat
            self.Ntruelat= float(token[9]) # N true lat
            self.slon    = float(token[10]) # standard longitude
            self.dx      = float(token[11]) # grid X increment in meters
            self.dy      = float(token[12]) # grid Y increment in meters
        
        else:
            raise Exception('not currently supported PDEF\n' + oneline)
    
    def get_projection(self):
        PROJ = self.proj
        
        if   PROJ is None:
            return ccrs.PlateCarree()
        elif PROJ in ['lcc', 'lccr']:
            return ccrs.LambertConformal(
                      central_latitude   = self.latref,
                      central_longitude  = self.lonref,
                      standard_parallels = (self.Struelat, self.Ntruelat),
                      false_easting  = self.iref * self.dx,
                      false_northing = self.jref * self.dy)
        elif PROJ == 'nps':
            return ccrs.NorthPolarStereo(
                      central_longitude = self.lonref,
                      true_scale_latitude = 60) # used by GrADS?
        elif PROJ == 'sps':
            return ccrs.SouthPolarStereo(
                      central_longitude = self.lonref,
                      true_scale_latitude = -60) # used by GrADS?

    def __str__(self):
        """
        Print this class as a string.
        """
        return '\n'.join(['%s: %s' % item for item in self.__dict__.items()])
    



class Coordinate(object):
    """
    Discrete sampled coordinate.  This is a simple wrapper
    for np.array for a coordinate.
    """
    def __init__(self, name, samples):
        """
        Constructor.
        
        Parameters
        ----------
        name : str
            The name of the coordinate.
        samples : np.array
            1D array for the discrete coordinate.
        """
        self.isLinear   = True
        self.isIncre    = True
        self.name       = name
        self.samples    = samples
        self.delSamples = None
        
        self.max = np.max(self.samples)
        self.min = np.min(self.samples)
        
        if len(samples) > 1:
            self.delSamples = np.diff(self.samples)
            
            if self.samples[-1] < self.samples[0]:
                self.isIncre=False
        else:
            self.delSamples = np.array([1])
    
    def length(self):
        return len(self.samples)
    
    def isPeriodic(self,period):
        # not physically but generally true
        if not self.isLinear: return False
        
        delta = self.delSamples[0]
        start = self.samples[-1] + delta - period
        
        if(abs((start - self.samples[0]) / delta > 1e-4)):
            return False
        
        return True
    
    def __str__(self):
        """
        Print this class as a string.
        """
        return str(self.samples)



class CtlVar(object):
    """
    A simple variable class used in .ctl file
    """
    __reBlank = re.compile('[\s\t]+')
    __reUnits = re.compile('\([^\(\)]+?\)')
    
    
    def __init__(self, oneLineStr):
        self.tcount = 0
        self.zcount = 0
        self.ycount = 0
        self.xcount = 0
        self.undef  = np.nan
        self.dependZ= True # whether the var depends on z
        
        self.unit   = ''
        self.name   = ''
        self.comment= ''
        
        self.index  = 0
        self.strPos = 0
        
        self.name, self.zcount, self.storage, self.comment = \
            CtlVar.__reBlank.split(oneLineStr, maxsplit=3)
        
        self.zcount = int(self.zcount)
        
        findMatch = CtlVar.__reUnits.findall(self.comment)
        
        if findMatch:
            self.unit    = findMatch[-1]
            self.comment = self.comment[:self.comment.index(self.unit)].strip()
        else:
            self.unit = ''

        if self.zcount == 0:
            self.zcount = 1
            self.dependZ= False
    
    def __str__(self):
        """
        Print this class as a string.
        """
        return '\n'.join(('%8s: %s' % item for item in self.__dict__.items()))
    
    def __repr__(self):
        """
        Print this class as a string.
        """
        return 'CtlVar: {0:8s} in shape (t={1:d}, z={2:d}, y={3:d}, x={4:d})'\
                .format(self.name,
                        self.tcount, self.zcount,
                        self.ycount, self.xcount)



"""
Some useful functions defined here
"""
def GrADStime_to_datetime(gradsTime):
    """
    Convert GrADS time string e.g., 00:00z01Jan2000 to datetime
    
    Parameters
    ----------
    gradsTime : str
        Grads time in str format e.g., 00:00z01Jan2000.
    
    Returns
    ----------
        re : datetime
    """
    lens = len(gradsTime)
    
    if   lens==15 or lens==14:
        time = datetime.strptime(gradsTime, "%H:%Mz%d%b%Y")
    elif lens==12 or lens==11:
        time = datetime.strptime(gradsTime, "%Hz%d%b%Y"   )
    elif lens== 9 or lens== 8:
        time = datetime.strptime(gradsTime, "%d%b%Y"      )
    elif lens== 7:
        time = datetime.strptime(gradsTime, "%b%Y"        )
    else:
        raise Exception('invalid length of GrADS date/time string')
    
    return time


def GrADStime_to_datetime64(gradsTime):
    """
    Convert GrADS time string e.g., 00:00z01Jan2000 to numpy.datetime64
    
    Parameters
    ----------
    gradsTime : str
        Grads time in str format e.g., 00:00z01Jan2000.
    
    Returns
    ----------
        re : datetime64
    """
    time = GrADStime_to_datetime(gradsTime)
    
    return datetime64(time.strftime('%Y-%m-%dT%H:%M:%S'))


def GrADS_increment_to_timedelta64(incre):
    """
    Convert GrADS time increment string to numpy.timedelta64
    
    Parameters
    ----------
    incre : str
        Grads time increment in str format e.g., 1dy.
    
    Returns
    ----------
        re : timedelta64
    """
    unit   = incre[-2:]
    amount = incre[:-2]

    unitDict = {
        'se': 's',
        'mn': 'm',
        'hr': 'h',
        'dy': 'D',
        'mo': 'M',
        'yr': 'Y'}

    return timedelta64(int(amount), unitDict[unit])



"""
Helper (private) methods are defined below
"""


