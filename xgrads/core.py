# -*- coding: utf-8 -*-
"""
Created on 2020.03.02

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
import os, sys, re
import numpy as np
from datetime import datetime
from numpy import datetime64, timedelta64


template_patterns = [
    r'%x1', r'%x3', r'%y2', r'%y4', r'%m1', r'%m2',
    r'%mc', r'%d1', r'%d2', r'%h1', r'%h2', r'%h3',
    r'%n2', r'%f2', r'%f3', r'%fn2', r'%fhn', r'%fdhn',
    r'%j3', r'%t1', r'%t2', r'%t3', r'%t4', r'%t5',
    r'%t6', r'%tm1', r'%tm2', r'%tm3', r'%tm4', r'%tm5', r'%tm6'
]

"""
Core classes are defined below
"""
class CtlDescriptor(object):
    """A class for a ctl file
    
    This class represents a descriptor file like the .ctl file for GrADS.
    It generally includes a multi-dimensional dataset on the Earth.
    
    Attributes
    ----------
    dsetPath: str
        dataset file
    descPath: str
        descriptor file
    indxPath: str
        index file (for GRIB file)
    stnmPath: str
        station map file (for station file)
    
    pdef: PDEF
        Projection-definition
    tdef: Coordinate
        time coordinate definition
    zdef: Coordinate
        z-coordinate definition
    ydef: Coordinate
        y-coordinate definition
    xdef: Coordinate
        x-coordinate definition
    vdef: list of CtlVar
        Variables definition
    edef: list of str
        Ensemble definition
    
    comments: list of str
        list of global string comments
    
    zrev: bool
        z-dimension reverse (i.e., from north to south)
    yrev: bool
        y-dimension reverse (i.e., from upper to lower levels)
    
    hasData: bool
        whether the corresponding binary data file exist
    vcount: int
        variable count
    
    dtype: numpy.dtype
        data type
    periodicX: bool
        whether xdef is periodic
    
    cal365Days: bool
        whether the calendar is always 365 days (no leap year)
    template: bool
        whether it is a template for multiple binary files
    sequential: bool
        whether it is a sequential file (Fortran style)
    byteOrder: str
        byte order, ['little', 'big'] for little-endian or big-endian
    storage: str
        storage type, '99' or '-1,20' or others
    
    totalZCount: int
        total number of horizontal slice
    zRecLength: int
        record length of a single horizontal slice
    tRecLength: int
        record length of a single time (including all variables)
    """
    def __init__(self, encoding='GBK', **kwargs):
        """Constructor

        One of the keyword argument `file` or `content` should be specified.
        
        Parameters
        ----------
        encoding: {'GBK', 'UTF-8'}, optional
            Encoding for the ctl file contents.
        file: str
            The ctl path/file name.
        content: str
            A string representation for the ctl file contents.

        Returns
        -------
        CtlDescriptor
            An object represents the ctl file
        """
        self.vcount = 0
        
        self.pdef = None
        self.tdef = None
        self.zdef = None
        self.ydef = None
        self.xdef = None
        self.vdef = None
        self.edef = None
        self.comments = {}
        
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
            
            if not '/' in abspath: # thanks to Baofeng Jiao from IAP
                abspath = './' + abspath
            
            if os.path.getsize(abspath) / (1024.0*1024.0) > 2:
                raise Exception('ctl file is too large (> 2 MB)')
            
            with open(abspath, 'r', encoding=encoding) as f:
                fileContent = f.readlines()
        
        elif kwargs.get('content'):
            abspath = None
            fileContent = kwargs['content'].splitlines()
        else:
            raise Exception('invalid key word '+
                            '(file='' or content='' is allowed)')
        
        for i, line in enumerate(fileContent):
            llower = line.lower()
            if (llower.startswith('dset') or
                llower.startswith('index') or
                llower.startswith('stnmap')
               ) and abspath is not None and '^' in line:
                dirname = os.path.dirname(abspath)
                if dirname[-1] != '/':
                    fileContent[i] = line.replace('^',
                                       os.path.dirname(abspath) + '/')
                else:
                    fileContent[i] = line.replace('^',
                                       os.path.dirname(abspath))
        
        # trim down all the spaces and tabs at the beginning and end
        fileContent = [line.strip() for line in fileContent]
        
        self.descPath=abspath
        self.parse(fileContent)
    
    def parse(self, fileContent):
        """Parse file content as a multi-line str"""
        dpath_str = None
        
        for oneline in fileContent:
            onelineL = oneline.strip().lower()
            
            if onelineL.startswith('dset'):
                dpath_str = oneline.split()[1]
            elif onelineL.startswith('index'):
                self._processIndex(oneline)
            elif onelineL.startswith('stnmap'):
                self._processStnmap(oneline)
            elif onelineL.startswith('dtype'):
                self.dtype = oneline[5:].strip()
            elif onelineL.startswith('pdef'):
                self._processPDEF(oneline)
            elif onelineL.startswith('title'):
                self.title = oneline.split()[1].strip()
            elif onelineL.startswith('undef'):
                self.undef = float(oneline.split()[1].strip())
            elif onelineL.startswith('options'):
                self._processOptions(oneline)
            elif onelineL.startswith('byteswapped'):
                self.byteOrder  = 'big' \
                    if sys.byteorder == 'little' else 'little'
            elif onelineL.startswith('xdef'):
                self._processXDef(oneline, fileContent)
            elif onelineL.startswith('ydef'):
                self._processYDef(oneline,fileContent)
            elif onelineL.startswith('zdef'):
                self._processZDef(oneline, fileContent)
            elif onelineL.startswith('tdef'):
                self._processTDef(oneline)
            elif onelineL.startswith('edef'):
                self._processEDef(oneline, fileContent)
            elif onelineL.startswith('vars'):
                self._processVars(oneline, fileContent)
            elif onelineL.startswith('@ global string comment'):
                self._processGlobalComments(oneline)
            elif onelineL.startswith('*') or oneline == '':
                continue
        
        if dpath_str == None:
            raise Exception('no valid dset is parsed')
        
        if self.template:
            self._processDSets(dpath_str)
        else:
            self._processDSet(dpath_str)
        
        if self.yrev:
            self.ydef.samples = np.flip(self.ydef.samples)
        
        if self.zrev:
            self.zdef = np.flip(self.zdef)
        
        if self.edef:
            strPos = 0
            for i, e in enumerate(self.edef):
                self.edef[i] = self.edef[i]._replace(strPos=strPos)
                strPos += e.tcount * self.tRecLength
    
    def _processDSets(self, dpath_str):
        strPos = dpath_str.find('%')
        
        if strPos == -1:
            raise Exception('template is used in ctl but no % in dset')
        
        matches = find_patterns(dpath_str)
        fmtO = ''.join(matches)
        
        tokens = [self._get_template_format(token) for token in matches]
        fmtN = ''.join(tokens)
        
        fileList = []
        
        times = self.tdef.samples
        base  = self._get_field(times[0])
        for l in range(len(times)):
            part = times[l].astype('datetime64[s]').item().strftime(fmtN)
            
            fname = dpath_str.replace(fmtO, part)
            fname = self._replace_forecast_template(fname, l, base)
            
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
        lineLower = oneline.lower()
        
        if 'yrev'             in lineLower: self.yrev       = True
        if 'zrev'             in lineLower: self.zrev       = True
        if 'template'         in lineLower: self.template   = True
        if 'sequential'       in lineLower: self.sequential = True
        if '365_day_calendar' in lineLower: self.cal365Days = True
        if 'big_endian'       in lineLower: self.byteOrder  = 'big'
        if 'byteswapped'      in lineLower: self.byteOrder  = \
            'big' if sys.byteorder == 'little' else 'little'
    
    def _processXDef(self, oneline, fileContent):
        tokens = oneline.split()
        xnum   = int(tokens[1])
        
        if   tokens[2].lower() == 'linear': xlnr = True
        elif tokens[2].lower() == 'levels': xlnr = False
        else: raise Exception('invalid type for xdef (linear or levels): ' +
                              tokens[2])

        if xlnr:
            start, intv = float(tokens[3]), float(tokens[4])
            self.xdef = Coordinate('xdef', np.linspace(start,
                                                       start + intv * (xnum-1),
                                                       xnum, dtype=np.float32))
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
        
        if   tokens[2].lower() == 'linear': ylnr = True
        elif tokens[2].lower() == 'levels': ylnr = False
        else: raise Exception('invalid type for ydef (linear or levels): ' +
                              tokens[2])
        
        if ylnr:
            start, intv = float(tokens[3]), float(tokens[4])
            self.ydef = Coordinate('ydef', np.linspace(start,
                                                       start + intv * (ynum-1),
                                                       ynum, dtype=np.float32))
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
        
        if   tokens[2].lower() == 'linear': zlnr = True
        elif tokens[2].lower() == 'levels': zlnr = False
        else: raise Exception('invalid type for zdef (linear or levels): ' +
                              tokens[2])
        
        if zlnr:
            start, intv = float(tokens[3]), float(tokens[4])
            self.zdef = Coordinate('zdef', np.linspace(start,
                                                       start + intv * (znum-1),
                                                       znum, dtype=np.float32))
            
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
                raise Exception('zdef not parsed correctly')
            
            self.zdef = Coordinate('zdef', np.array(values))

    def _processTDef(self, oneline):
        tokens = oneline.split()
        tnum   = int(tokens[1])

        if tokens[2].lower() != 'linear':
            raise Exception('nonlinear tdef is not supported')

        times = self._times_to_array(tokens[3].lower(), tokens[4].lower(), tnum)
        
        self.incre = GrADS_increment_to_timedelta64(tokens[4].lower())
        self.tdef  = Coordinate('tdef', times)

    def _processEDef(self, oneline, fileContent):
        if not self.tdef:
            raise Exception('edef should be after tdef')
        
        from collections import namedtuple
        
        Ensemble = namedtuple('Ensemble', ['name', 'tcount', 'tstart',
                                           'codes', 'strPos'])
        
        tokens = oneline.split()
        enum   = int(tokens[1])
        tdef   = self.tdef
        index  = fileContent.index(oneline) + 1
        
        if len(tokens) > 2:
            if tokens[2].lower() == 'names':
                NAMES = True
            else:
                raise Exception(f'Expected str \'NAMEs\', but found {tokens[2]}')
            
            if NAMES:
                enames = [i.lower() for i in tokens[3:]]
                count  = len(enames)
                
                while count < enum:
                    split = fileContent[index].split()
                    enames += [v for v in split]
                    count  += len(split)
                    index  += 1
                
                if count != enum:
                    raise Exception((f'edef not parsed correctly, count={count} '+
                                     'while enum={enum}'))
                
                self.edef = [Ensemble(name, tdef.length(), tdef.samples[0], None, 0)
                             for name in enames]
        else:
            tmp = []
            for i in range(enum):
                eline = fileContent[index + i].strip()
                
                if eline.lower() == 'endedef':
                    raise Exception(f'not enough ensembles, need {enum}')
                
                tkns  = eline.split()
                
                name = tkns[0]
                tcnt = int(tkns[1])
                tstr = GrADStime_to_datetime(tkns[2])
                code = None
                
                if len(tkns) >3:
                    code = tkns[3]
                    
                if tcnt > tdef.length():
                    raise Exception('Tcount in EDEF must be less than or equal' +
                                    'to the Tcount in TDEF')
                
                tmp.append(Ensemble(name, tcnt, tstr, code, 0))
            
            self.edef = tmp
    
    def _processVars(self, oneline, fileContent):
        if (self.dtype != 'station' and 
            not all([self.tdef, self.zdef, self.ydef, self.xdef])):
            raise Exception('vdef should be after x, y, z and t definitions')
        
        t = self.tdef.length()
        
        if self.dtype != 'station':
            y = self.ydef.length() if self.pdef is None else self.pdef.jsize
            x = self.xdef.length() if self.pdef is None else self.pdef.isize
        
            self.zRecLength = x * y * 4
            
            # add two bytes at the beginning and the end
            if self.sequential:
                self.zRecLength += 8
        
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
            
            # if v.storage in ['99', '0']:
            if v.storage in ['99', '0', '00', '000', '1', '11', '111']:
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
        
        if fileContent[start + vnum].strip().lower() != 'endvars':
            raise Exception('endvars is expected')
    
    def _processGlobalComments(self, oneline):
        cnt = oneline[24:].strip().split('=')
        
        self.comments[cnt[0].strip()] = cnt[1].strip()
        
    def _get_template_format(self, part):
        """Get time format string
        
        See the following URL for reference:
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
        str
            The format in python datetime
        """
        for c in template_patterns:
            if c in part:
                length = len(c)
                fmt = part[:length] # format in template_patterns
                rem = part[length:] # remaining str in part
                break
            else:
                Exception('unsupported format: ' + part)
        
        if   fmt == '%y2':
            return '%y' + rem
        elif fmt == '%y4':
            return '%Y' + rem
        elif fmt == '%m1':
            return '%m' + rem
        elif fmt == '%m2':
            return '%m' + rem
        elif fmt == '%mc':
            return '%b' + rem
        elif fmt == '%d1':
            return '%d' + rem
        elif fmt == '%d2':
            return '%d' + rem
        elif fmt == '%h1':
            return '%H' + rem
        elif fmt == '%h2':
            return '%H' + rem
        elif fmt == '%n2':
            return '%M' + rem
        elif fmt in ['%f3', '%f2', '%fn2', '%fhn', '%fdhn']:
            # this is not supported by strftime()
            return '_miniufo_' + part[1:] + rem
        else:
            raise Exception('unsupported format: ' + part)
    
    def _replace_forecast_template(self, fname, l, base):
        """Replace forecast str %f as a template in dset
    
        Parameters
        ----------
        fname: str
            A given string of binary file.
        l: int
            Index of file in a template.
        base: numpy.datetime64
            Base time.
    
        Returns
        -------
        str
            A string after replacing the %f template.
        """
        amount = l * self.incre
        
        days = base[..., 2]
        hour = base[..., 3]
        mins = base[..., 4]
        
        if fname.find('_miniufo_f3') != -1:
            dt = (amount.astype("timedelta64[h]") + hour).astype(int)
            fname = fname.replace('_miniufo_f3', '{0:03d}'.format(dt))
        
        if fname.find('_miniufo_f2') != -1:
            dt = (amount.astype("timedelta64[h]") + hour).astype(int)
            fname = fname.replace('_miniufo_f2', '{0:02d}'.format(dt))
        
        if fname.find('_miniufo_fn2') != -1:
            dt = (amount.astype("timedelta64[m]") + mins).astype(int)
            fname = fname.replace('_miniufo_fn2', '{0:02d}'.format(dt))
        
        if fname.find('_miniufo_fhn') != -1:
            dt_h = (amount.astype("timedelta64[h]") + hour).astype(int)
            dt_m = (amount.astype("timedelta64[m]") + mins).astype(int) % 60
            fname = fname.replace('_miniufo_fhn',
                                  '{0:02d}'.format(dt_h)+'{0:02d}'.format(dt_m))
        
        if fname.find('_miniufo_fdhn') != -1:
            dt_d = (amount.astype("timedelta64[D]") + days).astype(int)
            dt_h = (amount.astype("timedelta64[h]") + hour).astype(int) % 24
            dt_m = (amount.astype("timedelta64[m]") + mins).astype(int) % 60
            fname = fname.replace('_miniufo_fdhn',
                                  '{0:02d}'.format(dt_d)+
                                  '{0:02d}'.format(dt_h)+
                                  '{0:02d}'.format(dt_m))
        
        return fname

    def _get_field(self, datetime64):
        """Get fields of a datetime64
        
        Convert array of datetime64 to a calendar array of
        [year, month, day, hour, minute, seconds, microsecond]
        with these quantites indexed on the last axis.
    
        Parameters
        ----------
        datetime64: datetime64 array (...)
            numpy.ndarray of datetimes of arbitrary shape
    
        Returns
        -------
        cal: uint32 array (..., 7)
            calendar array with last axis representing year, month, day, hour,
            minute, second, microsecond
        """
        # allocate output 
        out = np.empty(datetime64.shape + (7,), dtype="u4")
        
        # decompose calendar floors
        Y, M, D, h, m, s = [datetime64.astype(f'M8[{x}]') for x in 'YMDhms']
        
        out[..., 0] = Y + 1970 # Gregorian Year
        out[..., 1] = (M - Y) + 1 # month
        out[..., 2] = (D - M) + 1 # dat
        out[..., 3] = (datetime64 - D).astype('m8[h]' ) # hour
        out[..., 4] = (datetime64 - h).astype('m8[m]' ) # minute
        out[..., 5] = (datetime64 - m).astype('m8[s]' ) # second
        out[..., 6] = (datetime64 - s).astype('m8[us]') # microsecond
        
        return out


    def _times_to_array(self, strTime, incre, tnum):
        """Change format of time
        
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
        -------
        numpy array of datetime64
            Times in datetime64 format
        """
        if 'mo' in incre:
            start = GrADStime_to_datetime(strTime)
            
            lst = []
            for l in range(tnum):
                y, m = start.year, start.month
                y, m = y+int((m+l-1)/12), int((m+l-1)%12)+1
                lst.append(start.replace(year=y, month=m))
            
            return np.asarray(lst, dtype='datetime64[s]').astype('datetime64[ns]')
            
        elif 'yr' in incre:
            start = GrADStime_to_datetime(strTime)
            
            lst = []
            for l in range(tnum):
                y = start.year + l
                lst.append(start.replace(year=y))
            
            return np.asarray(lst, dtype='datetime64[s]').astype('datetime64[ns]')
        
        else:
            start = GrADStime_to_datetime64(strTime)
            intv  = GrADS_increment_to_timedelta64(incre)

            if self.cal365Days:

                lst = []
                while len(lst) < tnum:
                    isfeb29 = (start.item().day == 29) and (start.item().month == 2)
                    if not isfeb29:
                        lst.append(start)
                    start += intv

                return np.asarray(lst).astype('datetime64[ns]')

            else:

                return np.arange(start, start + intv * tnum, intv).astype('datetime64[ns]')

    def __repr__(self):
        """Print this class as a string"""
        vdef = np.array(self.vdef)
        pdef = self.pdef.proj if self.pdef is not None else ''
        
        if self.edef is not None:
            edef = [f'Ensem: {e.name}, {e.tcount}, {e.tstart}, {e.codes}, {e.strPos}'
                    for e in self.edef]
        else:
            edef = 'None'
        
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
            '       edef: ' + str(edef)           + '\n'\
            '       pdef: ' + str(pdef)           + '\n'\
            '       vdef: ' + str(vdef)



class PDEF(object):
    """PDEF class
    
    Parse necessary info in PDEF.
    
    Reference: http://cola.gmu.edu/grads/gadoc/pdef.html
    
    Attributes
    ----------
    isize: int
        size of native grid in x direction
    jsize: int
        size of native grid in y direction
    proj: str
        type of projection
    lonref: str
        reference longitude
    """
    def __init__(self, oneline):
        """Constructor
        
        Parameters
        ----------
        oneline: str
            The ASCII line of PDEF in ctl file.
        """
        lineLower = oneline.lower()
        
        if 'nps' in lineLower or 'sps' in lineLower:
            token = lineLower.split()
            
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
            
        elif 'lccr' in lineLower or 'lcc' in lineLower:
            token = lineLower.split()
            
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

    def __repr__(self):
        """Print this class as a string"""
        return '\n'.join(['%s: %s' % item for item in self.__dict__.items()])
    


class Coordinate(object):
    """Discrete sampled coordinate
    
    This is a simple wrapper for np.array for a coordinate.
    
    Attributes
    ----------
    isLinear: bool
        Steps are even or uneven
    isIncre: int
        Is increasing or decreasing
    name: str
        The name of the coordinate
    samples: str
        Discretized coordinate samples
    delSamples: str
        Finite difference between samples
    """
    def __init__(self, name, samples):
        """Constructor
        
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
        """Print this class as a string"""
        return str(self.samples)



class CtlVar(object):
    """A simple variable class used in .ctl file
    
    Attributes
    ----------
    tcount: int
        T grid points
    zcount: int
        Z grid points
    ycount: int
        Y grid points
    xcount: int
        X grid points
    undef: float
        Undefined values
    dependZ: bool
        Whether the var depends on z
    unit: str
        Unit of the variable
    name: str
        Name of the variable
    comment: str
        A short comment
    index: str
        Index of this variable
    strPos: str
        Start position (in bytes) of this variable in the binary file
    """
    __reBlank = re.compile(r'[\s\t]+')
    __reUnits = re.compile(r'\([^\(\)]+?\)')
    
    
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
        
        if len(CtlVar.__reBlank.split(oneLineStr.strip(), maxsplit=3))== 3:
            self.name, self.zcount, self.storage = \
                CtlVar.__reBlank.split(oneLineStr.strip(), maxsplit=3)
            self.comment = self.name
        else:
            self.name, self.zcount, self.storage, self.comment = \
                CtlVar.__reBlank.split(oneLineStr.strip(), maxsplit=3)
        
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
    """Convert GrADS time string e.g., 00:00z01Jan2000 to datetime
    
    Parameters
    ----------
    gradsTime: str
        Grads time in str format e.g., 00:00z01Jan2000.
    
    Returns
    -------
    datetime
        GrADS time in datetime format
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
    """Convert GrADS time string e.g., 00:00z01Jan2000 to numpy.datetime64
    
    Parameters
    ----------
    gradsTime : str
        Grads time in str format e.g., 00:00z01Jan2000.
    
    Returns
    -------
    datetime64
        GrADS time in datetime64 format
    """
    time = GrADStime_to_datetime(gradsTime)
    
    return datetime64(time.strftime('%Y-%m-%dT%H:%M:%S'))


def GrADS_increment_to_timedelta64(incre):
    """Convert GrADS time increment string to numpy.timedelta64
    
    Parameters
    ----------
    incre: str
        Grads time increment in str format e.g., 1dy.
    
    Returns
    -------
    timedelta64
        GrADS time in datetime64 format
    """
    unit   = incre[-2:].lower()
    amount = incre[:-2]

    unitDict = {
        'se': 's',
        'mn': 'm',
        'hr': 'h',
        'dy': 'D',
        'mo': 'M',
        'yr': 'Y'}

    return timedelta64(int(amount), unitDict[unit])


def find_patterns(template_path):
    """Find template patterns in a given path
    
    Parameters
    ----------
    template_path: str
        A path containing template strings.
    
    Returns
    -------
    re: list
        matched strings
    """
    pattern = re.compile('|'.join(template_patterns))
    
    return pattern.findall(template_path)

