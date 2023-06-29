# -*- coding: utf-8 -*-
"""
Created on 2020.03.02

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
#%%
import numpy as np
import pytest
from xgrads import CtlDescriptor


def test_parse_str():
    content=\
        "dset ^Model10\n"\
        "*ccccccccccccccccc\n"\
        "title 10-deg resolution model\n"\
        "undef -9.99e8\n"\
        "xdef 36 linear   0 10\n"\
        "ydef 19 linear -90 10\n"\
        "zdef   1 linear 0 1\n"\
        "tdef   1 linear 00z01Jan2000 1dy\n"\
        "vars 1\n"\
        "test 1 99 test variable\n"\
        "endvars\n"
    
    ctl = CtlDescriptor(content=content)
    
    assert ctl.dsetPath == '^Model10'
    assert ctl.descPath == None
    assert ctl.title == '10-deg'
    assert ctl.undef == -9.99e8
    assert ctl.zrev == False
    assert ctl.yrev == False
    assert ctl.template == False
    assert ctl.periodicX == True
    assert ctl.cal365Days == False
    assert ctl.sequential == False
    assert ctl.byteOrder == 'little'
    assert (ctl.xdef.samples == np.linspace(0, 350, 36)).all()
    assert (ctl.ydef.samples == np.linspace(-90, 90, 19)).all()
    assert (ctl.zdef.samples == np.array([0.])).all()
    assert (ctl.tdef.samples == np.array(['2000-01-01T00:00:00'], dtype='datetime64[s]'))
    assert ctl.pdef == None


    content=\
        "dset D:/Data/ERAInterim/ERAWindAmplitudes.dat\n"\
        "undef -9999.0\n"\
        "title null\n"\
        "xdef 240 linear 0.0 1.5\n"\
        "ydef 121 linear -90.0 1.5\n"\
        "zdef 1 linear 1000.0 50.0\n"\
        "tdef 1 linear 00:00Z01JAN1993 24hr\n"\
        "vars 12\n"\
        "u10amp0 0 99 amplitude of harmonics of period 365.0\n"\
        "u10amp1 0 99 amplitude of harmonics of period 182.0\n"\
        "u10amp2 0 99 amplitude of harmonics of period 91.0\n"\
        "u10amp3 0 99 amplitude of harmonics of period 60.0\n"\
        "v10amp0 0 99 amplitude of harmonics of period 365.0\n"\
        "v10amp1 0 99 amplitude of harmonics of period 182.0\n"\
        "v10amp2 0 99 amplitude of harmonics of period 91.0\n"\
        "v10amp3 0 99 amplitude of harmonics of period 60.0\n"\
        "mamp0 0 99 amplitude of harmonics of period 365.0\n"\
        "mamp1 0 99 amplitude of harmonics of period 182.0\n"\
        "mamp2 0 99 amplitude of harmonics of period 91.0\n"\
        "mamp3 0 99 amplitude of harmonics of period 60.0\n"\
        "endvars\n"
    
    ctl = CtlDescriptor(content=content)
        
    assert ctl.dsetPath == 'D:/Data/ERAInterim/ERAWindAmplitudes.dat'
    assert ctl.descPath == None
    assert ctl.title == 'null'
    assert ctl.undef == -9999.0
    assert ctl.zrev == False
    assert ctl.yrev == False
    assert ctl.template == False
    assert ctl.periodicX == True
    assert ctl.cal365Days == False
    assert ctl.sequential == False
    assert ctl.byteOrder == 'little'
    assert (ctl.xdef.samples == np.linspace(0, 358.5, 240)).all()
    assert (ctl.ydef.samples == np.linspace(-90, 90, 121)).all()
    assert (ctl.zdef.samples == np.array([1000.])).all()
    assert (ctl.tdef.samples == np.array(['1993-01-01T00:00:00'], dtype='datetime64[s]'))
    assert ctl.pdef == None


def test_coding():
    with pytest.raises(UnicodeDecodeError):
        ctl = CtlDescriptor(file='./ctls/test0.ctl')
    
    ctl = CtlDescriptor(encoding='UTF-8', file='./ctls/test0.ctl')
    
    assert ctl.dsetPath == './ctls/intensity16070412.dat'


def test_all_sample_ctls():
    paths = [
        './ctls/19855pt.grb',
        './ctls/2017101712_1721.dat',
        './ctls/ElenaIsenTest.dat',
        './ctls/flux_stn.dat',
        './ctls/KeffUninterp.dat',
        './ctls/Shanshan.dat',
        './ctls/data.grads',
        ['./ctls/test8_2013010100.dat', './ctls/test8_2013010106.dat',
         './ctls/test8_2013010112.dat', './ctls/test8_2013010118.dat'],
        ['./ctls/test9_20130101.dat', './ctls/test9_20130102.dat',
         './ctls/test9_20130103.dat'],
        './ctls/test10.dat',
    ]
    for i in range(10):
        ctl = CtlDescriptor(file='./ctls/test'+str(i+1)+'.ctl')
        re = ctl.dsetPath == paths[i]
        
        if type(re) == bool:
            assert re
        else:
            assert re.all()


def test_ensemble_ctls():
    ctl1 = CtlDescriptor(file='./ctls/ecmf_medium_T2m1.ctl')
    
    assert len(ctl1.edef) == 3
    
    ens = ['c00', 'p01', 'p02']
    
    strPos = 0
    for name, en in zip(ens, ctl1.edef):
        assert name == en.name
        assert en.strPos == strPos
        strPos += en.tcount * ctl1.tRecLength
    
    ctl2 = CtlDescriptor(file='./ctls/ecmf_medium_T2m2.ctl')
    
    assert len(ctl2.edef) == 4
    
    enames = ['ensm', 'm01', 'm02', 'm03']
    tcount = [41, 40, 39, 38]
    codes  = [None, '3,1', None, '3,2']
    
    strPos = 0
    for name, tc, code, en in zip(enames, tcount, codes, ctl2.edef):
        assert name == en.name
        assert tc   == en.tcount
        assert code == en.codes
        assert en.strPos == strPos
        strPos += en.tcount * ctl1.tRecLength

