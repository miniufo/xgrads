# -*- coding: utf-8 -*-
"""
Created on 2020.03.02

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
import numpy as np
import xarray as xr
import sys
from xgrads import open_CtlDataset, open_mfdataset


def test_template():
    dset1 = open_CtlDataset('./ctls/test8.ctl')
    dset2 = open_CtlDataset('./ctls/test9.ctl')
    dset3 = xr.tutorial.open_dataset('air_temperature').load().astype('>f4')
    
    use_close = True if sys.version_info[0] == 3 and sys.version_info[1]>8 else False
    
    if use_close:
        print('1')
        for l in range(len(dset1.time)):
            xr.testing.assert_allclose(dset1.air[l], dset2.air[l])
            xr.testing.assert_allclose(dset1.air[l], dset3.air[l])        
    else:
        print('2')
        for l in range(len(dset1.time)):
            xr.testing.assert_equal(dset1.air[l], dset2.air[l])
            xr.testing.assert_equal(dset1.air[l], dset3.air[l])

    
    dset1 = open_mfdataset('./ctls/test8_*.ctl', parallel=True)
    dset2 = open_CtlDataset('./ctls/test8.ctl').load()
    dset33 = xr.tutorial.open_dataset('air_temperature').load().astype('>f4')
    
    if use_close:
        print('3')
        print(dset1)
        print(dset2)
        print(dset33)
        for l in range(len(dset1.time)):
            print('dset1 and 3 ', dset1.time[l].values, dset3.time[l].values)
            print('dset1 and 33 ', dset1.time[l].values, dset33.time[l].values)
            print(dset1.air[:,0,0].values, dset3.air[:10,0,0].values, dset33.air[:10,0,0].values)
            xr.testing.assert_allclose(dset1.air[l], dset3.air[l])        
    else:
        print('4')
        for l in range(len(dset1.time)):
            xr.testing.assert_equal(dset1.air[l], dset3.air[l])
    
    print('5')
    dset1 = open_mfdataset('./ctls/test9_*.ctl', parallel=True)
    dset2 = open_CtlDataset('./ctls/test9.ctl').load()
    
    for l in range(len(dset1.time)):
        xr.testing.assert_equal(dset1.air[l], dset2.air[l])
    
    print('6')
    # test blank line in ctls
    dset1 = open_CtlDataset('./ctls/test81.ctl')
    dset2 = open_CtlDataset('./ctls/test82.ctl')
    
    assert (dset1.x == dset2.x).all()
    print('7')
    assert (dset1.y == dset2.y).all()
    print('8')
    assert (dset1.air[0] == dset2.air).all()

def test_ensemble():
    dset1 = open_CtlDataset('./ctls/ecmf_medium_T2m1.ctl')
    
    expected = np.array([2.011963 , 1.1813354, 1.1660767])
    
    # check several ensemble values
    assert np.isclose(dset1.t2[:,-1,-1,-1], expected).all()
