# -*- coding: utf-8 -*-
"""
Created on 2023.06.11

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
#%% test objective analysis method
import numpy as np
import xarray as xr
from xgrads import oacressman


def test_oacressman():
    ds = xr.open_dataset('./ctls/stationPrecip.nc')
    
    lonV = np.linspace(70, 140, 141)
    latV = np.linspace(15,  55, 81)
    
    lon = xr.DataArray(lonV, dims='lon', coords={'lon': lonV})
    lat = xr.DataArray(latV, dims='lat', coords={'lat': latV})
    
    precip = xr.where(np.isnan(ds.precip), 0, ds.precip).load()
    
    res1, wei1 = oacressman(precip, ds.lons, ds.lats, 'stnID',
                           lon, lat, rads=[4], method='cressman')
    res2, wei2 = oacressman(precip, ds.lons, ds.lats, 'stnID',
                           lon, lat, rads=[4], method='exp')
    
    assert res1.shape == res2.shape == (31, 81, 141)
    assert wei1.shape == wei2.shape == (31, 81, 141)
    
    


