# -*- coding: utf-8 -*-
"""
Created on 2022.05.06

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
#%%
import numpy as np
from pyproj import CRS
from xgrads import open_CtlDataset, get_data_projection,\
                          get_coordinates_from_PDEF

def test_llc_proj_wrf():
    # load data
    ctl_path = './ctls/grid.d1.ctl'
    dset, ctl = open_CtlDataset(ctl_path, returnctl=True)
    
    Rearth = 6370000 # consistent with WRF
    
    # get lat/lon and see the difference w.r.t WRF output
    lats, lons = get_coordinates_from_PDEF(ctl, Rearth=Rearth)

    XLAT, XLON = dset.XLAT.squeeze(), dset.XLONG.squeeze()
    
    assert np.isclose(np.abs(lats-XLAT), 0, atol=5e-5).all()
    assert np.isclose(np.abs(lons-XLON), 0, atol=5e-5).all()
    
    
    # write proj info into NetCDF
    crs = get_data_projection(ctl, Rearth=Rearth)
    projcrs = CRS(crs.proj4_init)
    
    cf = projcrs.to_cf()
    crs_new = projcrs.from_cf(cf)
    
    assert crs_new == projcrs


