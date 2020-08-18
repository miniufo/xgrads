# -*- coding: utf-8 -*-
"""
Created on 2020.08.01

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
from pyproj import Proj


"""
Some map projection function useful for plotting
"""
def get_data_projection(ctl):
    """
    Return the data projection indicated in PDEF for plot using cartopy.

    Parameters
    ----------
    ctl : str or CtlDescriptor
        Either a string representing a `ctl` file or a CtlDescriptor.

    Returns
    -------
    xarray.Dataset
    """
    pdef = ctl.pdef
    
    if pdef is None:
        return ccrs.PlateCarree()
    else:
        PROJ = pdef.proj
        
        if   PROJ is None:
            return ccrs.PlateCarree()
        elif PROJ in ['lcc', 'lccr']:
            return ccrs.LambertConformal(
                      central_latitude   = pdef.latref,
                      central_longitude  = pdef.lonref,
                      standard_parallels = (pdef.Struelat, pdef.Ntruelat),
                      false_easting  = pdef.iref * pdef.dx,
                      false_northing = pdef.jref * pdef.dy)
        elif PROJ == 'nps':
            return ccrs.NorthPolarStereo(
                      central_longitude = pdef.lonref,
                      true_scale_latitude = 60) # used by GrADS?
        elif PROJ == 'sps':
            return ccrs.SouthPolarStereo(
                      central_longitude = pdef.lonref,
                      true_scale_latitude = -60) # used by GrADS?



def get_latlon_from_PDEF(ctl):
    """
    Calculate lat/lon coordinates of the PDEF data.

    Parameters
    ----------
    ctl : str or CtlDescriptor
        Either a string representing a `ctl` file or a CtlDescriptor.

    Returns
    -------
    xarray.Dataset
    """
    if isinstance(ctl, str):
        from .core import CtlDescriptor
        ctl = CtlDescriptor(file=ctl)
    
    pdef = ctl.pdef
    
    if pdef is None:
        lats, lons = xr.broadcast(ctl.lat, ctl.lon)
        
        return lats.ravel(), lons.ravel()
    
    else:
        PROJ = pdef.proj
        
        if   PROJ is None:
            lats, lons = xr.broadcast(ctl.lat, ctl.lon)
            
            return lats.ravel(), lons.ravel()
        
        elif PROJ in ['lcc', 'lccr']:
            p = Proj(proj='lcc', datum='WGS84',
                     lat_1=pdef.Struelat,
                     lat_2=pdef.Ntruelat,
                     lat_0=pdef.latref,
                     lon_0=pdef.lonref,
                     x_0=pdef.iref * pdef.dx,
                     y_0=pdef.jref * pdef.dy)
            
            ypos = np.linspace(0, (pdef.jsize-1) * pdef.dx, pdef.jsize)
            xpos = np.linspace(0, (pdef.isize-1) * pdef.dy, pdef.isize)
        
        elif PROJ == 'nps':
            inc = pdef.gridinc * 1000 # change unit from km to m
            
            p = Proj(proj='stere',
                     lon_0=pdef.lonref, lat_0=90, lat_ts=60,
                     x_0=0, y_0=0,
                     a=6371200, b=6371200)
            
            ypos = np.linspace(-(pdef.jpole), (pdef.jsize-pdef.jpole),
                                 pdef.jsize) * inc
            xpos = np.linspace(-(pdef.ipole), (pdef.isize-pdef.ipole),
                                 pdef.isize) * inc
        
        elif PROJ == 'sps':
            inc = pdef.gridinc * 1000 # change unit from km to m
            
            p = Proj(proj='stere',
                     lon_0=pdef.lonref, lat_0=-90, lat_ts=-60,
                     x_0=0, y_0=0,
                     a=6371200, b=6371200)
            
            ypos = np.linspace(-(pdef.jpole), (pdef.jsize-pdef.jpole),
                                 pdef.jsize) * inc
            xpos = np.linspace(-(pdef.ipole), (pdef.isize-pdef.ipole),
                                 pdef.isize) * inc
        else:
            raise Exception('unsupported projection ' + PROJ)
            
        xpos, ypos = np.meshgrid(xpos, ypos)
        
        lons, lats = p(xpos, ypos, inverse=True)
        
        return lats.ravel(), lons.ravel()


