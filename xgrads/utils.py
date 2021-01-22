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


def interp_to_latlon(var, ctl):
    """
    Interpolate the preprojected data onto lat/lon grids defined by ctl.

    Parameters
    ----------
    ctl : str or CtlDescriptor
        Either a string representing a `ctl` file or a CtlDescriptor.
    var: DataArray
        A variable defined at preprojected coordinates and
        need to be interpolated onto lat/lon grids.

    Returns
    -------
    re: xarray.DataArray
        Interpolated variable
    """
    if isinstance(ctl, str):
        from .core import CtlDescriptor
        ctl = CtlDescriptor(file=ctl)
    
    # get preprojected coordinates
    ypos, xpos = get_coordinates_from_PDEF(ctl, latlon=False)
    
    return var.interp(dict(y=ypos, x=xpos))
    


def get_coordinates_from_PDEF(ctl, latlon=True):
    """
    Calculate coordinates based on the PDEF information.

    Parameters
    ----------
    ctl : str or CtlDescriptor
        Either a string representing a `ctl` file or a CtlDescriptor.
    latlon : boolean
        Return lat/lon coordinates or preprojected coordinates.

    Returns
    -------
    y, x: xarray.DataArray, xarray.DataArray
    """
    if isinstance(ctl, str):
        from .core import CtlDescriptor
        ctl = CtlDescriptor(file=ctl)
    
    pdef = ctl.pdef
    
    if pdef is None:
        lats, lons = xr.broadcast(ctl.ydef.samples, ctl.xdef.samples)
        
        return lats, lons
    
    else:
        PROJ = pdef.proj
        
        if   PROJ is None:
            lats, lons = xr.broadcast(ctl.ydef.samples, ctl.xdef.samples)
            
            return lats, lons
        
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
        
        if latlon:
            xmesh, ymesh = np.meshgrid(xpos, ypos)
            
            reX, reY = p(xmesh, ymesh, inverse=True)
            
            reX = xr.DataArray(reX, dims=['y', 'x'],
                               coords={'x':xpos, 'y':ypos})
            reY = xr.DataArray(reY, dims=['y', 'x'],
                               coords={'x':xpos, 'y':ypos})
            
        else:
            lats, lons = ctl.ydef.samples, ctl.xdef.samples
            xmesh, ymesh = np.meshgrid(lons, lats)
            
            reX, reY = p(xmesh, ymesh, inverse=False)
            
            reX = xr.DataArray(reX, dims=['lat', 'lon'],
                               coords={'lat':lats, 'lon':lons})
            reY = xr.DataArray(reY, dims=['lat', 'lon'],
                               coords={'lat':lats, 'lon':lons})
        
        return reY, reX


