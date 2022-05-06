# -*- coding: utf-8 -*-
"""
Created on 2020.08.01

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
from pyproj import Proj, transform

_Rearth = 6371200

"""
Some map projection function useful for plotting
"""
def get_data_projection(ctl, globe=None, Rearth=_Rearth):
    """
    Return the data projection indicated in PDEF for plot using cartopy.

    Parameters
    ----------
    ctl: str or CtlDescriptor
        Either a string representing a `ctl` file or a CtlDescriptor.
    globe: ccrs.Globe
        Default Globe parameter (None) in cartopy uses ellipse=WGS84.  Some
        regional numerical model like WRF use a spherical earth with a radius
        of 6370 km.  If one want to plot the data output from WRF with PDEF,
        one uses a default globe:
            ```python
            globe = ccrs.Globe(ellipse='sphere',
                               semimajor_axis=6370000,
                               semiminor_axis=6370000)
            ```
        and then provided this globe to this function for an accurate plot.
        Thanks to singledoggy at https://github.com/miniufo/xgrads/issues/32

    Returns
    -------
    xarray.Dataset
    """
    pdef = ctl.pdef
    
    if globe == None: # default globe
        globe = ccrs.Globe(ellipse='sphere',
                           semimajor_axis=Rearth, semiminor_axis=Rearth)
    
    if pdef is None:
        return ccrs.PlateCarree(globe=globe)
    else:
        PROJ = pdef.proj
        
        if PROJ is None:
            return ccrs.PlateCarree(globe=globe)
        
        elif PROJ in ['lcc', 'lccr']:
            if 'MOAD_CEN_LAT' in ctl.comments:
                clat = float(ctl.comments['MOAD_CEN_LAT'])
                
                wgs_proj = Proj(proj='latlong', datum='WGS84')
                
                tmp_proj = Proj(ccrs.LambertConformal(
                    globe=globe,
                    central_latitude=clat,
                    central_longitude=pdef.slon,
                    standard_parallels=(pdef.Struelat, pdef.Ntruelat)
                ).proj4_init)
                
                e, n = transform(wgs_proj, tmp_proj, pdef.lonref, pdef.latref)
                
                return ccrs.LambertConformal(globe=globe,
                          central_latitude   = clat,
                          central_longitude  = pdef.slon,
                          standard_parallels = (pdef.Struelat, pdef.Ntruelat),
                          false_easting  = (pdef.iref - 1) * pdef.dx - e,
                          false_northing = (pdef.jref - 1) * pdef.dy - n)
            else:
                clat = pdef.latref
                
                return ccrs.LambertConformal(globe=globe,
                          central_latitude   = clat,
                          central_longitude  = pdef.slon,
                          standard_parallels = (pdef.Struelat, pdef.Ntruelat),
                          false_easting  = pdef.iref * pdef.dx,
                          false_northing = pdef.jref * pdef.dy)
            
        elif PROJ == 'nps':
            return ccrs.NorthPolarStereo(globe=globe,
                      central_longitude = pdef.lonref,
                      true_scale_latitude = 60) # used by GrADS?
        
        elif PROJ == 'sps':
            return ccrs.SouthPolarStereo(globe=globe,
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
    


def get_coordinates_from_PDEF(ctl, latlon=True, Rearth=_Rearth):
    """
    Calculate coordinates based on the PDEF information.

    Parameters
    ----------
    ctl : str or CtlDescriptor
        Either a string representing a `ctl` file or a CtlDescriptor.
    latlon : boolean
        Return lat/lon coordinates or preprojected coordinates.
    Rearth: float
        Radius of the earth, default is 6371200m, consistent with GrADS.
        But this function is usually used in WRF data, then one should use
        6370000m, consistent with WRF model.
        Thanks to singledoggy at https://github.com/miniufo/xgrads/issues/32

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
            # p = Proj(proj='lcc', datum='WGS84',
            #          lat_1=pdef.Struelat,
            #          lat_2=pdef.Ntruelat,
            #          lat_0=clat,
            #          lon_0=pdef.slon,
            #          x_0=pdef.iref * pdef.dx,
            #          y_0=pdef.jref * pdef.dy,
            #          a=_Rearth, b=_Rearth)
            p = Proj(get_data_projection(ctl, Rearth=Rearth).proj4_init)
            
            ypos = np.linspace(0, (pdef.jsize-1) * pdef.dx, pdef.jsize)
            xpos = np.linspace(0, (pdef.isize-1) * pdef.dy, pdef.isize)
        
        elif PROJ == 'nps':
            inc = pdef.gridinc * 1000 # change unit from km to m
            
            p = Proj(proj='stere',
                     lon_0=pdef.lonref, lat_0=90, lat_ts=60,
                     x_0=0, y_0=0,
                     a=_Rearth, b=_Rearth)
            
            ypos = np.linspace(-(pdef.jpole), (pdef.jsize-pdef.jpole),
                                 pdef.jsize) * inc
            xpos = np.linspace(-(pdef.ipole), (pdef.isize-pdef.ipole),
                                 pdef.isize) * inc
        
        elif PROJ == 'sps':
            inc = pdef.gridinc * 1000 # change unit from km to m
            
            p = Proj(proj='stere',
                     lon_0=pdef.lonref, lat_0=-90, lat_ts=-60,
                     x_0=0, y_0=0,
                     a=_Rearth, b=_Rearth)
            
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


