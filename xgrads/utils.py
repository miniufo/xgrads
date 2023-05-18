# -*- coding: utf-8 -*-
"""
Created on 2020.08.01

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
import xarray as xr
import numpy as np
# import numba as nb

_Rearth = 6371200

"""
Some map projection function useful for plotting
"""
def get_data_projection(ctl, globe=None, Rearth=_Rearth):
    """Get Projection

    Return the data projection indicated in PDEF for plot using cartopy.

    Parameters
    ----------
    ctl: str or CtlDescriptor
        Either a string representing a `ctl` file or a CtlDescriptor.
    globe: ccrs.Globe
        Default Globe parameter (None) in cartopy uses ellipse=WGS84.  Some
        regional numerical model like WRF use a spherical earth with a radius
        of 6370 km.  If one want to plot the data output from WRF with PDEF,
        one uses a default spherical globe with a radius of 6370000 and then
        provided this globe to this function for an accurate plot.
        
        Thanks to singledoggy at
        https://github.com/miniufo/xgrads/issues/32
        https://github.com/miniufo/xgrads/issues/37

    Returns
    -------
    re: cartopy.ccrs
        Cartopy coordinate reference system object for plotting.
    """
    import cartopy.crs as ccrs
    
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
            from pyproj import Proj, transform

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
    """Interpolate PDEF data onto lat/lon grids

    Interpolate the preprojected data onto lat/lon grids defined by ctl.

    Parameters
    ----------
    ctl: str or CtlDescriptor
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
    ctl: str or CtlDescriptor
        Either a string representing a `ctl` file or a CtlDescriptor.
    latlon: boolean
        Return lat/lon coordinates on PDEF grids if True, from PDEF grids;
        Return PDEF coordinates if False, from lat/lon grids.
    Rearth: float
        Radius of the earth, default is 6371200m, consistent with GrADS.
        But this function is usually used in WRF data, then one should use
        6370000m, consistent with WRF model.
        Thanks to singledoggy at
        https://github.com/miniufo/xgrads/issues/32
        https://github.com/miniufo/xgrads/issues/37

    Returns
    -------
    y, x: xarray.DataArray, xarray.DataArray
        A tuple of (y, x) or (lat, lon)
    """
    if isinstance(ctl, str):
        from .core import CtlDescriptor
        ctl = CtlDescriptor(file=ctl)
    
    pdef = ctl.pdef
    
    if pdef is None:
        lats, lons = xr.broadcast(ctl.ydef.samples, ctl.xdef.samples)
        
        return lats, lons
    
    else:
        from pyproj import Proj #, transform
        
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



def oacressman(dataS, lonS, latS, lonG, latG, rads=[10, 7, 4, 2, 1], undef=np.nan):
    r"""Objective analysis using Cressmen method [1]_
    
    The function tries to reproduce `oacres` in GrADS.  Note that:
        
    1. The oacres function can be quite slow to execute, depending on grid
       and station data density;
    2. The scaling of the grid must be linear in lat-lon;
    3. The Cressman Analysis scheme can be unstable if the grid density is
       substantially higher than the station data density (ie, far more
       grid points than station data points). In such cases, the analysis
       can produce extrema in the grid values that are not realistic. It
       is thus suggested that you examine the results of oacres and compare
       them to the station data to insure they meet your needs.
   

    .. [1] Cressmen G. P., 1959: An operational objective analysis system,
       Monthly Weather Review, 87, 367–374.
    
    
    Parameters
    ----------
    dataS: DataArray
        A 1D DataArray representing the station data.
    lonS: DataArray
        A 1D array for longitudes of the stations.
    latS: DataArray
        A 1D array for latitudes of the stations.
    lonG: DataArray or numpy.array or list
        A 1D array for longitudes of the grid (should be linear).
    latG: DataArray or numpy.array or list
        A 1D array for latitudes of the grid (should be linear).
    rads: float or list of floats
        Radii at which the analysis is performed.  Default rads are the same
        as those of GrADS.
    undef: float
        Undefined value is set if a grid has no station nearby.

    Returns
    -------
    dataGrid: xarray.DataArray
        Result of the objective analysis
    """
    dimS = dataS.dims[-1] # assume the rightmost dim as station
    
    if isinstance(lonG, list):
        lonG = np.deg2rad(xr.DataArray(lonG, dims='lon', coords={'lon': lonG}))
    if isinstance(latG, list):
        latG = np.deg2rad(xr.DataArray(latG, dims='lat', coords={'lat': latG}))
    
    lonS = np.deg2rad(lonS)
    latS = np.deg2rad(latS)
    
    dataG = latG  + lonG  # allocate one 2D slice
    dataG = dataG - dataG + dataS.mean(dimS) # initialize to mean value
    
    yc, xc = dataG.shape
    sc     = len(dataS[dimS])
    
    radRs = np.deg2rad((latG[-1] - latG[0]) / (yc - 1)) * np.array(rads)
    
    stntag, stndis, stnwei, grdtag, grddis, grdwei = \
        __cWeights(dataS, lonS, latS, dimS, lonG, latG, radRs)
    
    return dataG


"""
Helper (private) methods are defined below
"""
def __cWeights(dataS, lonS, latS, dimS, lonG, latG, rad, method='cressman'):
    r"""calculate weights and store them
    
    Parameters
    ----------
    dataS: numpy.ndarray
        A 1D DataArray representing the station data.
    lonS: numpy.ndarray
        A 1D array for longitudes of the stations.
    latS: numpy.ndarray
        A 1D array for latitudes of the stations.
    dimS: str
        Dimension name for station data.
    lonG: numpy.ndarray
        A 1D array for longitudes of the grid (should be linear).
    latG: numpy.ndarray
        A 1D array for latitudes of the grid (should be linear).
    rad: float
        Radius (in radian) at which the analysis is performed.

    Returns
    -------
    dataGrid: xarray.DataArray
        Result of the objective analysis
    """
    yc, xc = len(latG), len(lonG)
    sc     = len(dataS[dimS])
    
    radR = 1.0
    
    stntag = np.empty([sc], dtype=object)
    stndis = np.empty([sc], dtype=object)
    stnwei = np.empty([sc], dtype=object)
    
    grdtag = np.empty([yc, xc], dtype=object)
    grddis = np.empty([yc, xc], dtype=object)
    grdwei = np.empty([yc, xc], dtype=object)
    
    # find grids
    for s in range(sc):
        ########## find grid points near to a station ##########
        for j in range(yc):
            for i in range(xc):
                if np.abs(lonS[s]-lonG[i])<=radR*3 and \
                   np.abs(latS[s]-latG[j])<=radR:
                    sdis = __geodist(lonG[i], latG[j], lonS[s], latS[s])
                    
                    if sdis<=rad:
                        rad2 = rad  ** 2.0
                        dis2 = sdis ** 2.0
                        stntag[s].append([j, i])
                        stndis[s].append(sdis)
                        stnwei[s].append((rad2-dis2)/(rad2+dis2))
    
    # find stations
    for j in range(yc):
        for i in range(xc):
            ########## find stations near to a grid point ##########
            for s in range(sc):
                if np.abs(lonS[s]-lonG[i])<=radR*3 and \
                   np.abs(latS[s]-latG[j])<=radR:
                    sdis = __geodist(lonG[i], latG[j], lonS[s], latS[s])
                    
                    if sdis<=rad:
                        rad2 = rad  ** 2.0
                        dis2 = sdis ** 2.0
                        grdtag[j,i].append(s)
                        grddis[j,i].append(sdis)
                        grdwei[j,i].append((rad2-dis2)/(rad2+dis2))
    
    return stntag, stndis, stnwei, grdtag, grddis, grdwei






def __geodist(lon1, lon2, lat1, lat2):
    """Calculate great-circle distance on a sphere.

    Parameters
    ----------
    lon1: float
        Longitude for point 1 in radian.
    lon2: float
        Longitude for point 2 in radian.
    lat1: float
        Latitude  for point 1 in radian.
    lat2: float
        Latitude  for point 2 in radian.

    Returns
    -------
    dis: float
        Great circle distance (in radian)
    """
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat/2.0)**2.0 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2.0
    
    dis = 2.0 * np.arcsin(np.sqrt(a))
    
    return dis



