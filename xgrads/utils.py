# -*- coding: utf-8 -*-
"""
Created on 2020.08.01

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
import xarray as xr
import numpy as np
import numba as nb
from numba.typed import List

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
            from pyproj import Transformer
            
            if 'MOAD_CEN_LAT' in ctl.comments:
                clat = float(ctl.comments['MOAD_CEN_LAT'])
                
                tmp_proj = ccrs.LambertConformal(
                    globe=globe,
                    central_latitude=clat,
                    central_longitude=pdef.slon,
                    standard_parallels=(pdef.Struelat, pdef.Ntruelat)
                ).proj4_init
                
                transformer = Transformer.from_crs('WGS84', tmp_proj)
                e, n = transformer.transform(pdef.latref, pdef.lonref, errcheck=True)
                
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


def oacressman(dataS, lonS, latS, dimS, lonG, latG,
               rads=[10, 7, 4, 2, 1], undef=-9.99e8, method='cressman'):
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
    dimS: str
        Name of the station dimension.
    lonG: DataArray or numpy.array or list
        A 1D array for longitudes of the grid (should be linear).
    latG: DataArray or numpy.array or list
        A 1D array for latitudes of the grid (should be linear).
    rads: float or list of floats
        Radii at which the analysis is performed.  Default rads are the same
        as those of GrADS.
    undef: float
        Undefined value is set if a grid has no station nearby.
    method: str
        Methods of calculating weights, should one of ['cressman', 'exp'].
    
    Returns
    -------
    dataGrid: xarray.DataArray
        Result of the objective analysis.
    wei: xarray.DataArray
        Normalized total weights.
    """
    if isinstance(lonG, list):
        lonG = xr.DataArray(lonG, dims='lon', coords={'lon': lonG})
    if isinstance(latG, list):
        latG = xr.DataArray(latG, dims='lat', coords={'lat': latG})
    
    lonS = np.deg2rad(lonS)
    latS = np.deg2rad(latS)
    lonG = np.deg2rad(lonG)
    latG = np.deg2rad(latG)
    
    meanS = dataS.mean(dimS)
    dataG = latG  + lonG  # allocate one 2D slice
    dataG = dataG - dataG + meanS # initialize to mean value
    weisG = dataG - dataG # all set zeros
    
    radRs = ((latG[-1] - latG[0]) / (len(latG) - 1)).values * np.array(rads)
    
    tmp1, tmp2 = dataG, weisG
    for rad, RR in zip(radRs, rads):
        print(f'processing {np.rad2deg(rad):4.1f}-deg ({RR} grids) radius...')
        tmp1, tmp2 = xr.apply_ufunc(__slice, dataS, lonS, latS, lonG, latG,
                                    tmp1, tmp2, rad, undef,
                                    kwargs={'method':method},
                                    dask='allow',
                                    input_core_dims=[[dimS], [dimS], [dimS],
                                                     ['lon'], ['lat'],
                                                     ['lat','lon'],
                                                     ['lat','lon'], [], []],
                                    vectorize=True,
                                    output_core_dims=[['lat','lon'],
                                                      ['lat','lon']])
    
    re = tmp1.where(tmp1!=meanS) / len(radRs)
    re = re.where(~np.isnan(re), undef)
    
    return re.rename(dataG.name), (tmp2/tmp2.max()).rename('wei')


"""
Helper (private) methods are defined below
"""
def __slice(dataS, lonS, latS, lonG, latG, dataG, weis,
            rad, undef, method='cressman'):
    #print(dataS.shape, lonS.shape, latS.shape, lonG.shape,
    #      latG.shape, dataG.shape, rad, undef)
    if method == 'cressman':
        func_wei = __weight_cressman
    else:
        func_wei = __weight_exp
    
    stntag, stnwei, grdtag, grdwei = __cWeights(dataS, lonS, latS, lonG, latG,
                                                rad, func_wei)
    
    tmp = __interp_to_stations(dataG, stntag, stnwei, dataS, undef)
    
    return __interp_to_grids(dataS-tmp, grdtag, grdwei, dataG, weis, rad, undef)


@nb.jit(nopython=True, cache=False)
def __cWeights(dataS, lonS, latS, lonG, latG, rad, func_wei):
    r"""calculate weights and store them
    
    Parameters
    ----------
    dataS: numpy.ndarray
        A 1D DataArray representing the station data.
    lonS: numpy.ndarray
        A 1D array for longitudes of the stations.
    latS: numpy.ndarray
        A 1D array for latitudes of the stations.
    lonG: numpy.ndarray
        A 1D array for longitudes of the grid (should be linear).
    latG: numpy.ndarray
        A 1D array for latitudes of the grid (should be linear).
    rad: float
        Radius (in radian) at which the analysis is performed.
    
    Returns
    -------
    stntag: list
        Indices of grid points near a station
    stnwei: list
        Weights of grid points near a station
    grdtag: list
        Indices of stations near a grid point
    grdwei: list
        Weights of stations near a grid point
    """
    yc, xc = len(latG), len(lonG)
    sc     = len(dataS)
    
    stntag = List() # station indices within the radius
    stnwei = List() # weight coefficient of the stations
    grdtag = List() # grid indices within the radius
    grdwei = List() # weight coefficient of the grids
    
    # find grids
    for s in range(sc):
        tmp1 = [[0, 0]] # with a default one for numba to determine the type
        tmp2 = [np.nan]
        ########## find grid points near to a station ##########
        for j in range(yc):
            for i in range(xc):
                if np.abs(lonS[s]-lonG[i])<=rad*2 and \
                    np.abs(latS[s]-latG[j])<=rad:
                    sdis = __geodist(lonG[i], lonS[s], latG[j], latS[s])
                    
                    if sdis < rad:
                        rad2 = rad  ** 2.0
                        dis2 = sdis ** 2.0
                        tmp1.append([j, i])
                        tmp2.append(func_wei(rad2, dis2))
        
        stntag.append(np.array(tmp1))
        stnwei.append(np.array(tmp2))
    
    # find stations
    for j in range(yc):
        for i in range(xc):
            tmp1 = [0]  # with a default one for numba to determine the type
            tmp2 = [np.nan]
            ########## find stations near to a grid point ##########
            for s in range(sc):
                if np.abs(lonS[s]-lonG[i])<=rad*2 and \
                    np.abs(latS[s]-latG[j])<=rad:
                    sdis = __geodist(lonG[i], lonS[s], latG[j], latS[s])
                    
                    if sdis < rad:
                        rad2 = rad  ** 2.0
                        dis2 = sdis ** 2.0
                        tmp1.append(s)
                        tmp2.append(func_wei(rad2, dis2))
            
            grdtag.append(np.array(tmp1))
            grdwei.append(np.array(tmp2))
    
    return stntag, stnwei, grdtag, grdwei


@nb.jit(nopython=True, cache=False)
def __interp_to_stations(dataG, stntag, stnwei, dataS, undef=-9.99e8):
    x = dataS.shape[0]
    
    re = np.zeros(dataS.shape)
    
    for i in range(x):
        sum, wsum = 0, 0
        gc = len(stntag[i])
        
        if gc > 1:
            tags = stntag[i]
            weis = stnwei[i]
            for m in range(1, tags.shape[0]): # skip first value
                w = weis[m]
                jj, ii = tags[m]

                if dataG[jj, ii] != undef:
                    sum  += w * dataG[jj, ii]
                    wsum += w
            if wsum != 0:
                re[i] = sum / wsum
            else:
                re[i] = undef
        else:
            re[i] = undef
    
    return re


@nb.jit(nopython=True, cache=False)
def __interp_to_grids(dataS, grdtag, grdwei, dataG, weisG, rad, undef=-9.99e8):
    y, x = dataG.shape

    re = dataG.copy()
    we = weisG.copy()
    
    for j in range(y):
        for i in range(x):
            sum, wsum = 0, 0
            idx = j * x + i
            sc  = len(grdtag[idx])
            
            if sc > 1:
                tags = grdtag[idx]
                weis = grdwei[idx]
                for m in range(1, tags.shape[0]): # skip first value
                    w = weis[m]
                    t = tags[m]

                    if dataS[t] != undef:
                        sum  += w * dataS[t]
                        wsum += w
                
                if wsum != 0:
                    we[j,i] += wsum
                    if re[j,i] != undef:
                        re[j,i] += sum / wsum
                    else:
                        re[j,i] = sum / wsum
    
    return re, we


@nb.jit(nopython=True, cache=False)
def __geodist(lon1, lon2, lat1, lat2):
    """Calculate great-circle distance on a sphere of radius 1.
    
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


@nb.jit(nopython=True, cache=False)
def __weight_cressman(rad2, dis2):
    return (rad2-dis2)/(rad2+dis2)


@nb.jit(nopython=True, cache=False)
def __weight_exp(rad2, dis2):
    return np.exp(-dis2/(2.0*rad2))


