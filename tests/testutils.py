# -*- coding: utf-8 -*-
"""
Created on 2020.08.01

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
#%% test objective analysis method
import numpy as np
import numba as nb
import xarray as xr


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
        lonG = xr.DataArray(lonG, dims='lon', coords={'lon': lonG})
    if isinstance(latG, list):
        latG = xr.DataArray(latG, dims='lat', coords={'lat': latG})
    
    lonS = np.deg2rad(lonS)
    latS = np.deg2rad(latS)
    lonG = np.deg2rad(lonG)
    latG = np.deg2rad(latG)
    
    dataG = latG  + lonG  # allocate one 2D slice
    dataG = dataG - dataG + dataS.mean(dimS) # initialize to mean value
    
    yc, xc = dataG.shape[-2:]
    sc     = len(dataS[dimS])
    
    radRs = (latG[-1] - latG[0]) / (yc - 1) * np.array(rads)
    
    stntag, stndis, stnwei, grdtag, grddis, grdwei = \
        __cWeights(dataS.values, lonS.values, latS.values,
                   lonG.values, latG.values, radRs.values)
    
    stntag = np.array(stntag, dtype=list).reshape([sc,])
    stndis = np.array(stndis, dtype=list).reshape([sc,])
    stnwei = np.array(stnwei, dtype=list).reshape([sc,])

    grdtag = np.array(grdtag, dtype=list).reshape([yc,xc])
    grddis = np.array(grddis, dtype=list).reshape([yc,xc])
    grdwei = np.array(grdwei, dtype=list).reshape([yc,xc])
    
    return stntag, stndis, stnwei, grdtag, grddis, grdwei


"""
Helper (private) methods are defined below
"""
@nb.jit(nopython=True, cache=False)
def __cWeights(dataS, lonS, latS, lonG, latG, rad):
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
    dataGrid: xarray.DataArray
        Result of the objective analysis
    """
    yc, xc = len(latG), len(lonG)
    sc     = len(dataS)
    
    stntag, stndis, stnwei = [], [], []
    grdtag, grddis, grdwei = [], [], []
    
    # find grids
    for s in range(sc):
        tmp1 = [[0, 0]] # with a default one for numba to determine the type
        tmp2 = [np.nan]
        tmp3 = [np.nan]
        ########## find grid points near to a station ##########
        for j in range(yc):
            for i in range(xc):
                if np.abs(lonS[s]-lonG[i])<=rad*2 and \
                    np.abs(latS[s]-latG[j])<=rad:
                    sdis = __geodist(lonG[i], lonS[s], latG[j], latS[s])
                    
                    if sdis<=rad:
                        # print('{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}'
                        #       .format(np.rad2deg(sdis),
                        #               np.rad2deg(rad),
                        #               np.rad2deg(lonG[i]), np.rad2deg(latG[j]),
                        #               np.rad2deg(lonS[s]), np.rad2deg(latS[s])))
                        
                        rad2 = rad  ** 2.0
                        dis2 = sdis ** 2.0
                        tmp1.append([j, i])
                        tmp2.append(sdis)
                        tmp3.append((rad2-dis2)/(rad2+dis2))
        
        stntag.append(tmp1)
        stndis.append(tmp2)
        stnwei.append(tmp3)
    
    # find stations
    for j in range(yc):
        for i in range(xc):
            tmp1 = [0]  # with a default one for numba to determine the type
            tmp2 = [np.nan]
            tmp3 = [np.nan]
            ########## find stations near to a grid point ##########
            for s in range(sc):
                if np.abs(lonS[s]-lonG[i])<=rad*2 and \
                    np.abs(latS[s]-latG[j])<=rad:
                    sdis = __geodist(lonG[i], lonS[s], latG[j], latS[s])
                    
                    if sdis<=rad:
                        rad2 = rad  ** 2.0
                        dis2 = sdis ** 2.0
                        tmp1.append(s)
                        tmp2.append(sdis)
                        tmp3.append((rad2-dis2)/(rad2+dis2))
            
            grdtag.append(tmp1)
            grddis.append(tmp2)
            grdwei.append(tmp3)
    
    return stntag, stndis, stnwei, grdtag, grddis, grdwei



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



ds = xr.open_dataset('./xgrads/ctls/stationPrecip.nc')

lonV = np.linspace(73, 135, 63)
latV = np.linspace(15, 55, 41)

lonG = xr.DataArray(lonV, dims='lonG', coords={'lonG': lonV})
latG = xr.DataArray(latV, dims='latG', coords={'latG': latV})

stntag, stndis, stnwei, grdtag, grddis, grdwei =\
    oacressman(ds.precip.mean('time'), ds.lons, ds.lats, lonG, latG, rads=4)






