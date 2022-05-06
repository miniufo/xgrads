# -*- coding: utf-8 -*-
"""
Created on 2022.05.06

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
#%%
from xgrads.xgrads import open_CtlDataset, get_data_projection,\
                          get_coordinates_from_PDEF
import numpy as np

# load data
ctl_path = 'E:/OneDrive/Python/MyPack/xgrads/ctls/grid.d1.ctl'
dset, ctl = open_CtlDataset(ctl_path, returnctl=True)

Rearth = 6370000 # consistent with WRF

crs = get_data_projection(ctl, Rearth=Rearth)

#%% plot the data and compare to WRF coords
import proplot as pplt

fig, axes = pplt.subplots(nrows=1, ncols=2, figsize=(11, 6), proj=crs)

x = dset.XLONG.squeeze()
y = dset.XLAT.squeeze()

ax = axes[0]
m=ax.pcolormesh(dset.x, dset.y, dset.land.squeeze(), transform=crs)
ax.set_title('crs')
ax.colorbar(m, loc='b')


ax = axes[1]
m=ax.pcolormesh(x, y, dset.land.squeeze())
ax.set_title('LatLon')
ax.colorbar(m, loc='b')

axes.format(abc='(a)', lonlim=[80.3, 124.1], latlim=[8,50.5])



#%% get lat/lon and see the difference w.r.t WRF output
lats, lons = get_coordinates_from_PDEF(ctl, Rearth=Rearth)

XLAT, XLON = dset.XLAT.squeeze(), dset.XLONG.squeeze()

fig, axes = pplt.subplots(nrows=1, ncols=2, figsize=(11, 6))

ax = axes[0]
m=ax.pcolormesh(XLON, XLAT, np.abs(lats-XLAT))
ax.set_title('error of lats (m)')
ax.colorbar(m, loc='b')

ax = axes[1]
m=ax.pcolormesh(XLON, XLAT, np.abs(lons-XLON))
ax.set_title('error of lons (m)')
ax.colorbar(m, loc='b')

axes.format(abc='(a)')


