# -*- coding: utf-8 -*-
"""
Created on 2020.08.01

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""

#%%
from xgrads.xgrads import open_CtlDataset


dset, ctl = open_CtlDataset('d:/EMI_2014_monthly.ctl', returnctl=True)

data = dset.emi_index[0]
data.where(data!=ctl.undef).plot(figsize=(9,5), cmap='jet')


#%%

import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
from xgrads.xgrads.utils import get_data_projection

# this data projection is defined by PDEF, and will
# be used by cartopy for plotting.
data_proj = get_data_projection(ctl)

# Note that data projection is uniquely defined by PDEF.
# But we can plot the data in different map projections.
# Here choose three for demonstration.
map_proj_pcr = ccrs.PlateCarree(central_longitude=105)
map_proj_lcc = ccrs.LambertConformal(central_longitude=105)
map_proj_orth = ccrs.Orthographic(central_longitude=105)

plt.figure(figsize=(18,10))
ax = plt.subplot(131, projection = map_proj_pcr)
ax.contourf(data.x, data.y, data, transform=data_proj, cmap='jet')
ax.coastlines('50m')
ax.add_feature(cf.BORDERS)
ax.set_title('PlateCarree projection (similar to GrADS)')

ax = plt.subplot(132, projection = map_proj_lcc)
ax.contourf(data.x, data.y, data, transform=data_proj, cmap='jet')
ax.coastlines('50m')
ax.add_feature(cf.BORDERS)
ax.set_title('Lambert conformal projection')

ax = plt.subplot(133, projection = map_proj_orth)
ax.contourf(data.x, data.y, data, transform=data_proj, cmap='jet')
ax.coastlines('50m')
ax.add_feature(cf.BORDERS)
ax.set_title('Orthographic projection')
ax.set_global()


#%%
from xgrads.xgrads.utils import get_coordinates_from_PDEF, get_data_projection


lat, lon = get_coordinates_from_PDEF(ctl)
data_proj = get_data_projection(ctl)

data = dset.emi_index[0]


plt.figure(figsize=(18,10))

ax = plt.subplot(131, projection=map_proj_pcr)
ax.scatter(lon, lat, c=data, cmap='jet', transform=ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cf.BORDERS)
ax.set_title('PlateCarree projection (similar to GrADS)')

ax = plt.subplot(132, projection = ccrs.NorthPolarStereo(central_longitude=180))
ax.scatter(lon, lat, c=data, cmap='jet', transform=ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cf.BORDERS)
ax.set_title('Lambert conformal projection')

ax = plt.subplot(133, projection = map_proj_orth)
ax.scatter(lon, lat, c=data, cmap='jet', transform=ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cf.BORDERS)
ax.set_title('Orthographic projection')
ax.set_global()


#%%
from xgrads.xgrads.utils import interp_to_latlon


emi = interp_to_latlon(dset.emi_index, ctl)

emi[0].plot(figsize=(14,7), cmap='jet')


