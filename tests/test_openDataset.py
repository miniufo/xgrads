# -*- coding: utf-8 -*-
"""
Created on 2020.03.02

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
from xgrads.xgrads import open_CtlDataset, open_mfdataset
import xarray as xr

dset = open_CtlDataset('D:/Data/ULFI/output/2017101712_1721.ctl')
dset.sst[-1].where(dset.sst[-1]!=-999).plot()

#%%
dset = open_CtlDataset('D:/Data/ERAInterim/ElenaIsen/ElenaTest.ctl')
dset.ut[-1].plot()

#%%
dset = open_CtlDataset('D:/Data/ERAInterim/Keff/PV/KeffUninterp.ctl')

#%%
dset = open_CtlDataset('D:/Data/ERAInterim/Keff/PV/PV.ctl')


#%%
dset = open_CtlDataset('D:/Data/Haima/HaimaPTOrig.ctl')

#%%
dset = open_CtlDataset('D:/Data/MITgcm/flt/float/Stat.ctl')

#%%
dset1 = open_CtlDataset('./xgrads/ctls/test8.ctl')
dset2 = open_CtlDataset('./xgrads/ctls/test9.ctl')
dset3 = xr.tutorial.open_dataset('air_temperature')

for l in range(len(dset1.time)):
    xr.testing.assert_equal(dset1.air[l], dset2.air[l])
    xr.testing.assert_equal(dset1.air[l], dset3.air[l])

#%%
dset1 = open_mfdataset('./xgrads/ctls/test8_*.ctl', parallel=True)
dset2 = open_CtlDataset('./xgrads/ctls/test8.ctl').load()
dset3 = xr.tutorial.open_dataset('air_temperature').load()

for l in range(len(dset1.time)):
    xr.testing.assert_equal(dset1.air[l], dset3.air[l])


#%%
dset1 = open_mfdataset('./xgrads/ctls/test9_*.ctl', parallel=True)
dset2 = open_CtlDataset('./xgrads/ctls/test9.ctl').load()
dset3 = xr.tutorial.open_dataset('air_temperature').load()

for l in range(len(dset1.time)):
    xr.testing.assert_equal(dset1.air[l], dset2.air[l])
    
#%%
dset1, ctl = open_CtlDataset('./xgrads/ctls/test10.ctl', returnctl=True)

#%%
from xgrads.xgrads import CtlDescriptor

ctl1 = CtlDescriptor(file='./xgrads/ctls/test11.ctl')
dset1, ctl = open_CtlDataset('./xgrads/ctls/test11.ctl', returnctl=True)
