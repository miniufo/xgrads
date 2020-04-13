# -*- coding: utf-8 -*-
"""
Created on 2020.03.02

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
from xgrads.xgrads import open_CtlDataset
import xarray as xr

# dset = open_CtlDataset('D:/Data/ULFI/output/2017101712_1721.ctl')
# dset = open_CtlDataset('D:/Data/Elena/ElenaIsenTest.ctl')
# dset = open_CtlDataset('D:/Data/ERAInterim/Keff/PV/KeffUninterp.ctl')
# dset = open_CtlDataset('D:/Data/ERAInterim/Keff/PV/PV.ctl')
# dset = open_CtlDataset('D:/Data/Haima/HaimaPTOrig.ctl')
# dset = open_CtlDataset('D:/Data/MITgcm/flt/float/Stat.ctl')
dset = open_CtlDataset('./xgrads/ctls/test8.ctl')


dset = open_CtlDataset('./xgrads/ctls/test9.ctl')

dset2 = xr.tutorial.open_dataset('air_temperature')

for l in range(len(dset.time)):
    xr.testing.assert_equal(dset.air[l], dset2.air[l])
