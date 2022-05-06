# -*- coding: utf-8 -*-
"""
Created on 2020.03.02

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
#%%
from xgrads.xgrads import CtlDescriptor


content=\
    "dset ^Model10\n"\
    "*ccccccccccccccccc\n"\
    "title 10-deg resolution model\n"\
    "undef -9.99e8\n"\
    "xdef 36 linear   0 10\n"\
    "ydef 19 linear -90 10\n"\
    "zdef   1 linear 0 1\n"\
    "tdef   1 linear 00z01Jan2000 1dy\n"\
    "vars 1\n"\
    "test 1 99 test variable\n"\
    "endvars\n"

print(CtlDescriptor(content=content))

print(CtlDescriptor(encoding='UTF-8', file='./xgrads/ctls/test0.ctl'))

for i in range(10):
    print(CtlDescriptor(file='./xgrads/ctls/test'+str(i+1)+'.ctl'))
    print('\n')

#%%
from xgrads.xgrads import CtlDescriptor, open_CtlDataset

content=\
    "dset D:/Data/ERAInterim/ERAWindAmplitudes.dat\n"\
    "undef -9999.0\n"\
    "title null\n"\
    "xdef 240 linear 0.0 1.5\n"\
    "ydef 121 linear -90.0 1.5\n"\
    "zdef 1 linear 1000.0 50.0\n"\
    "tdef 1 linear 00:00Z01JAN1993 24hr\n"\
    "vars 12\n"\
    "u10amp0 0 99 amplitude of harmonics of period 365.0\n"\
    "u10amp1 0 99 amplitude of harmonics of period 182.0\n"\
    "u10amp2 0 99 amplitude of harmonics of period 91.0\n"\
    "u10amp3 0 99 amplitude of harmonics of period 60.0\n"\
    "v10amp0 0 99 amplitude of harmonics of period 365.0\n"\
    "v10amp1 0 99 amplitude of harmonics of period 182.0\n"\
    "v10amp2 0 99 amplitude of harmonics of period 91.0\n"\
    "v10amp3 0 99 amplitude of harmonics of period 60.0\n"\
    "mamp0 0 99 amplitude of harmonics of period 365.0\n"\
    "mamp1 0 99 amplitude of harmonics of period 182.0\n"\
    "mamp2 0 99 amplitude of harmonics of period 91.0\n"\
    "mamp3 0 99 amplitude of harmonics of period 60.0\n"\
    "endvars\n"

ctl = CtlDescriptor(content=content)
print(ctl)

dset = open_CtlDataset(ctl)


#%%
from xgrads.xgrads import CtlDescriptor, open_CtlDataset

# ctl = CtlDescriptor(file='d:/Z.ctl')
dset = open_CtlDataset('d:/Z.ctl')


#%%
from xgrads.xgrads import CtlDescriptor

ctl = CtlDescriptor(file='E:/OneDrive/Python/MyPack/xgrads/ctls/grid.d1.ctl')


