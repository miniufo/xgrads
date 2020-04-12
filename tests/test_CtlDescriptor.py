# -*- coding: utf-8 -*-
"""
Created on 2020.03.02

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""
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

for i in range(10):
    print(CtlDescriptor(file='./xgrads/ctls/test'+str(i+1)+'.ctl'))
    print('\n')
