dset ^./Test/%y4%m2/test_%y4%m2.bin
options little_endian template
title test
undef -999.9
xdef 360 linear 0.0 1.0
ydef 180 linear -89.5 1.0
zdef 1 levels 1000.0
tdef 2 linear 1jan1994 1mo
vars 1
v 0 99 test_var
endvars
