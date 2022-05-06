dset ^test10.dat
title CUACE_emi_index data
undef -9999.
pdef 360 320 lcc 35.00 103.50 180.50 160.50 30.00 60.00 103.50 15000. 15000.
xdef 1220 linear 62.33 .0676
ydef 682 linear 11.18 .0676
zdef 1 levels 1
tdef 1 linear JAN2019 1mo
vars 1
emi_index 0 99 pm u2/m3
endvars
