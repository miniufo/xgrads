dset ^test9_%y4%m2%d2.dat
options template yrev big_endian
undef -9.99e33
xdef 53 LINEAR 200  2.5
ydef 25 LINEAR  15  2.5
zdef  1 LEVELS 1000
tdef 12 LINEAR 01JAN2013 6hr
vars 2
air  0 99 air temperature
air2 0 99 air temperature **2
endvars
