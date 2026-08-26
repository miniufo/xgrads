## 2. Parse CTL files

------------------
### Introduction
One of the very fundamental purposes of [xgrads](https://github.com/miniufo/xgrads) is to properly parse a given CTL file.  Now [xgrads](https://github.com/miniufo/xgrads) supports most commonly-used CTL files except those with `dtype` (data type) being `station`, `GRIB`, or `NetCDF`.  These types of datasets, especially `GRIB` and `NetCDF` are all supported in [xarray](https://docs.xarray.dev/en/stable/), so we don't need extra efforts to support these kinds of datasets.  `station`-type of data is not very commonly-used now because interpolation can be easily done by [xarray](https://docs.xarray.dev/en/stable/)


------------------
### Practice
#### Standard CTL file
For a standard CTL file `ElenaIsenTest.ctl` like:
```text
dset ^ElenaIsenTest.dat
undef -9.99E8
title ERAInterim
xdef 1 linear 231.0 0.5
ydef 91 linear -9.0 0.5
zdef 34 linear 265.0 5.0
tdef 32 linear 00:00z28Aug1985 6hr
vars 7
utm        34    99 tangential velocity (m s^-1)
vrm        34    99 radial velocity (m s^-1)
aam        34    99 traditional absolute angular momentum (m^2 s^-1)
pvm        34    99 vertical component of potential voticity (m^2 K kg^-1 s^-1)
sgm        34    99 pseudo density (kg m^-2 K^-1)
mgm        34    99 Montgomery potential [m**2 s**-2] (profile)
prm        34    99 Pressure [Pa] (profile)
endvars
```

the parsing of the CTL file is as simple as:
```python
from xgrads import CtlDescriptor

ctl = CtlDescriptor(file='ElenaIsenTest.ctl')
print(ctl)
```

and the output is:
```python
>>> print(ctl)
   dsetPath: ./ElenaIsenTest.dat
   descPath: ./ElenaIsenTest.ctl
   indxPath: 
   stnmPath: 
      title: ERAInterim
      undef: -999000000.0
       zrev: False
       yrev: False
      dtype: 
   template: False
  periodicX: True
 cal365Days: False
 sequential: False
  byteOrder: little
       xdef: [231.]
       ydef: [-9.  -8.5 -8.  -7.5 -7.  -6.5 -6.  -5.5 -5.  -4.5 -4.  -3.5 -3.  -2.5
 -2.  -1.5 -1.  -0.5  0.   0.5  1.   1.5  2.   2.5  3.   3.5  4.   4.5
  5.   5.5  6.   6.5  7.   7.5  8.   8.5  9.   9.5 10.  10.5 11.  11.5
 12.  12.5 13.  13.5 14.  14.5 15.  15.5 16.  16.5 17.  17.5 18.  18.5
 19.  19.5 20.  20.5 21.  21.5 22.  22.5 23.  23.5 24.  24.5 25.  25.5
 26.  26.5 27.  27.5 28.  28.5 29.  29.5 30.  30.5 31.  31.5 32.  32.5
 33.  33.5 34.  34.5 35.  35.5 36. ]
       zdef: [265. 270. 275. 280. 285. 290. 295. 300. 305. 310. 315. 320. 325. 330.
 335. 340. 345. 350. 355. 360. 365. 370. 375. 380. 385. 390. 395. 400.
 405. 410. 415. 420. 425. 430.]
       tdef: ['1985-08-28T00:00:00' '1985-08-28T06:00:00' '1985-08-28T12:00:00'
 '1985-08-28T18:00:00' '1985-08-29T00:00:00' '1985-08-29T06:00:00'
 '1985-08-29T12:00:00' '1985-08-29T18:00:00' '1985-08-30T00:00:00'
 '1985-08-30T06:00:00' '1985-08-30T12:00:00' '1985-08-30T18:00:00'
 '1985-08-31T00:00:00' '1985-08-31T06:00:00' '1985-08-31T12:00:00'
 '1985-08-31T18:00:00' '1985-09-01T00:00:00' '1985-09-01T06:00:00'
 '1985-09-01T12:00:00' '1985-09-01T18:00:00' '1985-09-02T00:00:00'
 '1985-09-02T06:00:00' '1985-09-02T12:00:00' '1985-09-02T18:00:00'
 '1985-09-03T00:00:00' '1985-09-03T06:00:00' '1985-09-03T12:00:00'
 '1985-09-03T18:00:00' '1985-09-04T00:00:00' '1985-09-04T06:00:00'
 '1985-09-04T12:00:00' '1985-09-04T18:00:00']
       pdef: 
       vdef: [CtlVar: utm      in shape (t=32, z=34, y=91, x=1)
 CtlVar: vrm      in shape (t=32, z=34, y=91, x=1)
 CtlVar: aam      in shape (t=32, z=34, y=91, x=1)
 CtlVar: pvm      in shape (t=32, z=34, y=91, x=1)
 CtlVar: sgm      in shape (t=32, z=34, y=91, x=1)
 CtlVar: mgm      in shape (t=32, z=34, y=91, x=1)
 CtlVar: prm      in shape (t=32, z=34, y=91, x=1)]
```


#### CTL file with non-ASCII characters
For a CTL file `intensity16070412.ctl` containing Chinese or other non-ASCII characters:
```text
dset ^intensity16070412.dat
undef -9999.0
title 2016年07月04日12时ecmf模式数据
xdef    1 linear 75.0 0.25
ydef    1 linear -10.0 0.25
zdef    1 levels 1000.0
tdef   21 linear 12:00z04Jul2016 6hr
vars    2
prsNEPARTAK   0    99 time series of NEPARTAK's minimum central pressure
wndNEPARTAK   0    99 time series of NEPARTAK's maximum sustained wind speed
endvars
```

the parsing of the CTL file should add a kwarg `encoding=UTF-8` as:
```python
from xgrads import CtlDescriptor

ctl = CtlDescriptor(file='ElenaIsenTest.ctl', encoding='UTF-8')
print(ctl)
```

and the output is:
```python
>>> print(ctl)
   dsetPath: ./intensity16070412.dat
   descPath: ./intensity16070412.ctl
   indxPath: 
   stnmPath: 
      title: 2016年07月04日12时ecmf模式数据
      undef: -9999.0
       zrev: False
       yrev: False
      dtype: 
   template: False
  periodicX: True
 cal365Days: False
 sequential: False
  byteOrder: little
       xdef: [75.]
       ydef: [-10.]
       zdef: [1000.]
       tdef: ['2016-07-04T12:00:00' '2016-07-04T18:00:00' '2016-07-05T00:00:00'
 '2016-07-05T06:00:00' '2016-07-05T12:00:00' '2016-07-05T18:00:00'
 '2016-07-06T00:00:00' '2016-07-06T06:00:00' '2016-07-06T12:00:00'
 '2016-07-06T18:00:00' '2016-07-07T00:00:00' '2016-07-07T06:00:00'
 '2016-07-07T12:00:00' '2016-07-07T18:00:00' '2016-07-08T00:00:00'
 '2016-07-08T06:00:00' '2016-07-08T12:00:00' '2016-07-08T18:00:00'
 '2016-07-09T00:00:00' '2016-07-09T06:00:00' '2016-07-09T12:00:00']
       pdef: 
       vdef: [CtlVar: prsNEPARTAK in shape (t=21, z=1, y=1, x=1)
 CtlVar: wndNEPARTAK in shape (t=21, z=1, y=1, x=1)
 CtlVar: ULFI1_250 in shape (t=21, z=1, y=1, x=1)
 CtlVar: ULFI1_200 in shape (t=21, z=1, y=1, x=1)
 CtlVar: ULFI1_150 in shape (t=21, z=1, y=1, x=1)
 CtlVar: ULFI2_250 in shape (t=21, z=1, y=1, x=1)
 CtlVar: ULFI2_200 in shape (t=21, z=1, y=1, x=1)
 CtlVar: ULFI2_150 in shape (t=21, z=1, y=1, x=1)]
```

Otherwise, an error will be thrown as:
```python
UnicodeDecodeError: 'gbk' codec can't decode byte 0xb4 in position 54: illegal multibyte sequence
```

#### CTL with template
For a single CTL file `test8.ctl` containing `template` for describing multiple binary data files:
```text
dset ^test8_%y4%m2%d2%h2.dat
options template yrev big_endian
undef -9.99e33
xdef 53 LINEAR 200  2.5
ydef 25 LINEAR  15  2.5
zdef  1 LEVELS 1000
tdef  4 LINEAR 01JAN2013 6hr
vars 1
air 0 99 air temperature
endvars
```

[xgrads](https://github.com/miniufo/xgrads) will automatically recognize multiple data files that satisfying the pattern of the `dset` and output as:
```python
>>> print(ctl)
   dsetPath: ['./test8_2013010100.dat'
 './test8_2013010106.dat'
 './test8_2013010112.dat'
 './test8_2013010118.dat']
   descPath: ./test8.ctl
   indxPath: 
   stnmPath: 
      title: 
      undef: -9.99e+33
       zrev: False
       yrev: True
      dtype: 
   template: True
  periodicX: True
 cal365Days: False
 sequential: False
  byteOrder: big
       xdef: [200.  202.5 205.  207.5 210.  212.5 215.  217.5 220.  222.5 225.  227.5
 230.  232.5 235.  237.5 240.  242.5 245.  247.5 250.  252.5 255.  257.5
 260.  262.5 265.  267.5 270.  272.5 275.  277.5 280.  282.5 285.  287.5
 290.  292.5 295.  297.5 300.  302.5 305.  307.5 310.  312.5 315.  317.5
 320.  322.5 325.  327.5 330. ]
       ydef: [75.  72.5 70.  67.5 65.  62.5 60.  57.5 55.  52.5 50.  47.5 45.  42.5
 40.  37.5 35.  32.5 30.  27.5 25.  22.5 20.  17.5 15. ]
       zdef: [1000.]
       tdef: ['2013-01-01T00:00:00' '2013-01-01T06:00:00' '2013-01-01T12:00:00'
 '2013-01-01T18:00:00']
       pdef: 
       vdef: [CtlVar: air      in shape (t=4, z=1, y=25, x=53)]
```

#### CTL with PDEF
For a CTL file `test10.ctl` containing `pdef` for describing e.g., [Lambert Conformal Conic](https://en.wikipedia.org/wiki/Lambert_conformal_conic_projection) projection dataset (usually output by numerical model like [WRF](https://www.mmm.ucar.edu/weather-research-and-forecasting-model)):
```text
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

```

[xgrads](https://github.com/miniufo/xgrads) will automatically recognize two sets of grid: one is native projected grid defining the data and the other is lat/lon grid that used for plotting.  The output of `ctl` and `ctl.pdef` is:
```python
>>> print(ctl)
   dsetPath: ./test10.dat
   descPath: ./test10.ctl
   indxPath: 
   stnmPath: 
      title: CUACE_emi_index
      undef: -9999.0
       zrev: False
       yrev: False
      dtype: 
   template: False
  periodicX: True
 cal365Days: False
 sequential: False
  byteOrder: little
       xdef: [ 62.33    62.3976  62.4652 ... 144.5992 144.6668 144.7344]
       ydef: [11.18   11.2476 11.3152 11.3828 11.4504 11.518  11.5856 11.6532 11.7208
 11.7884 11.856  11.9236 11.9912 12.0588 12.1264 12.194  12.2616 12.3292
 12.3968 12.4644 12.532  12.5996 12.6672 12.7348 12.8024 12.87   12.9376
 ......
 56.81   56.8776 56.9452 57.0128 57.0804 57.148  57.2156]
       zdef: [1.]
       tdef: ['2019-01-01T00:00:00']
       pdef: lcc
       vdef: [CtlVar: emi_index in shape (t=1, z=1, y=320, x=360)]


>>> print(ctl.pdef)
isize: 360
jsize: 320
proj: lcc
latref: 35.0
lonref: 103.5
iref: 180.5
jref: 160.5
Struelat: 30.0
Ntruelat: 60.0
slon: 103.5
dx: 15000.0
dy: 15000.0
```

Note that we can use function `get_data_projection()` in `utils.py` to get a `cartopy.crs` for plotting the data.  Detailed information can be found in [this notebook](https://github.com/miniufo/xgrads/blob/master/notebooks/Plot%20preprojected%20data%20defined%20in%20PDEF.ipynb).


#### CTL content as a string
Finally, if one write CTL file content as a `str` in python code, one can also parse it using the kwarg `content=CTL_content` as:
```python
CTL_content = \
    "dset ^binary.dat\n" \
    "* this is a comment line\n" \
    "title 10-deg resolution model\n" \
    "undef -9.99e8\n" \
    "xdef 36 linear   0 10\n" \
    "ydef 19 linear -90 10\n" \
    "zdef  1 linear   0  1\n" \
    "tdef  1 linear 00z01Jan2000 1dy\n" \
    "vars  1\n" \
    "test  1 99 test variable\n" \
    "endvars\n"

ctl = CtlDescriptor(content=CTL_content)
print(ctl)
```

with the output as:
```python
>>> print(ctl)
   dsetPath: ^binary.dat
   descPath: None
   indxPath: 
   stnmPath: 
      title: 10-deg
      undef: -999000000.0
       zrev: False
       yrev: False
      dtype: 
   template: False
  periodicX: True
 cal365Days: False
 sequential: False
  byteOrder: little
       xdef: [  0.  10.  20.  30.  40.  50.  60.  70.  80.  90. 100. 110. 120. 130.
 140. 150. 160. 170. 180. 190. 200. 210. 220. 230. 240. 250. 260. 270.
 280. 290. 300. 310. 320. 330. 340. 350.]
       ydef: [-90. -80. -70. -60. -50. -40. -30. -20. -10.   0.  10.  20.  30.  40.
  50.  60.  70.  80.  90.]
       zdef: [0.]
       tdef: ['2000-01-01T00:00:00']
       pdef: 
       vdef: [CtlVar: test     in shape (t=1, z=1, y=19, x=36)]
```
