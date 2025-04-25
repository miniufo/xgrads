# xgrads

[![DOI](https://zenodo.org/badge/244529165.svg)](https://zenodo.org/badge/latestdoi/244529165)
![GitHub](https://img.shields.io/github/license/miniufo/xgrads)
[![Documentation Status](https://readthedocs.org/projects/xgrads/badge/?version=latest)](https://xgrads.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/xgrads.svg)](https://badge.fury.io/py/xgrads)
[![Publish to PyPI](https://github.com/miniufo/xgrads/actions/workflows/PyPI-publish.yml/badge.svg)](https://github.com/miniufo/xgrads/actions/workflows/PyPI-publish.yml)
[![pytest](https://github.com/miniufo/xgrads/actions/workflows/tests.yml/badge.svg)](https://github.com/miniufo/xgrads/actions/workflows/tests.yml)
[![Build Status](https://app.travis-ci.com/miniufo/xgrads.svg?branch=master)](https://app.travis-ci.com/miniufo/xgrads)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f86676a904ef4ebfa4c2f9cda46264f5)](https://app.codacy.com/gh/miniufo/xgrads/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

![3D plot](https://raw.githubusercontent.com/miniufo/xgrads/master/pics/3D.png)


## 1. Introduction
The Grid Analysis and Display System ([GrADS](http://cola.gmu.edu/grads/) or [OpenGrADS](http://www.opengrads.org/)) is a widely used software for easy access, manipulation, and visualization of earth science data.  It uses a [descriptor (or control) file with a suffix `.ctl`](http://cola.gmu.edu/grads/gadoc/descriptorfile.html) to  describe a raw binary 4D dataset.  The `ctl` file is similar to the header information of a [NetCDF](https://www.unidata.ucar.edu/software/netcdf/docs/file_structure_and_performance.html) file, containing all the information about dimensions, attributes, and variables except for the variable data.

This python package [`xgrads`](https://github.com/miniufo/xgrads) is designed for parse and read the `.ctl` file commonly used by [GrADS](http://cola.gmu.edu/grads/).  Right now it can parse various kinds of `.ctl` files.  However, only the commonly used raw binary 4D datasets can be read using [`dask`](https://dask.org/) and return as a [`xarray.Dataset`](http://xarray.pydata.org/en/stable/)  Other types of binary data, like `dtype` is `station` or`grib`, may be supported in the future.

---
## 2. How to install
**Requirements**
`xgrads` is developed under the environment with `xarray` (=version 0.15.0), `dask` (=version 2.11.0), `numpy` (=version 1.15.4), `cartopy` (=version 0.17.0), and `pyproj` (=version 1.9.6).  Older versions of these packages are not well tested.

**Install via pip**
```
pip install xgrads
```

**Install via conda**
```
conda install -c conda-forge xgrads
```

**Install from github**
```
git clone https://github.com/miniufo/xgrads.git
cd xgrads
python setup.py install
```


---
## 3. Examples
### 3.1 Parse a `.ctl` file
Parsing a `.ctl` file is pretty simple using the following code:
```python
from xgrads import CtlDescriptor

ctl = CtlDescriptor(file='test.ctl')

# print all the info in ctl file
print(ctl)
```

If you have already load the ASCII content in the `.ctl` file, you can do it as:
```python
content = \
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

ctl = CtlDescriptor(content=content)

# print all the info
print(ctl)
```
---

### 3.2 Read binary data into a `xarray.Dataset`
Reading a `.ctl` related binary data file is also pretty simple using the following code:
```python
from xgrads import open_CtlDataset

dset = open_CtlDataset('test.ctl')

# print all the info in ctl file
print(dset)
```

Then you have the `dset` as a `xarray.Dataset`.  This is similar to [`xarray.open_dataset`](http://xarray.pydata.org/en/stable/generated/xarray.open_dataset.html) that use [`dask`](https://dask.org/) to chunk (buffer) parts of the whole dataset in physical memory if the whole dataset is too large to fit in.

If there are many `.ctl` files in a folder, we can also open all of them in a single call of `open_mfdataset` as:
```python
from xgrads import open_mfDataset

dset = open_mfDataset('./folder/*.ctl')

# print all the info in ctl file
print(dset)
```
assuming that every `.ctl` file has similar data structure except the time step is different.  This is similar to [`xarray.open_mfdataset`](http://xarray.pydata.org/en/v0.12.3/generated/xarray.open_mfdataset.html).

---

### 3.3 Convert a GrADS dataset to a NetCDF dataset
With the above functionality, it is easy to convert a `.ctl` ([GrADS](http://cola.gmu.edu/grads/)) dataset to a `.nc` ([NetCDF](https://www.unidata.ucar.edu/software/netcdf/docs/file_structure_and_performance.html)) dataset:
```python
from xgrads import open_CtlDataset

open_CtlDataset('input.ctl').to_netcdf('output.nc')
```
