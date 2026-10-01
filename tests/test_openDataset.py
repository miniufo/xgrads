import numpy as np
import xarray as xr
from xgrads import open_CtlDataset, open_mfdataset


def test_template1():
    dset1 = open_CtlDataset('./ctls/test8.ctl')
    dset2 = open_CtlDataset('./ctls/test9.ctl')
    for level in range(len(dset1.time)):
        xr.testing.assert_allclose(dset1.air[level], dset2.air[level])


def test_open_mfdataset_template8():
    combined = open_mfdataset('./ctls/test8_*.ctl', parallel=False).load()
    expected = open_CtlDataset('./ctls/test8.ctl').load()
    xr.testing.assert_allclose(combined.air, expected.air)


def test_open_mfdataset_template9():
    combined = open_mfdataset('./ctls/test9_*.ctl', parallel=False).load()
    expected = open_CtlDataset('./ctls/test9.ctl').load()
    xr.testing.assert_allclose(combined.air, expected.air)


def test_template4():
    dset1 = open_CtlDataset('./ctls/test81.ctl')
    dset2 = open_CtlDataset('./ctls/test82.ctl')
    assert (dset1.x == dset2.x).all()
    assert (dset1.y == dset2.y).all()
    assert (dset1.air[0] == dset2.air).all()


def test_ensemble():
    dset1 = open_CtlDataset('./ctls/ecmf_medium_T2m1.ctl')
    expected = np.array([2.011963, 1.1813354, 1.1660767])
    assert np.isclose(dset1.t2[:, -1, -1, -1], expected).all()
