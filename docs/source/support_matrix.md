# Supported GrADS features

The table below describes the intended support scope. Features marked
experimental should be tested with representative files before production use.

| Feature | Status | Notes |
| --- | --- | --- |
| Regular binary 4D datasets | Supported | Read through xarray and dask |
| `TDEF`, `XDEF`, `YDEF`, `ZDEF` | Supported | Linear and common level definitions |
| Template datasets | Supported | Use `open_mfdataset` for multiple files |
| Ensemble definitions (`EDEF`) | Supported | Covered by test fixtures |
| Byte order and sequential records | Supported | Depends on descriptor options |
| `PDEF` projected grids | Supported | Supported projections are documented by the API |
| Station data | Experimental | Validate against the target station format |
| GRIB data | Not supported | GrADS parsing is not a GRIB decoder |
| NetCDF conversion | Supported | Uses xarray's `to_netcdf` |

If a descriptor uses an unsupported option, please open an issue with a small
descriptor and a reproducible sample or a synthetic replacement for the data.
