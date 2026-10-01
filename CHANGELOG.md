# Changelog

## Unreleased

- Modernized Python packaging with `pyproject.toml`.
- Expanded CI coverage to Python 3.9–3.13 on Ubuntu and Windows.
- Added documentation build verification.
- Added multi-file dataset tests and removed the dependency on the remote xarray tutorial dataset.
- Updated the conda recipe and documentation dependencies.
- Fixed the documentation build, which was failing on a notebook cell that
  pandoc parsed as a YAML metadata block, and enabled strict warning checks.
- Declared the optional dependencies the test suite needs on a clean install
  (`cartopy`, `h5netcdf`, `h5py`).
