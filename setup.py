from pathlib import Path
import re
from setuptools import find_packages, setup


ROOT = Path(__file__).parent
version = re.search(r'__version__ = "(.*?)"', (ROOT / 'xgrads/__init__.py').read_text(encoding='utf-8')).group(1)

setup(
    name='xgrads',
    version=version,
    description='Parse and read ctl files commonly used by GrADS.',
    long_description=(ROOT / 'README.md').read_text(encoding='utf-8'),
    long_description_content_type='text/markdown',
    url='https://github.com/miniufo/xgrads',
    author='miniufo',
    author_email='miniufo@163.com',
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],
    keywords='grads opengrads xarray dask',
    packages=find_packages(exclude=['docs', 'tests', 'ctls', 'notebooks', 'pics', 'private']),
    install_requires=['numpy', 'xarray', 'dask', 'pyproj', 'numba'],
    extras_require={'test': ['pytest', 'pytest-cov', 'h5netcdf']},
)
