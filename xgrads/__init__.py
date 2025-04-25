# -*- coding: utf-8 -*-
from .core import CtlDescriptor
from .io import open_CtlDataset, open_mfdataset
from .utils import interp_to_latlon, get_coordinates_from_PDEF, \
                   get_data_projection, oacressman

__version__ = "0.2.7"
