import os
import pathlib
from astropy.table import Table
from astropy.coordinates import SkyCoord
from ligo.skymap.io import read_sky_map
from ligo.skymap.postprocess import crossmatch as ligoskymap_crossmatch

from .utils import convert_to_astropy_unit, parse_skymap_str

# The relative path to the actual catalog in csv format
_catalog_path = "data/catalog.csv"

# Construct an astropy Table object from the csv file
catalog = Table.read(
    os.path.join(
        pathlib.Path(__file__).resolve().parent,
        *_catalog_path.split('/')
    ), format="ascii.csv")
convert_to_astropy_unit(catalog)
