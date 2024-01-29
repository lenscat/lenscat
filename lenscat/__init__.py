import os
import pathlib
from astropy.table import Table
from .catalog import Catalog

# The relative path to the actual catalog in csv format
_catalog_path = "data/catalog.csv"

# Construct a Catalog object from the csv file
catalog = Catalog.read(
    os.path.join(
        pathlib.Path(__file__).resolve().parent,
        *_catalog_path.split('/')
    ), format="ascii.csv")

