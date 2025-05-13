import os
import pathlib
from astropy.table import Table
from .catalog import Catalog, CrossmatchResult

# The relative path to the actual catalog in csv format
_catalog_path = "data/catalog.csv"
# The absolute path to the actual catalog in csv format
catalog_filepath = os.path.join(
    pathlib.Path(__file__).resolve().parent,
    *_catalog_path.split('/')
)

# Construct a Catalog object from the csv file
catalog = Catalog.read(catalog_filepath, format="ascii.csv")

