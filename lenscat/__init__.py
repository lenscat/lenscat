import os
import pathlib
from astropy.table import Table

_catalog_path = "data/catalog.csv"

catalog = Table.read(
    os.path.join(
        pathlib.Path(__file__).resolve().parent,
        *_catalog_path.split('/')
    ), format="ascii.csv")