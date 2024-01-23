import os
import pathlib
import copy
from astropy.table import Table
from ligo.skymap.io import read_sky_map

from .utils import convert_to_astropy_unit, parse_skymap_str, crossmatch_with_catalog

# The relative path to the actual catalog in csv format
_catalog_path = "data/catalog.csv"

# Construct an astropy Table object from the csv file
catalog = Table.read(
    os.path.join(
        pathlib.Path(__file__).resolve().parent,
        *_catalog_path.split('/')
    ), format="ascii.csv")
convert_to_astropy_unit(catalog)

def crossmatch(skymap):
    if type(skymap) == Table:
        _skymap = skymap
    elif type(skymap) == str:
        # Preserve multiorder if possible
        _skymap = read_sky_map(parse_skymap_str(skymap), moc=True)
    else:
        raise ValueError(f"Does not support {skymap}")
    
    result = crossmatch_with_catalog(_skymap, catalog)
    table_with_result = copy.copy(catalog)
    table_with_result.add_column(
        result.searched_prob,
        name="searched probability"
    )
    table_with_result.add_column(
        result.searched_area,
        name="searched area"
    )

    return table_with_result