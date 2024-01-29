import os
import pathlib
import copy
from astropy.table import Table
from ligo.skymap.io import read_sky_map
import astropy.units as u

from .utils import convert_to_astropy_unit, parse_skymap_str, crossmatch_with_catalog, search_catalog

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
    # Searched areas are in sqdeg
    table_with_result["searched area"].unit = u.deg**2

    # Sort by searched prob, then by seared area
    table_with_result.sort(["searched probability", "searched area"])

    return table_with_result

def search(RA_range=None, DEC_range=None, zlens_range=None, lens_type=None, grading=None):
    return search_catalog(
        catalog,
        RA_range=RA_range,
        DEC_range=DEC_range,
        zlens_range=zlens_range,
        lens_type=lens_type,
        grading=grading
    )
