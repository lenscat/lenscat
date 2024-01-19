import os
import pathlib
from astropy.table import Table
import copy
import re
import astropy.units as u

def convert_to_astropy_unit(table):
    """
    Loop over the column names of the given table,
    if a column name contains a pair of square bracket,
    convert that column with the proper astropy unit if possible. 
    After the conversion, the square bracket will be removed from the name.
    Note that changes are made **in-place**.

    For example, suppose we have the following column names in a table:
    ["Name", "RA [deg]", "DEC [deg]"]

    After invoking this function, the table will have the following column names
    ["Name", "RA", "DEC"]
    and the second and third columns will have a unit of deg
    """
    _regex = r"\[([A-Za-z0-9_]+)\]"
    column_names = copy.deepcopy(table.columns)
    for colname in column_names:
        m = re.search(_regex, colname)

        if m is not None:
            unit_name = m.group(1)
            # Search for this in astropy.units
            if unit_name in dir(u):
                unit = eval(f"u.{unit_name}")
                # Assign unit
                table[colname].unit = unit
                # Rename column
                new_name = re.sub(_regex, '', colname).strip()
                table.rename_column(colname, new_name)

# The relative path to the actual catalog in csv format
_catalog_path = "data/catalog.csv"

# Construct an astropy Table object from the csv file
catalog = Table.read(
    os.path.join(
        pathlib.Path(__file__).resolve().parent,
        *_catalog_path.split('/')
    ), format="ascii.csv")
convert_to_astropy_unit(catalog)
