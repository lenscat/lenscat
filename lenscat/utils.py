import os
import copy
import re
import decimal
import numpy as np
from matplotlib import pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
import requests
import ligo.skymap.plot

def convert_to_astropy_unit(table):
    """
    TODO: convert this to a proper docstring

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

def parse_skymap_str(skymap_str):
    # If skymap_str points to a valid filepath,
    # do nothing and return
    if os.path.exists(skymap_str):
        return skymap_str

    _supported_GWTC1_event_list = [
        "GW150914",
        "GW151012",
        "GW151226",
        "GW170104",
        "GW170608",
        "GW170729",
        "GW170809",
        "GW170814",
        "GW170817",
        "GW170818",
        "GW170823",
    ]
    if skymap_str in _supported_GWTC1_event_list:
        # skymap_str points to an event in GWTC-1
        dcc_url_template = "https://dcc.ligo.org/public/0157/P1800381/007/{}_skymap.fits.gz"
        url = dcc_url_template.format(skymap_str)
        if requests.options(url).ok:
            return url
    elif skymap_str.startswith('S'):
        # Maybe skymap_str is the name of a superevent
        gracedb_skymap_url_template = "https://gracedb.ligo.org/api/superevents/{}/files/{}.multiorder.fits"

        # Try using bilby skymap if possible
        url = gracedb_skymap_url_template.format(skymap_str, "Bilby")
        if requests.options(url).ok:
            return url
        
        # If bilby skymap is not available, try bayestar next
        url = gracedb_skymap_url_template.format(skymap_str, "bayestar")
        if requests.options(url).ok:
            return url

    # Exhausted all possible resolutions, give up
    raise ValueError(f"Does not recognize {skymap_str}")

def convert_from_ICRS_to_healpy_convention(RA, DEC):
    """
    TODO: convert this to a proper docstring

    Convert (RA, DEC) complying to the ICRS standard to
    the healpy convention where \theta runs from 0 at the north
    pole to \pi at the south pole, and $\phi$ runs from 0 to 2\pi

    Input: RA, DEC
    Output: theta, phi
    """
    return np.pi/2 - np.deg2rad(DEC), np.deg2rad(RA)

def plot_catalog(catalog, RA_unit="deg", filename="catalog.png", dark_theme=False):
    # Filter catalog by type
    galaxy_lenses = catalog.filter_by_type("galaxy")
    cluster_lenses = catalog.filter_by_type("cluster")

    fig = plt.figure(dpi=300)
    if RA_unit == "deg":
        _projection = "astro degrees mollweide"
    elif RA_unit == "hms":
        _projection = "astro hours mollweide"
    else:
        raise ValueError(f"Does not understand {RA_unit}")
    ax = plt.axes(projection=_projection)
    ax.grid()
    ax.scatter_coord(
        SkyCoord(ra=galaxy_lenses["RA"], dec=galaxy_lenses["DEC"]),
        color="deepskyblue",
        s=5,
        label="Galaxy",
        alpha = 0.95
    )
    ax.scatter_coord(
        SkyCoord(ra=cluster_lenses["RA"], dec=cluster_lenses["DEC"]),
        color="violet",
        s=5,
        label="Cluster",
        alpha = 0.95
    )
    ax.set_title("") # Override default title
    ax.legend(
        labelcolor="black",
        facecolor="whitesmoke",
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.2),
    )

    if dark_theme:
        ax.tick_params(colors="white")

    if filename is None:
        return fig
    else:
        plt.tight_layout()
        plt.savefig(filename, transparent=dark_theme)

def get_precision(x):
    return 10**(decimal.Decimal(str(x)).as_tuple().exponent)

def check_possible_duplicates(catalog):
    # Make a report listing all possible duplicates
    num_duplicates = 0

    # Loop over all entries in the catalog
    for i in range(len(catalog)):
        RA_center = catalog[i]["RA"]
        DEC_center = catalog[i]["DEC"]

        # Determine the search window width based on the precision of the coordinates
        RA_search_width = get_precision(RA_center)
        DEC_search_width = get_precision(DEC_center)
        
        # Perform a simple search on RA and DEC
        res = catalog.search(
            RA_range=(RA_center-RA_search_width/2, RA_center+RA_search_width/2),
            DEC_range=(DEC_center-DEC_search_width/2, DEC_center+DEC_search_width/2)
        )

        # If there are more than one entry in the search result,
        # then there are possible duplicates
        if len(res) > 1:
            num_duplicates += 1
            print(f"Possible duplicates found for {catalog[i]['name']} centered at ({RA_center}, {DEC_center}):")
            print(res)
            print("")

    print(f"Total number of possible duplicates: {num_duplicates}")
