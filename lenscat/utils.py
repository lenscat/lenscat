import os
import copy
import re
import decimal
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import requests
import ligo.skymap.plot
from ligo.skymap.io import read_sky_map
from ligo.skymap.postprocess.util import find_greedy_credible_levels
import healpy as hp

def convert_to_astropy_unit(table: Table) -> None:
    """
    Convert columns in the given astropy table to the proper astropy unit if possible.

    Parameters
    ----------
    table : astropy.table.Table
        The table containing the columns to be converted.

    Notes
    -----
    - This function modifies the input table in-place.
    - Columns with names containing a pair of square brackets will be converted to the proper astropy unit.
    - After the conversion, the square brackets will be removed from the column names.

    Examples
    --------
    Suppose we have the following column names in a table:
    ["Name", "RA [deg]", "DEC [deg]"]

    After invoking this function, the table will have the following column names:
    ["Name", "RA", "DEC"]
    and the second and third columns will have a unit of deg.
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
    _supported_GWTC2_event_list = [
        "GW190408_181802",
        "GW190412",
        "GW190413_052954",
        "GW190413_134308",
        "GW190421_213856",
        "GW190424_180648",
        "GW190425",
        "GW190426_152155",
        "GW190503_185404",
        "GW190512_180714",
        "GW190513_205428",
        "GW190514_065416",
        "GW190517_055101",
        "GW190519_153544",
        "GW190521",
        "GW190521_074359",
        "GW190527_092055",
        "GW190602_175927",
        "GW190620_030421",
        "GW190630_185205",
        "GW190701_203306",
        "GW190706_222641",
        "GW190707_093326",
        "GW190708_232457",
        "GW190719_215514",
        "GW190720_000836",
        "GW190727_060333",
        "GW190728_064510",
        "GW190731_140936",
        "GW190803_022701",
        "GW190814",
        "GW190828_063405",
        "GW190828_065509",
        "GW190909_114149",
        "GW190910_112807",
        "GW190915_235702",
        "GW190924_021846",
        "GW190929_012149",
        "GW190930_133541",
    ]
    _supported_GWTC3_event_list = [
        "GW191103_012549",
        "GW191105_143521",
        "GW191109_010717",
        "GW191113_071753",
        "GW191126_115259",
        "GW191127_050227",
        "GW191129_134029",
        "GW191204_110529",
        "GW191204_171526",
        "GW191215_223052",
        "GW191216_213338",
        "GW191219_163120",
        "GW191222_033537",
        "GW191230_180458",
        "GW200105_162426",
        "GW200112_155838",
        "GW200115_042309",
        "GW200128_022011",
        "GW200129_065458",
        "GW200202_154313",
        "GW200208_130117",
        "GW200208_222617",
        "GW200209_085452",
        "GW200210_092254",
        "GW200216_220804",
        "GW200219_094415",
        "GW200220_061928",
        "GW200220_124850",
        "GW200224_222234",
        "GW200225_060421",
        "GW200302_015811",
        "GW200306_093714",
        "GW200308_173609",
        "GW200311_115853",
        "GW200316_215756",
        "GW200322_091133",
    ]
    if skymap_str in _supported_GWTC1_event_list:
        # skymap_str points to an event in GWTC-1
        dcc_url_template = "https://dcc.ligo.org/public/0157/P1800381/007/{}_skymap.fits.gz"
        url = dcc_url_template.format(skymap_str)
        if requests.options(url).ok:
            return url
    elif skymap_str in _supported_GWTC2_event_list:
        # Need to download *all* GWTC-2 skymaps first
        # NOTE This is not going to be very fast
        if not os.path.exists(".cache/all_skymaps.tar"):
            os.makedirs(".cache", exist_ok=True)
            response = requests.get("https://dcc.ligo.org/public/0169/P2000223/007/all_skymaps.tar")
            with open(".cache/all_skymaps.tar", 'wb') as f:
                f.write(response.content)
        # Extract from the tarball
        import tarfile
        _filename_template = "{}_PublicationSamples.fits"
        filename = _filename_template.format(skymap_str)
        tarfile.open(".cache/all_skymaps.tar").extract("all_skymaps/{}".format(filename), path=".cache")
        return f".cache/all_skymaps/{filename}"
    elif skymap_str in _supported_GWTC3_event_list:
        # Need to download *all* GWTC-3 skymaps first
        # NOTE This is not going to be very fast
        if not os.path.exists(".cache/IGWN-GWTC3p0-v2-PESkyLocalizations.tar.gz"):
            os.makedirs(".cache", exist_ok=True)
            response = requests.get(
                "https://zenodo.org/api/records/8177023/files/IGWN-GWTC3p0-v2-PESkyLocalizations.tar.gz/content"
            )
            with open(".cache/IGWN-GWTC3p0-v2-PESkyLocalizations.tar.gz", 'wb') as f:
                f.write(response.content)
        # Extract from the tarball
        import tarfile
        _filename_template = "IGWN-GWTC3p0-v2-{}_PEDataRelease_cosmo_reweight_C01:Mixed.fits"
        filename = _filename_template.format(skymap_str)
        tarfile.open(".cache/IGWN-GWTC3p0-v2-PESkyLocalizations.tar.gz").extract(
            "IGWN-GWTC3p0-v2-PESkyLocalizations/{}".format(filename), path=".cache"
        )
        return f".cache/IGWN-GWTC3p0-v2-PESkyLocalizations/{filename}"
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
    elif requests.options(skymap_str).ok:
        return skymap_str
    elif skymap_str.startswith("http"):
        return skymap_str # Try anyway

    # Exhausted all possible resolutions, give up
    raise ValueError(f"Does not recognize {skymap_str}")

def parse_skymap(skymap, moc=False):
    if type(skymap) == str:
        # Preserve multiorder if possible
        _skymap = read_sky_map(parse_skymap_str(skymap), moc=moc)
    else:
        raise ValueError(f"Does not understand {skymap}")

    return _skymap

def convert_from_ICRS_to_healpy_convention(RA, DEC):
    """
    Convert (RA, DEC) complying to the ICRS standard to
    the healpy convention where theta runs from 0 at the north
    pole to pi at the south pole, and phi runs from 0 to 2*pi

    Parameters
    ----------
    RA : float
        Right Ascension in degrees.
    DEC : float
        Declination in degrees.

    Returns
    -------
    theta : float
        Polar angle in radians, ranging from 0 at the north pole to pi at the south pole.
    phi : float
        Azimuthal angle in radians, ranging from 0 to 2*pi.
    """
    return np.pi/2 - np.deg2rad(DEC), np.deg2rad(RA)

def get_ra_dec_from_skymap(skymap):
    """
    Get the RA and DEC of the pixel with the maximum probability in a skymap

    Parameters
    ----------
    skymap : array_like
        The skymap containing the probability values.

    Returns
    -------
    RA : float
        The right ascension (RA) of the maximum probability pixel in degrees.
    DEC : float
        The declination (DEC) of the maximum probability pixel in degrees.
    """
    index_of_max = np.argmax(skymap)
    nside = hp.npix2nside(len(skymap))
    theta, phi = hp.pix2ang(nside, index_of_max)
    return np.rad2deg(phi), np.rad2deg(np.pi/2-theta)

def make_transparent_colormap(colormap: str) -> ListedColormap:
    """
    Produce a transparent colormap from the given colormap.

    Parameters
    ----------
    colormap : str
        The name of the colormap to be made transparent.

    Returns
    -------
    ListedColormap
        The transparent colormap.
    """
    cmap = plt.get_cmap(colormap)
    cmap_transparent = cmap(np.arange(cmap.N))
    alphas = np.linspace(0, 1, cmap.N)
    bkgrd = np.asarray([1., 1., 1.,])
    for j in range(cmap.N):
        cmap_transparent[j,:-1] = cmap_transparent[j,:-1]*alphas[j] + bkgrd*(1. - alphas[j])
    # Mess with the alpha levels for the first 20% of the colors
    _idxs = np.arange(0, int(0.2*cmap.N), 1)
    cmap_transparent[_idxs,-1] = np.linspace(0, 1, len(_idxs))
    cmap_transparent = ListedColormap(cmap_transparent)

    return cmap_transparent

def plot_catalog(catalog, RA_unit="deg", filename="catalog.png", dark_theme=False, **axes_kwargs):
    # Filter catalog by type
    galaxy_lenses = catalog.filter_by_type("galaxy")
    cluster_lenses = catalog.filter_by_type("cluster")
    group_lenses = catalog.filter_by_type("group")

    if RA_unit == "deg":
        _projection = "astro degrees mollweide"
    elif RA_unit == "hms":
        _projection = "astro hours mollweide"
    else:
        raise ValueError(f"Does not understand {RA_unit}")
    _kwargs = {
        "projection": _projection
    }
    _kwargs.update(axes_kwargs)

    fig = plt.figure(dpi=300)
    ax = plt.axes(**_kwargs)
    ax.grid()
    if len(galaxy_lenses) > 0:
        ax.scatter_coord(
            SkyCoord(ra=galaxy_lenses["RA"], dec=galaxy_lenses["DEC"]),
            color="deepskyblue",
            s=5,
            label="Galaxy",
            alpha = 0.95
        )
    if len(group_lenses) > 0:
        ax.scatter_coord(
            SkyCoord(ra=group_lenses["RA"], dec=group_lenses["DEC"]),
            color="orange",
            s=5,
            label="Group",
            alpha = 0.95
        )
    if len(cluster_lenses) > 0:
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
        loc="upper right",
        bbox_to_anchor=(1.3, 1.0),
    )

    if dark_theme:
        ax.tick_params(colors="white")

    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, transparent=dark_theme, bbox_inches="tight")

    return fig

def plot_crossmatch(
        skymap,
        crossmatch_res,
        searched_prob_threshold=1.0,
        RA_unit="deg",
        max_zoom_radius=15,
        filename="crossmatch_result.png",
        dark_theme=False
):
    # Check if there is "searched probability" in crossmatch_res
    if "searched probability" not in crossmatch_res.colnames:
        raise ValueError("No searched probability found in table")
    if "searched area" not in crossmatch_res.colnames:
        raise ValueError("No searched area found in table")

    assert searched_prob_threshold >= 0 and searched_prob_threshold <= 1, \
        "searched_prob_threshold must be between 0 and 1"

    # Filter res by searched_prob_threshold
    res = crossmatch_res[crossmatch_res["searched probability"] <= searched_prob_threshold]  
    if len(res) == 0:
        raise ValueError("No entries satisfy the searched probability threshold")

    # Read skymap
    _skymap = parse_skymap(skymap, moc=False)

    # Estimate the zoom radius needed
    # Maybe this is a sign that we should use mollewide projection instead
    zoom_radius = 2*np.sqrt(res[-1]["searched area"]) # A fudge factor of 2
    # Center the plot at the maxP estimate for ra and dec
    ra_maxP, dec_maxP = get_ra_dec_from_skymap(_skymap[0])

    # If the zoom radius is above max_zoom_radius, abort
    if zoom_radius > max_zoom_radius:
        # Use mollweide projection
        fig = plot_catalog(
            res,
            RA_unit=RA_unit,
            filename=None,
            dark_theme=dark_theme,
        )
        ax = fig.gca()
        ax.imshow_hpx(_skymap[0], cmap=make_transparent_colormap('cylon'))
    else:
        axes_kwargs = {
            "center": '{}d {}d'.format(int(ra_maxP), int(dec_maxP)),
            "radius": "{} deg".format(int(zoom_radius)),
        }

        if RA_unit == "deg":
            axes_kwargs["projection"] = "astro degrees zoom"
        elif RA_unit == "hms":
            axes_kwargs["projection"] = "astro hours zoom"

        fig = plot_catalog(
            res,
            filename=None,
            dark_theme=dark_theme,
            **axes_kwargs
        )

        ax = fig.gca()
        label_color = 'w' if dark_theme else 'k'

        # Plot skymap
        ax.imshow_hpx(_skymap[0], cmap=make_transparent_colormap('cylon'))
        # Plot credible level contours
        contour = ax.contour_hpx(
            (find_greedy_credible_levels(_skymap[0]), "ICRS"),
            nested=False,
            colors='grey',
            linewidths=0.5,
            zorder=2,
            alpha=0.6,
            levels=0.15*np.array([1, 2, 3, 4, 5, 6]),
        )
        ax.clabel(contour, inline=True, fontsize=6)

        res.sort("DEC") # Sort by DEC
        align_left_horizontally = 1
        for row in res:
            _alignment = "left" if align_left_horizontally == 1 else "right"
            ax.text_coord(
                SkyCoord(ra=row["RA"]*u.deg, dec=row["DEC"]*u.deg),
                row["name"],
                fontsize="small",
                horizontalalignment=_alignment,
                color=label_color,
            )
            align_left_horizontally *= -1

        ax.set_xlabel("right ascension", color=label_color)
        ax.set_ylabel("declination", color=label_color)

    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, transparent=dark_theme, bbox_inches="tight")
    
    return fig

def get_precision(x):
    """
    Get the precision of a number.

    Parameters
    ----------
    x : float/int
        The number to be checked.
    
    Returns
    -------
    float
        The precision of the number.

    Examples
    --------
    >>> get_precision(0.0001)
    0.0001
    >>> get_precision(1.001)
    0.001
    >>> get_precision(10)
    1
    >>> get_precision(1.0)
    0.1
    """
    return 10**(decimal.Decimal(str(x)).as_tuple().exponent)

def check_possible_duplicates(catalog):
    """
    Check for possible duplicates in the given catalog.

    Parameters
    ----------
    catalog : lenscat.catalog.Catalog
        The catalog to be checked for possible duplicates.
    
    Notes
    -----
    - This function prints out a report listing all possible duplicates.
    - A possible duplicate is defined as an entry in the catalog that has the same
        RA and DEC as another entry in the catalog within measurement error/precision.
    """
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
