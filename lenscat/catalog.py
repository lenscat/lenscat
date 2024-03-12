import copy
import re
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from ligo.skymap.postprocess import crossmatch as ligoskymap_crossmatch
from .utils import convert_to_astropy_unit, parse_skymap
from .utils import plot_catalog, plot_crossmatch

class _Catalog(Table):
    _allowed_type = ["galaxy", "cluster"]
    _allowed_grading = ["confirmed", "probable"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        convert_to_astropy_unit(self)
        self.hide_ref() # By default, hide the ref column

    def hide_ref(self):
        try:
            self.pprint_exclude_names.add("ref")
        except:
            pass # Fail silently
        return self

    def show_ref(self):
        try:
            self.pprint_exclude_names.remove("ref")
        except:
            pass # Fail silently
        return self

    def filter_by_RA(self, RA_min=0.0, RA_max=360.0):
        assert RA_min < RA_max, "RA_max must be greater than RA_min"
        return self[(self["RA"] >= RA_min) & (self["RA"] <= RA_max)]

    def filter_by_DEC(self, DEC_min=-90, DEC_max=90):
        assert DEC_min < DEC_max, "DEC_max must be greater than DEC_min"
        return self[(self["DEC"] >= DEC_min) & (self["DEC"] <= DEC_max)]

    def filter_by_zlens(self, zlens_min=0.0, zlens_max=np.inf):
        assert zlens_min < zlens_max, "zlens_max must be greater than zlens_min"
        # zlens is unfortunately heterogeneous
        # Loop over each row in the catalog to see if zlens is something
        # that looks like a number
        _regex = r"^(0|[1-9]\d*)(\.\d+)?(.*)$"
        regularized_zlens_arr = []
        for zlens_str in self["zlens"]:
            m = re.search(_regex, str(zlens_str))
            if m is None:
                regularized_zlens_arr.append(np.nan)
            else:
                if m.group(2) is not None:
                    regularized_zlens_arr.append(float(m.group(1)+m.group(2)))
                else:
                    regularized_zlens_arr.append(float(m.group(1)))
        regularized_zlens_arr = np.array(regularized_zlens_arr)

        return self[(regularized_zlens_arr >= zlens_min) & (regularized_zlens_arr <= zlens_max)]

    def filter_by_type(self, lens_type):
        assert lens_type in self._allowed_type, f"{lens_type} is not one of the recognized types: "+", ".join(self._allowed_type)
        return self[self["type"] == lens_type]

    def filter_by_grading(self, grading):
        assert grading in self._allowed_grading, f"{grading} is not one of the recognized gradings: "+", ".join(self._allowed_grading)
        return self[self["grading"] == grading]

    def search(self, RA_range=None, DEC_range=None, zlens_range=None, lens_type=None, grading=None):
        filtered_catalog = self # Create a 'view' of the catalog table
        if RA_range is not None:
            filtered_catalog = filtered_catalog.filter_by_RA(*RA_range)
        if DEC_range is not None:
            filtered_catalog = filtered_catalog.filter_by_DEC(*DEC_range)
        if lens_type is not None:
            filtered_catalog = filtered_catalog.filter_by_type(lens_type)
        if grading is not None:
            filtered_catalog = filtered_catalog.filter_by_grading(grading)
    
        # Always filter by zlens last
        if zlens_range is not None:
            filtered_catalog = filtered_catalog.filter_by_zlens(*zlens_range)

        return filtered_catalog

    def plot(self, **kwargs):
        return plot_catalog(self, **kwargs)

class CrossmatchResult(_Catalog):
    def __init__(self, *args, **kwargs):
        self.skymap = kwargs.pop("skymap", None)
        super().__init__(*args, **kwargs)

    def plot(self, **kwargs):
        return plot_crossmatch(
            self.skymap,
            self,
            **kwargs
        )

class Catalog(_Catalog):
    def crossmatch(self, skymap):
        _skymap = parse_skymap(skymap, moc=True)
    
        coordinates = SkyCoord(self["RA"], self["DEC"])
        result = ligoskymap_crossmatch(_skymap, coordinates)
        crossmatch_result = CrossmatchResult(self.columns[:], skymap=skymap)
        crossmatch_result.add_column(
            result.searched_prob,
            name="searched probability"
        )
        crossmatch_result.add_column(
            result.searched_area,
            name="searched area"
        )
        # Searched areas are in sqdeg
        crossmatch_result["searched area"].unit = u.deg**2

        # Sort by searched prob, then by seared area
        crossmatch_result.sort(["searched probability", "searched area"])

        # Regularize searched probability
        for idx in range(len(crossmatch_result)):
            if crossmatch_result[idx]["searched probability"] > 1.0:
                crossmatch_result[idx]["searched probability"] = 1.0

        return crossmatch_result
