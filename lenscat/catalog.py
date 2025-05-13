import copy
import re
import numpy as np
import astropy.units as u
from astropy.table import Table, TableAttribute
from astropy.coordinates import SkyCoord
from ligo.skymap.postprocess import crossmatch as ligoskymap_crossmatch
from .utils import convert_to_astropy_unit, parse_skymap
from .utils import plot_catalog, plot_crossmatch

class _Catalog(Table):
    _allowed_type = ["galaxy", "cluster", "group"]
    _allowed_grading = ["confident", "probable"]

    def __init__(self, *args, **kwargs):
        kwargs.update({
            "copy": False,
        })
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

    def filter_by_RA(self, RA_min=0.0, RA_max=360.0, return_copy=True):
        assert RA_min < RA_max, "RA_max must be greater than RA_min"
        assert RA_min >= 0.0 and RA_max <= 360.0, "RA must be between 0 and 360"
        # Use for loop instead of slicing to minimize memory usage
        filter = np.zeros(len(self), dtype=bool)
        for idx in range(len(self)):
            if self[idx]["RA"] >= RA_min and self[idx]["RA"] <= RA_max:
                filter[idx] = True
        if return_copy:
            return self[filter]
        else:
            return filter

    def filter_by_DEC(self, DEC_min=-90, DEC_max=90, return_copy=True):
        assert DEC_min < DEC_max, "DEC_max must be greater than DEC_min"
        assert DEC_min >= -90.0 and DEC_max <= 90.0, "DEC must be between -90 and 90"
        # Use for loop instead of slicing to minimize memory usage
        filter = np.zeros(len(self), dtype=bool)
        for idx in range(len(self)):
            if self[idx]["DEC"] >= DEC_min and self[idx]["DEC"] <= DEC_max:
                filter[idx] = True
        if return_copy:
            return self[filter]
        else:
            return filter

    def filter_by_zlens(self, zlens_min=0.0, zlens_max=np.inf, return_copy=True):
        assert zlens_min < zlens_max, "zlens_max must be greater than zlens_min"
        assert zlens_min >= 0.0, "zlens_min must be greater than 0"
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

        # Use for loop instead of slicing to minimize memory usage
        filter = np.zeros(len(self), dtype=bool)
        for idx in range(len(self)):
            if regularized_zlens_arr[idx] >= zlens_min and regularized_zlens_arr[idx] <= zlens_max:
                filter[idx] = True
        if return_copy:
            return self[filter]
        else:
            return filter

    def filter_by_type(self, lens_type, return_copy=True):
        assert lens_type in self._allowed_type, f"{lens_type} is not one of the recognized types: "+", ".join(self._allowed_type)
        filter = (self["type"] == lens_type)
        if return_copy:
            return self[filter]
        else:
            return filter

    def filter_by_grading(self, grading, return_copy=True):
        assert grading in self._allowed_grading, f"{grading} is not one of the recognized gradings: "+", ".join(self._allowed_grading)
        filter = (self["grading"] == grading)
        if return_copy:
            return self[filter]
        else:
            return filter

    def search(self, RA_range=None, DEC_range=None, zlens_range=None, lens_type=None, grading=None):
        filter = np.ones(len(self), dtype=bool)
        if RA_range is not None:
            filter *= self.filter_by_RA(*RA_range, return_copy=False)
        if DEC_range is not None:
            filter *= self.filter_by_DEC(*DEC_range, return_copy=False)
        if lens_type is not None:
            filter *= self.filter_by_type(lens_type, return_copy=False)
        if grading is not None:
            filter *= self.filter_by_grading(grading, return_copy=False)
    
        # Always filter by zlens last
        if zlens_range is not None:
            filter *= self.filter_by_zlens(*zlens_range, return_copy=False)

        # Create a new copy at the very end
        return self[filter]

    def plot(self, **kwargs):
        return plot_catalog(self, **kwargs)

class CrossmatchResult(_Catalog):
    skymap = TableAttribute(default=None)

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
