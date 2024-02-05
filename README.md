# $\texttt{lenscat}$
![license](https://img.shields.io/github/license/lenscat/lenscat)
[![GitHub release](https://img.shields.io/github/v/release/lenscat/lenscat.svg)](https://github.com/lenscat/lenscat/releases)
[![Upload Python Package](https://github.com/lenscat/lenscat/actions/workflows/python-publish.yml/badge.svg)](https://github.com/lenscat/lenscat/actions/workflows/python-publish.yml)
[![Create Release and Tag](https://github.com/lenscat/lenscat/actions/workflows/create-release-tag.yml/badge.svg)](https://github.com/lenscat/lenscat/actions/workflows/create-release-tag.yml)

A public and community-contributed catalog of known strong gravitational lenses. 

![Known Lenses](https://raw.githubusercontent.com/lenscat/lenscat/main/catalog.png)

## Quickstart

The catalog is available as a _plain csv file_ under [lenscat/data/catalog.csv](https://github.com/lenscat/lenscat/blob/main/lenscat/data/catalog.csv). Alternatively, one can interact with the catalog using a [web app](https://lenscat.streamlit.app/) (mobile-friendly).

We also provide a _python package_ `lenscat`, available in pypi. Simply do
```bash
pip install lenscat
```
to install the latest version. Here we adopt the _continuous deployment_ paradigm (similar to [`astroquery`](https://github.com/astropy/astroquery)). Whenever there is a change in the catalog content, or a major change in the code, a new release will be available instantaneously on both [GitHub](https://github.com/lenscat/lenscat/releases) and [PyPI](https://pypi.org/project/lenscat/).

The code converts the catalog in the csv file into a custom `Catalog` object that is inherited from the `Table` object in `astropy`. To access the catalog, simply run
```python
In [1]: import lenscat; lenscat.catalog
Out[1]:
<Catalog length=4587>
     name         RA       DEC     zlens     type   grading
                 deg       deg
    str20      float64   float64   str15     str7     str9
------------- ---------- -------- -------- ------- ---------
   J0011-0845    2.83435  -8.7643        -  galaxy confirmed
   J0013+5119   3.348077  51.3183        -  galaxy confirmed
 PSJ0028+0631    7.09369   6.5317        -  galaxy confirmed
 PSJ0030-1525     7.5636 -15.4177 measured  galaxy confirmed
   J0047+2514 11.9465943  25.2411        -  galaxy confirmed
  HE0047-1756    12.6158 -17.6693    0.407  galaxy confirmed
          ...        ...      ...      ...     ...       ...
235614+023115  359.05923  2.52107    0.372  galaxy  probable
235730+010133 359.377643  1.02596    0.638  galaxy  probable
235811+003309   359.5475   0.5527    0.639 cluster  probable
235853+012406 359.721708 1.401816    0.481  galaxy  probable
235933+020823   359.8897   2.1398     0.43 cluster confirmed
235948-005913  359.95245 -0.98702    0.758  galaxy  probable
235952+004154   359.9698   0.6985    0.267 cluster  probable
```
Note that the code will try to assign the unit for each of the columns inferred from its name, and that it will hide the 'ref' column *by default*. One can show or hide the 'ref' column by calling `.show_ref()` and `.hide_ref()` on the `Catalog` object respectively.

Every `Catalog` object supports two features: basic searching with `.search()` and crossmatching with a skymap with `.crossmatch()`. Note that these function will return a `Catalog` object, and hence *they can be composed together* (e.g., `.crossmatch().search()`).

### Basic searching
This feature is implemented as `.search()`. One can search/filter by any combination of
- ranges of right ascension (specified as `RA_range=(RA_min, RA_max)`)
- ranges of declination (specified as `DEC_range=(DEC_min, DEC_max)`)
- ranges of lens redshift if available (specified as `zlens_range=(zlens_min, zlens_max)`)
- type of the lenses (specified as `lens_type`)
- grading of the lenses (specified as `grading`)

For example, to get a list of the cluster-scale lenses which are confirmed and with a redshift $z_{\mathrm{lens}} \geq 1$ together with the reference, run
```python
In [1]: import lenscat, numpy; lenscat.catalog.search(grading="confirmed", lens_type="cluster", zlens_range=(1,numpy.inf)).show_ref()
Out[1]:
<Catalog length=3>
     name         RA       DEC    zlens   type   grading                                        ref
                 deg       deg
    str20      float64   float64  str15   str7     str9                                        str171
------------- --------- --------- ----- ------- --------- ------------------------------------------------------------------
021118-042729 32.827087 -4.458069  1.02 cluster confirmed https://arxiv.org/abs/2004.00634 More et al. 2012 More et al. 2016
023100-062139   37.7516   -6.3608  1.17 cluster confirmed                                   https://arxiv.org/abs/2002.01611
220859+020655  332.2495    2.1153  1.04 cluster confirmed                                   https://arxiv.org/abs/2002.01611
```

### Crossmatching with a skymap
This feature is implemented as `.crossmatch()`. This function is simply a wrapper to the `crossmatch()` function in `ligo.skymap` which performs the cross-matching of a gravitational-wave (GW) skymap with a given list of coordinates. For example, to cross-match the GW skymap of GW190814 (download from [here](https://gracedb.ligo.org/apiweb/superevents/S190814bv/files/GW190814_PublicationSamples.multiorder.fits)) with only the confirmed lenses in the `lenscat` catalog, simply run
```python
In [1]: import lenscat; lenscat.catalog.search(grading="confirmed").crossmatch("GW190814_PublicationSamples.multiorder.fits")
Out[1]:
<Catalog length=761>
       name             RA         DEC       zlens     type   grading  searched probability   searched area
                       deg         deg                                                             deg2
      str20          float64     float64     str15     str7     str9         float64             float64
------------------ ----------- ------------ -------- ------- --------- -------------------- ------------------
  MACSJ0035.4-2015        8.85 -20.27083333    0.352 cluster confirmed   0.9881400284668479  46.21360633943737
          0047-281 12.42458333 -27.87380556    0.484  galaxy confirmed   0.9899384024067811  48.32495441567177
             A2813     10.8625 -20.61694444   0.2924 cluster confirmed   0.9924501903529119  51.82638259178721
     DESJ0145-3541    26.44493    -35.69093     0.49  galaxy confirmed   0.9998766459538251  82.04095381939021
     DESJ0138-2844    24.59567    -28.73555     0.44  galaxy confirmed   0.9999999999999968 151.91214978321142
     DESJ0130-3744    22.51201    -37.74938     0.68  galaxy confirmed   0.9999999999999968  173.7862404115262
AGEL 014253-183116 10.72041667 -18.52105556  0.63627  galaxy confirmed   0.9999999999999968 187.58212970468517
```

## Format

['name'] = Names of galaxies/galaxy clusters \
['RA [deg]'] = RA in dergees \
['DEC [deg]'] = DEC in degrees \
['zlens'] = Lens redshift (if known) \
['type'] = Type of lens (i.e. galaxy or galaxy cluster) \
['grading'] = Grading whether it is a confirmed lens or a probable lens (see individual references for internal grading systems) \
['ref'] = Reference to the corresponding catalog or study

## References

This catalog contains the known strong lenses from the following studies:

  - GLQ Database:
    https://research.ast.cam.ac.uk/lensedquasars/index.html

  - CLASH (Postman+2012):
    https://archive.stsci.edu/prepds/clash/

  - MUSES Cluster Followups (Richards+2020):
    https://cral-perso.univ-lyon1.fr/labo/perso/johan.richard/MUSE_data_release/

  - RELICS
    https://relics.stsci.edu/clusters.html

  - 37 Clusters from SDSS Giant Arcs Survey
    https://iopscience.iop.org/article/10.3847/1538-4365/ab5f13

  - An Extended Catalog of Galaxy–Galaxy Strong Gravitational Lenses Discovered in DES Using Convolutional Neural Networks
    https://iopscience.iop.org/article/10.3847/1538-4365/ab26b6#apjsab26b6t5

  - The AGEL Survey: Spectroscopic Confirmation of Strong Gravitational Lenses in the DES
    and DECaLS Fields Selected Using Convolutional Neural Networks
    https://arxiv.org/ftp/arxiv/papers/2205/2205.05307.pdf

  - LSD Survey
    https://web.physics.ucsb.edu/~tt/LSD/

  - (COSMOS) LensFlow: A Convolutional Neural Network in Search of Strong Gravitational Lenses
    https://ui.adsabs.harvard.edu/abs/2018ApJ...856...68P/abstract

  - SLACS. XIII. Galaxy-scale strong lens candidates
    https://ui.adsabs.harvard.edu/abs/2019yCat..18510048S/abstract

  - RINGFINDER: Automated Detection of Galaxy-scale Gravitational Lenses in Ground-based Multi-filter Imaging Data
    https://iopscience.iop.org/article/10.1088/0004-637X/785/2/144

  - Survey of Gravitationally-lensed Objects in HSC Imaging (SuGOHI) Candidate List
    https://www-utap.phys.s.u-tokyo.ac.jp/~oguri/sugohi/

## See also
[Master Lens Database](https://test.masterlens.org/index.php)

## Acknowledgements
This project was supported by the research grant no. VIL37766 and no. VIL53101 from Villum Fonden, and the DNRF Chair program grant no. DNRF162 by the Danish National Research Foundation.

This project has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 101131233.

We would also like to thank Jonah Kanner for introducing us the amazing [streamlit](https://streamlit.io) service that hosts the [web app](https://lenscat.streamlit.app/) for `lenscat`.
