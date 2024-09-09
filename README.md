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
import lenscat; lenscat.catalog
```
and this will show a formatted table to the output. For example,
```output_catalog
<Catalog length=4587>
     name         RA      DEC    zlens   type   grading 
                 deg      deg                           
    str20      float64  float64  str15   str7     str9  
------------- --------- -------- ----- ------- ---------
   J0011-0845   2.83435  -8.7643     -  galaxy confident
   J0013+5119  3.348077  51.3183     -  galaxy confident
          ...       ...      ...   ...     ...       ...
235948-005913 359.95245 -0.98702 0.758  galaxy  probable
235952+004154  359.9698   0.6985 0.267 cluster  probable
 
```
Note that the code will try to assign the unit for each of the columns inferred from its name, and that it will hide the 'ref' column *by default*. One can show or hide the 'ref' column by calling `.show_ref()` and `.hide_ref()` on the `Catalog` object respectively.

Every `Catalog` object supports three features: basic searching with `.search()`, crossmatching with a skymap with `.crossmatch()`, and visualizing with `.plot()`. Note that these function will return a `Catalog` object, and hence *they can be composed together* (e.g., `.crossmatch().search()`).

### Basic searching
This feature is implemented as `.search()`. One can search/filter by any combination of
- ranges of right ascension (specified as `RA_range=(RA_min, RA_max)`)
- ranges of declination (specified as `DEC_range=(DEC_min, DEC_max)`)
- ranges of lens redshift if available (specified as `zlens_range=(zlens_min, zlens_max)`)
- type of the lenses (specified as `lens_type`)
- grading of the lenses (specified as `grading`)

For example, to get a list of the cluster-scale lenses which are confidently identified and with a redshift $z_{\mathrm{lens}} \geq 1$ together with the reference, run
```python
import lenscat, numpy
lenscat.catalog.search(grading="confident", lens_type="cluster", zlens_range=(1,numpy.inf)).show_ref()
```
The output would be something like
```output_search_confident_cluster
<Catalog length=3>
     name         RA       DEC    zlens   type   grading                                                 ref                                               
                 deg       deg                                                                                                                             
    str20      float64   float64  str15   str7     str9                                                 str227                                             
------------- --------- --------- ----- ------- --------- -------------------------------------------------------------------------------------------------
021118-042729 32.827087 -4.458069  1.02 cluster confident https://arxiv.org/abs/2004.00634 https://arxiv.org/abs/1109.1821 https://arxiv.org/abs/1504.05587
023100-062139   37.7516   -6.3608  1.17 cluster confident                                                                  https://arxiv.org/abs/2002.01611
220859+020655  332.2495    2.1153  1.04 cluster confident                                                                  https://arxiv.org/abs/2002.01611
 
```

### Crossmatching with a skymap
This feature is implemented as `.crossmatch()`. This function is simply a wrapper to the `crossmatch()` function in `ligo.skymap` which performs the cross-matching of a gravitational-wave (GW) skymap with a given list of coordinates. For example, to cross-match the GW skymap of GW170814 (download from [here](https://dcc.ligo.org/public/0157/P1800381/007/GW170814_skymap.fits.gz)) with only galaxy-scale lenses in the `lenscat` catalog, simply run
```python
import lenscat
lenscat.catalog.search(lens_type="galaxy").crossmatch("GW170814_skymap.fits.gz")
```
Running this will give
```output_search_GW170814
<CrossmatchResult length=3818>
     name         RA       DEC    zlens  type  grading  searched probability   searched area  
                 deg       deg                                                      deg2      
    str20      float64   float64  str15  str7    str9         float64             float64     
------------- --------- --------- ----- ------ -------- -------------------- -----------------
DESJ0303-4626   45.9507 -46.44066  1.37 galaxy probable  0.11857081625736535  2.59328622407046
DESJ0311-4232  47.86322 -42.53863  0.37 galaxy probable   0.2301796464718608 5.619333233953008
          ...       ...       ...   ...    ...      ...                  ...               ...
145459+015846 223.74648  1.979505 0.375 galaxy probable   0.9999999999999997 41252.96124941707
145612+011941 224.05141  1.328216 0.555 galaxy probable   0.9999999999999997 41252.96124941707
 
```
The cross-matching can be done to the sky localization from any type of transients as long as it is in the FITS format. For example, to cross-match the localization of GRB 240229A (download from [here](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2024/bn240229588/quicklook/glg_healpix_all_bn240229588.fit)), simply run
```python
import lenscat
lenscat.catalog.crossmatch("glg_healpix_all_bn240229588.fit")
```
In this case, the output would be
```output_search_GRB
<CrossmatchResult length=4587>
      name           RA         DEC     zlens   type   grading  searched probability   searched area   
                    deg         deg                                                         deg2       
     str20        float64     float64   str15   str7     str9         float64             float64      
--------------- ----------- ----------- ------ ------ --------- -------------------- ------------------
SDSSJ1320+1644*   200.24778    16.73437  0.899 galaxy confident   0.1614180609749184 6.9241725729921235
 SDSSJ1330+1750   202.63079    17.84456 0.2074 galaxy confident   0.6132034472687292  44.48256319619201
            ...         ...         ...    ...    ...       ...                  ...                ...
      B0128+437    22.80585    43.97032 1.145? galaxy confident                  1.0  41038.10207619734
      H1417+526 4.402083333 52.44444444   0.81 galaxy confident                  1.0  41038.10207619734
 
```

To generate a visualization of a crossmatching result, simply invoke `.plot()` to a crossmatching result. For example,
```python
import lenscat
lenscat.catalog.crossmatch("GW170814_skymap.fits.gz").plot(searched_prob_threshold=0.7)
```
will generate a figure like this
![image](https://github.com/lenscat/lenscat/assets/55488840/12d398e4-6f58-40e5-8cea-edc6bf19d701)


## Format

<table width="300">
  <tr>
    <th width="20%">Column name</th>
    <th width="80%">Description</th>
  </tr>
  <tr>
    <td><code>name</code></td>
    <td>Names of galaxies/galaxy clusters</td>
  </tr>
  <tr>
    <td><code>RA [deg]</code></td>
    <td>Right ascension in dergees</td>
  </tr>
  <tr>
     <td><code>DEC [deg]</code></td>
     <td>Declination in degress</td>
  </tr>
  <tr>
     <td><code>zlens</code></td>
     <td>Lens redshift (if known)</td>
  </tr>
  <tr>
     <td><code>type</code></td>
     <td>Type of lens (i.e. galaxy or galaxy cluster). Either<code>galaxy</code> or <code>cluster</code>.</td>
  </tr>
  <tr>
     <td><code>grading</code></td>
     <td>Grading whether it is a confidently identified lens or a probable lens (see individual references for internal grading systems). Either <code>confident</code> or <code>probable</code>.</td>
  </tr>
  <tr>
     <td><code>ref</code></td>
     <td>Reference to the corresponding catalog or study</td>
  </tr>
</table>

## References

This catalog contains the known strong lenses from the following studies:

  - [GLQ Database](https://research.ast.cam.ac.uk/lensedquasars/index.html)

  - [CLASH (Postman+2012)](https://archive.stsci.edu/prepds/clash/)

  - [MUSES Cluster Followups (Richards+2020)](https://cral-perso.univ-lyon1.fr/labo/perso/johan.richard/MUSE_data_release/)

  - [RELICS](https://relics.stsci.edu/clusters.html)

  - [37 Clusters from SDSS Giant Arcs Survey](https://iopscience.iop.org/article/10.3847/1538-4365/ab5f13)

  - [An Extended Catalog of Galaxy–Galaxy Strong Gravitational Lenses Discovered in DES Using Convolutional Neural Networks](https://iopscience.iop.org/article/10.3847/1538-4365/ab26b6#apjsab26b6t5)

  - [The AGEL Survey: Spectroscopic Confirmation of Strong Gravitational Lenses in the DES
    and DECaLS Fields Selected Using Convolutional Neural Networks](https://arxiv.org/ftp/arxiv/papers/2205/2205.05307.pdf)

  - [LSD Survey](https://web.physics.ucsb.edu/~tt/LSD/)

  - [(COSMOS) LensFlow: A Convolutional Neural Network in Search of Strong Gravitational Lenses](https://ui.adsabs.harvard.edu/abs/2018ApJ...856...68P/abstract)

  - [SLACS. XIII. Galaxy-scale strong lens candidates](https://ui.adsabs.harvard.edu/abs/2019yCat..18510048S/abstract)

  - [RINGFINDER: Automated Detection of Galaxy-scale Gravitational Lenses in Ground-based Multi-filter Imaging Data](https://iopscience.iop.org/article/10.1088/0004-637X/785/2/1440)

  - [Survey of Gravitationally-lensed Objects in HSC Imaging (SuGOHI) Candidate List](https://www-utap.phys.s.u-tokyo.ac.jp/~oguri/sugohi/)

## See also
[Master Lens Database](https://test.masterlens.org/index.php)

## Acknowledgements
This project was supported by the research grant no. VIL37766 and no. VIL53101 from Villum Fonden, and the DNRF Chair program grant no. DNRF162 by the Danish National Research Foundation.

This project has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 101131233.

We would also like to thank Jonah Kanner for introducing us the amazing [streamlit](https://streamlit.io) service that hosts the [web app](https://lenscat.streamlit.app/) for `lenscat`.
