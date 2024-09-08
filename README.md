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
     name         RA        DEC      zlens     type   grading 
                 deg        deg                               
    str20      float64    float64    str15     str7     str9  
------------- ---------- ---------- -------- ------- ---------
   J0011-0845    2.83435    -8.7643        -  galaxy confident
   J0013+5119   3.348077    51.3183        -  galaxy confident
 PSJ0028+0631    7.09369     6.5317        -  galaxy confident
 PSJ0030-1525     7.5636   -15.4177 measured  galaxy confident
   J0047+2514 11.9465943    25.2411        -  galaxy confident
  HE0047-1756    12.6158   -17.6693    0.407  galaxy confident
DESJ0053-2012    13.4349   -20.2091 observed  galaxy confident
   J0102+2445   15.69675    24.7544   0.272?  galaxy confident
DESJ0112-1650  18.141188 -16.840955     0.54  galaxy confident
          ...        ...        ...      ...     ...       ...
235428-001256  358.61846   -0.21578    0.262  galaxy  probable
235503+014443 358.765312   1.745288    0.522  galaxy  probable
235535+021044   358.8969     2.1791    0.519 cluster  probable
235614+023115  359.05923    2.52107    0.372  galaxy  probable
235730+010133 359.377643    1.02596    0.638  galaxy  probable
235811+003309   359.5475     0.5527    0.639 cluster  probable
235853+012406 359.721708   1.401816    0.481  galaxy  probable
235933+020823   359.8897     2.1398     0.43 cluster confident
235948-005913  359.95245   -0.98702    0.758  galaxy  probable
235952+004154   359.9698     0.6985    0.267 cluster  probable
 
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
     name         RA       DEC    zlens   type   grading                                        ref                                       
                 deg       deg                                                                                                            
    str20      float64   float64  str15   str7     str9                                        str171                                     
------------- --------- --------- ----- ------- --------- --------------------------------------------------------------------------------
021118-042729 32.827087 -4.458069  1.02 cluster confident https://arxiv.org/abs/2004.00634 More et al. 2012[FIXME] More et al. 2016[FIXME]
023100-062139   37.7516   -6.3608  1.17 cluster confident                                                 https://arxiv.org/abs/2002.01611
220859+020655  332.2495    2.1153  1.04 cluster confident                                                 https://arxiv.org/abs/2002.01611
 
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
     name          RA        DEC    zlens   type   grading  searched probability   searched area   
                  deg        deg                                                        deg2       
    str20       float64    float64  str15   str7     str9         float64             float64      
-------------- ---------- --------- ------ ------ --------- -------------------- ------------------
 DESJ0303-4626    45.9507 -46.44066   1.37 galaxy  probable  0.11857081625736535   2.59328622407046
 DESJ0311-4232   47.86322 -42.53863   0.37 galaxy  probable   0.2301796464718608  5.619333233953008
 DESJ0310-4647   47.63526 -46.78398   0.71 galaxy  probable  0.36778302134840013 10.261676209026643
 DESJ0301-4426    45.4638 -44.44055   0.76 galaxy  probable   0.4683381098641989 14.661410864782189
 DESJ0304-4921   46.06729 -49.35725   0.34 galaxy confident   0.6465740359340766 26.791826830724982
 DESJ0300-5001   45.09019 -50.02469   0.53 galaxy confident   0.7082286031333002 33.860252998986475
 DESJ0320-4159   50.11682 -41.98881   0.66 galaxy  probable   0.7333678230113636 37.689530255261715
 DESJ0235-4818   38.86317 -48.30583   0.47 galaxy  probable   0.7502684940168021   40.6008300870876
 DESJ0242-4811   40.68463 -48.19355   0.48 galaxy  probable   0.7899164033300675  48.52166386376813
           ...        ...       ...    ...    ...       ...                  ...                ...
 145130+042054 222.878531  4.348363    0.7 galaxy  probable   0.9999999999999997  41252.96124941707
 150408+015744 226.037262  1.962391  0.525 galaxy  probable   0.9999999999999997  41252.96124941707
 145402+013809  223.51216  1.636013  0.537 galaxy  probable   0.9999999999999997  41252.96124941707
SDSSJ1425+0951  216.39804      9.85 0.1583 galaxy  probable   0.9999999999999997  41252.96124941707
SDSSJ1531+0652  232.80183   6.87033 0.2085 galaxy  probable   0.9999999999999997  41252.96124941707
ULASJ1529+1038  232.41203  10.63454   0.4? galaxy confident   0.9999999999999997  41252.96124941707
SDSSJ1455+1447  223.75782  14.79301  0.42* galaxy confident   0.9999999999999997  41252.96124941707
SDSSJ1515+1511   228.9106  15.19309  0.742 galaxy confident   0.9999999999999997  41252.96124941707
 145459+015846  223.74648  1.979505  0.375 galaxy  probable   0.9999999999999997  41252.96124941707
 145612+011941  224.05141  1.328216  0.555 galaxy  probable   0.9999999999999997  41252.96124941707
 
```
The cross-matching can be done to the sky localization from any type of transients as long as it is in the FITS format. For example, to cross-match the localization of GRB 240229A (download from [here](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2024/bn240229588/quicklook/glg_healpix_all_bn240229588.fit)), simply run
```python
import lenscat
lenscat.catalog.crossmatch("glg_healpix_all_bn240229588.fit")
```
In this case, the output would be
```output_search_GRB
<CrossmatchResult length=4587>
      name           RA          DEC     zlens    type   grading  searched probability   searched area   
                    deg          deg                                                          deg2       
     str20        float64      float64   str15    str7     str9         float64             float64      
--------------- ------------ ----------- ------ ------- --------- -------------------- ------------------
SDSSJ1320+1644*    200.24778    16.73437  0.899  galaxy confident   0.1614180609749184 6.9241725729921235
 SDSSJ1330+1750    202.63079    17.84456 0.2074  galaxy confident   0.6132034472687292  44.48256319619201
 SDSSJ1304+2001    196.18166    20.01805   0.4?  galaxy confident   0.6545106150094973  51.40673576918417
 SDSSJ1330+1810    202.57772    18.17581  0.373  galaxy confident   0.6730233044611307  54.97373376133156
 SDSSJ1258+1657    194.58017    16.95489   0.4?  galaxy confident   0.6890963497394108  58.33090834217615
 SDSSJ1254+1857     193.6681    18.95333  0.555  galaxy confident   0.7941818270703923  90.22406686019981
 SDSSJ1322+1052    200.65152     10.8778  0.55?  galaxy confident   0.8414948789612942 114.35375916002097
 SDSSJ1339+1310    204.77999     13.1775  0.607  galaxy confident   0.8806457604955152 143.09956650850438
SDSS J1329+2243  202.3916667 22.71673333 0.4427 cluster confident   0.8841570484842778 146.24691767804637
            ...          ...         ...    ...     ...       ...                  ...                ...
 SDSSJ0746+4403      116.721    44.06425  0.513  galaxy confident                  1.0 40823.242903023274
 SDSSJ0754+4202    118.50279    42.04672 0.3692  galaxy  probable                  1.0 40823.242903023274
 SDSSJ0753+3839    118.27954    38.66086 0.0408  galaxy  probable                  1.0 40823.242903023274
 SDSSJ0821+4542     125.4943    45.71233  0.349  galaxy confident                  1.0 40823.242903023274
      Abell 611  120.2367917 36.05669444  0.288 cluster confident                  1.0 40823.242903023274
MACS0744.9+3927       116.22 39.45677778  0.686 cluster confident                  1.0 40823.242903023274
      H1543+535  10.83708333 53.86444444  0.497  galaxy confident                  1.0  41038.10207619734
      Q0957+561 0.3366666667 55.89705556  0.356  galaxy confident                  1.0  41038.10207619734
      B0128+437     22.80585    43.97032 1.145?  galaxy confident                  1.0  41038.10207619734
      H1417+526  4.402083333 52.44444444   0.81  galaxy confident                  1.0  41038.10207619734
 
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
