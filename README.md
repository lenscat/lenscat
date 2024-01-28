# $\texttt{lenscat}$
![license](https://img.shields.io/github/license/lenscat/lenscat)
[![GitHub release](https://img.shields.io/github/v/release/lenscat/lenscat.svg)](https://github.com/lenscat/lenscat/releases)
[![Upload Python Package](https://github.com/lenscat/lenscat/actions/workflows/python-publish.yml/badge.svg)](https://github.com/lenscat/lenscat/actions/workflows/python-publish.yml)

A public and community-contributed catalog of known strong gravitational lenses. 

![Known Lenses](https://raw.githubusercontent.com/lenscat/lenscat/main/catalog.png)

## Quickstart

The catalog is available as a _plain csv file_ under [lenscat/data/catalog.csv](https://github.com/lenscat/lenscat/blob/main/lenscat/data/catalog.csv).

We also provide a _python package_ `lenscat`, available in pypi. Simply do
```bash
pip install lenscat
```
to install the latest version. The code converts the catalog in the csv file into a `Table` object in `astropy`. To access the table, simply run
```python
In [1]: import lenscat; lenscat.catalog
Out[1]:
<Table length=1530>
     name          RA        DEC      zlens    type   grading                                     ref
                  deg        deg
    str20       float64    float64    str15    str7     str9                                     str76
-------------- ---------- ---------- -------- ------ --------- --------------------------------------------------------------------------
    J0011-0845    2.83435    -8.7643        - galaxy confirmed                    https://research.ast.cam.ac.uk/lensedquasars/index.html
    J0013+5119   3.348077    51.3183        - galaxy confirmed                    https://research.ast.cam.ac.uk/lensedquasars/index.html
  PSJ0028+0631    7.09369     6.5317        - galaxy confirmed                    https://research.ast.cam.ac.uk/lensedquasars/index.html
  PSJ0030-1525     7.5636   -15.4177 measured galaxy confirmed                    https://research.ast.cam.ac.uk/lensedquasars/index.html
    J0047+2514 11.9465943    25.2411        - galaxy confirmed                    https://research.ast.cam.ac.uk/lensedquasars/index.html
   HE0047-1756    12.6158   -17.6693    0.407 galaxy confirmed                    https://research.ast.cam.ac.uk/lensedquasars/index.html
 DESJ0053-2012    13.4349   -20.2091 observed galaxy confirmed                    https://research.ast.cam.ac.uk/lensedquasars/index.html
    J0102+2445   15.69675    24.7544   0.272? galaxy confirmed                    https://research.ast.cam.ac.uk/lensedquasars/index.html
```

The code also implements `crossmatch()`, a wrapper to the `crossmatch()` function in `ligo.skymap`, to cross-match a gravitational-wave (GW) skymap with the catalog. For example, to cross-match the GW skymap of GW190814 (download from [here](https://gracedb.ligo.org/apiweb/superevents/S190814bv/files/GW190814_PublicationSamples.multiorder.fits)) with the catalog, simply run
```python
In [1]: import lenscat; lenscat.crossmatch("GW190814_PublicationSamples.multiorder.fits")
Out[1]:
<Table length=1530>
       name             RA         DEC      ... searched probability   searched area
                       deg         deg      ...                             deg2
      str20          float64     float64    ...       float64             float64
------------------ ----------- ------------ ... -------------------- ------------------
     DESJ0133-3137    23.36279    -31.61788 ...   0.8993434968830626 18.766081348393186
  MACSJ0035.4-2015        8.85 -20.27083333 ...   0.9881400284668479  46.21360633943737
          0047-281 12.42458333 -27.87380556 ...   0.9899384024067811  48.32495441567177
             A2813     10.8625 -20.61694444 ...   0.9924501903529119  51.82638259178721
     DESJ0102-2911    15.73954    -29.18939 ...    0.993395993986835   53.3476023237325
     DESJ0124-2918    21.11886    -29.31561 ...   0.9951314306436486  56.65232105175158
     DESJ0058-2317    14.52023    -23.28713 ...    0.997642284491974  63.57649362474393
     DESJ0113-2924    18.48794     -29.4109 ...   0.9996945313931366  77.47729462355434
```

## Format

['name'] = Names of galaxies/galaxy clusters \
['RA [deg]'] = RA in dergees \
['DEC [deg]'] = DEC in degrees \
['zlens'] = Lens redshift (if known) \
['type'] = Type of lens (i.e. galaxy or galaxy cluster) \
['grading'] = Grading whether it is a confirmed lens or a probable lens (see individual references for internal ranking systems) \
['ref'] = Reference to the corresponding catalogue or study

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
