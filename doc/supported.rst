.. _supported_observatories:

Supported Observatories
~~~~~~~~~~~~~~~~~~~~~~~

The pipeline is currently set up to work on data from the following
observatories/instruments:

+--------------------------+--------------------+----------------+
| Telescope                | Instrument         | PP Keyword     |
+==========================+====================+================+
| Apache Point ARC 3.5m    | AGILE              | ARC35AGILE     |
+--------------------------+--------------------+----------------+
| Apache Point ARC 3.5m    | ARCTIC             | ARC35ARCTIC    |
+--------------------------+--------------------+----------------+
| Calar Alto 1.23m         | DLR-MkIII          | CA123DLRMKIII  |
+--------------------------+--------------------+----------------+
| CTIO 0.9m                | CFCCD              | CTIO09         |
+--------------------------+--------------------+----------------+
| CTIO 1.0m                | Y4KCam             | CTIO10         |
+--------------------------+--------------------+----------------+
| CTIO 1.3m                | ANDICAM (CCD)      | CTIO13CCD      |
+--------------------------+--------------------+----------------+
| Discovery Channel        | Large Monolithic   | DCTLMI         |
| Telescope                | Imager             |                |
+--------------------------+--------------------+----------------+
| Lowell 31"               | NASACAM            | LOWELL31       |
+--------------------------+--------------------+----------------+
| Lowell 42"               | NASA42             | LOWELL42       |
+--------------------------+--------------------+----------------+
| Magellan                 | IMACS              | MAGIMACS       |
+--------------------------+--------------------+----------------+
| Vatican Advanced         | VATT4k             | VATT4K         |
| Technology Telescope     |                    |                |
+--------------------------+--------------------+----------------+
| WIYN 0.9m                | Half Degree Imager | WIYN09HDI      |
+--------------------------+--------------------+----------------+

If you would like to use the pipeline for other observatories, please
contact me.

.. _supported_catalogs:

Supported Reference Catalogs 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PP is currently able to use the following catalogs for astrometric and
photometric calibration:

+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| Catalog Name           | Catalog Type             | Registration? | Filter Bands             | Comments                                                   |
+========================+==========================+===============+==========================+============================================================+
| 2MASS (`2MASS`_)       | astrometric/photometric  | yes           | J, H, Ks, K*             | all-sky NIR catalog, good astrometry                       |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| URAT-1 (`URAT-1`_)     | astrometric/photometric  | yes (SCAMP    | g, r, i,  B, V, R*,      | good coverage over the Northern hemisphere, photometry from|
|                        |                          | >= trunk.r345)| I*                       | APASS (see below)                                          |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| Sloan Digital Sky      | astrometric/photometric  | yes           | u, g, r, i, z, U*, B*,   | excellent photometry, Northern hemisphere, patchy coverage |
| Survey Release 9       |                          |               | V*, R*, I*               |                                                            | 
| (`SDSS-R9`_)           |                          |               |                          |                                                            |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| AAVSO Photometric All  | photometric              | no            | g, r, i, B, V, R*,       | good coverage, good photometry for stars with V<17         | 
| Sky Survey Release 9   |                          |               | I*                       |                                                            |
| (`APASS9`_)            |                          |               |                          |                                                            |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| Carlsberg Meridian     | photometric              | no            | r, J, H, Ks              | -40 deg < declination < 50 deg; J, H, Ks from 2MASS (above)|
| Catalogue 15 (`CMC15`_)|                          |               |                          | limiting magnitude r < 17 mag                              |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| PPMXL (`PPMXL`_)       | astrometric              | yes           | J, H, Ks, B, R           | accurate astrometry; J, H, Ks from 2MASS; B, R from USNO-B1|
|                        |                          |               |                          | (see below). Note that proper motions are not yet taken    |
|                        |                          |               |                          | into account (this will be integrated in the future).      |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| USNO-B1 (`USNO-B1`_)   | astrometric              | yes           | B, R, I                  | questionable photometric and astrometric accuracy, but     |
|                        |                          |               |                          | deepest catalog available (V < 21 mag)                     |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+

The catalog name in brackets is the identifier used by PP; e.g., if
you want to use URAT-1 for the registration of your images, use option
`-cat URAT-1` with the ``pp_register`` command. The "Registration?"
column identifies if this catalog is supported by SCAMP in order to
use it for image registration. All filter bands marked with an
asterisk (*) are obtained through transformations from other bands.


If you are interested in using catalogs other than those listed,
please let me know.


Supported Catalog Transformations
---------------------------------

PP supports the following catalog transformations:

* ugriz -> BVRI: `Chonis & Gaskell 2008`_
* JHKs (2MASS) -> JHK (UKIRT): `Hodgkin et al. 2009`_

Independent checks indicate that these transformations are reliable and accurate. More quantitative results coming soon...


.. _Chonis & Gaskell 2008: http://adsabs.harvard.edu/abs/2008AJ....135..264C
.. _Hodgkin et al. 2009: http://adsabs.harvard.edu/abs/2009MNRAS.394..675H


.. _2MASS: http://www.ipac.caltech.edu/2mass/
.. _URAT-1: http://cdsads.u-strasbg.fr/cgi-bin/nph-bib_query?2015AJ....150..101Z&db_key=AST&nosetcookie=1
.. _SDSS-R9: http://www.sdss3.org/dr9/
.. _APASS9: http://www.aavso.org/apass
.. _CMC15: http://cdsarc.u-strasbg.fr/viz-bin/Cat?I/327
.. _PPMXL: http://cdsads.u-strasbg.fr/cgi-bin/nph-bib_query?2010AJ....139.2440R&db_key=AST&nosetcookie=1
.. _USNO-B1: http://cdsads.u-strasbg.fr/cgi-bin/nph-bib_query?2003AJ....125..984M&db_key=AST&nosetcookie=1
