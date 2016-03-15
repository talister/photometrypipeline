Photometry Pipeline 
===================

The Photometry Pipeline (PP) is a Python software package for
automated photometric analysis of imaging data from small to
medium-sized observatories. It uses `Source Extractor`_ and `SCAMP`_ to
register and photometrically calibrate images based on catalogs that
are available online; photometry is measured using Source Extractor
aperture photometry. PP has been designed for asteroid observations,
but can be used with any kind of imaging data.


Please note that this project is still under development. 

Currently Supported Observatories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* CTIO/ANDICAM (CCD camera)
* DCT/LMI
* VATT/VATT4k


Currently Supported Reference Catalogs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* SDSS-R9
* URAT-1
* 2MASS



Installation
------------

PP requires `git`_ for the installation, a number of non-standard
Python modules:

* `numpy`_
* `scipy`_
* `astropy`_
* `matplotlib`_
* `callhorizons`_

and some freely available software:

* `Source Extractor`_ (requires the latest development version r345 to
  use URAT-1)
* `SCAMP`_  

PP is available from `github`_. You can get the source code by typing
into your terminal::

  git clone https://github.com/mommermi/photometrypipeline


Documentation
-------------

See the documentation for more information: `documentation`_


License and Contact
-------------------

The Photometry Pipeline is distributed under the GNU GPLv3 license.

Copyright (C) 2016  Michael Mommert 

Feel free to contact me in case of questions or suggestions: michael
(at) mommert . nau . edu


.. _github: https://github.com/mommermi/photometrypipeline
.. _git: http://www.git-scm.com/
.. _numpy: http://www.numpy.org/
.. _scipy: https://www.scipy.org/
.. _astropy: http://www.astropy.org/
.. _matplotlib: http://matplotlib.org/
.. _callhorizons: https://pypi.python.org/pypi/CALLHORIZONS
.. _Source Extractor: http://www.astromatic.net/software/sextractor
.. _SCAMP: http://www.astromatic.net/software/scamp
.. _documentation: http://mommermi.github.io/pp/index.html


