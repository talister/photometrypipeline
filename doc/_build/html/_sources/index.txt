.. Photometry Pipeline documentation master file, created by
   sphinx-quickstart on Mon Mar  7 11:53:14 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Photometry Pipeline Documentation
=================================

Introduction
------------

The Photometry Pipeline (PP) is a Python 2.7 software package for
automated photometric analysis of imaging data from small to
medium-sized observatories. It uses `Source Extractor`_ and `SCAMP`_
to register and photometrically calibrate images based on catalogs
that are available online; photometry is measured using Source
Extractor aperture photometry. PP has been designed for asteroid
observations, but can be used with any kind of imaging data.


Scope
-----

PP has been designed to provide automated photometry for the majority
of data coming from small to medium-sized observatories. It is not
intended to provide high-accuracy photometry, nor is it designed to
work on extremely sparse or crowded fields. PP requires a field of
view of a few arcminutes to ensure that sufficient background stars
are available for registration and photometric calibration. For a
7'x7' field, photometric uncertaintes are of the order of 0.03 mag (if
sufficient SDSS-R9 or URAT-1 stars are in the field). Feel free to try
PP on your data, but please be aware that it has its limits.



Contents
--------

.. toctree::
   :maxdepth: 2

   install
   quickstart
   supported
   functions
   diagnostics
   problems


License and Contact
-------------------

The Photometry Pipeline is distributed under the GNU GPLv3 license.

Copyright (C) 2016  Michael Mommert 

Feel free to contact me in case of questions or suggestions: michael
(at) mommert . nau . edu



Acknowledgments
---------------

PP is supported by NASA grants NNX15AE90G and NNX14AN82G and has been
developed in the framework of the Mission Accessible Near-Earth
Objects Survey (`MANOS`_).


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _Source Extractor: http://www.astromatic.net/software/sextractor
.. _SCAMP: http://www.astromatic.net/software/scamp
.. _MANOS: http://manosobs.wordpress.com/
