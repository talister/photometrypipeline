Installation and Setup
======================

Installation
------------

Getting PP
~~~~~~~~~~

PP is available from `github`_. You can get the source code by typing
into your terminal::

  git clone https://github.com/mommermi/photometrypipeline

This will create a ``photometrypipeline/`` directory in your current
directory. 

Installing Additional Software
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PP requires `git`_ for the installation, a number of non-standard
Python modules (available from the `Python Package Index`_ through
`pip`_):

* `numpy`_
* `scipy`_
* `astropy`_
* `matplotlib`_
* `callhorizons`_

and some freely available software:

* `Source Extractor`_ (requires the latest development version r345 to
  use URAT-1)
* `SCAMP`_  


Setup
-----

In order to be able to use PP anywhere on your machine, you
have to add ``photometrypipeline/`` to your `PYTHONPATH` variable, and
you have to create a `PHOTPIPEDIR` variable on your system that points
to the same directory (include both commands in your ``.bashrc``,
``.cshrc``, or ``.profile`` file.)

In order to get the latest version of PP, simply change into
``photometrypipeline/`` and type::

  git pull


.. _telescope_setup:

Telescope Setup
~~~~~~~~~~~~~~~

all the stuff that goes into ``setup/telescopes.py``...


.. _github: https://github.com/mommermi/photometrypipeline
.. _git: http://www.git-scm.com/
.. _Python Package Index: https://pypi.python.org/pypi
.. _pip: https://pypi.python.org/pypi/pip/
.. _numpy: http://www.numpy.org/
.. _scipy: https://www.scipy.org/
.. _astropy: http://www.astropy.org/
.. _matplotlib: http://matplotlib.org/
.. _callhorizons: https://pypi.python.org/pypi/CALLHORIZONS
.. _Source Extractor: http://www.astromatic.net/software/sextractor
.. _SCAMP: http://www.astromatic.net/software/scamp
