.. _diagnostics:

Diagnostics
===========

By default, PP creates extensive diagnostic output to provide the user
with a comfortable way to assess the quality of the derived
results. Focusing on human readability, the diagnostic output is
accessible with a web browser.


Overview
--------

HTML Output
~~~~~~~~~~~

The diagnostic output for one specific data set can be accessed by
opening the ``diagnostics.html`` website that is created in that
directory, e.g., using::

  firefox diagnostics.html

Data presented in this website are all stored in a sub-directory named
``.diagnostics/`` in the respective data directory (note that the dot
prefix makes this a hidden directory). 

In case PP has been run using the `all` option on a number of
underlying data directories, it will create ``diagnostics.html`` and
``.diagnostics/`` in every sub-directory containing data. In order to
provide easy access to the individual ``diagnostics.html`` file, PP
will also create a website named ``summary.html`` in the root
directory from which is has been called, linking to all created
``diagnostics.html`` websites.


LOG File
~~~~~~~~

All tasks performed by PP are documented using Python's logging
system. The LOG file is linked from the ``diagnostics.html`` files
(see top of the page). If you run PP using the `all` option, there is
one LOG file for all data sets (unfortunately, I haven't figured out
how to change the logging filename during runtime...).

The LOG file can be used to check proper performance of the
pipeline. Its main purpose is to simplify debbuging in case something
went wrong.


Diagnostics Details
-------------------

TBD


