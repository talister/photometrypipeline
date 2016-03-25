.. _quickstart:

Quickstart
==========


Prerequisits
------------

Image data should be properly reduced before using the pipeline for
best results, including cropping the data section. Bias subtraction
and flat fielding improves photometry results but is not absolutely
necessary. PP's ability to provide astrometric and photometric
calibration puts some constraints on the way data is stored: data from
separated fields, as well as data using different instrument settings
(e.g., different binning modes) should be stored in individual
directories, which in turn should be separated by filters::

  -+- all_data -+
                |
                +- field_1 -+- filter_1
                |           +- filter_2
                |           +- filter_3
                |
                +- field_2 -+- filter_1
                |           +- filter_2
               ...

Separated fields are defined as having gaps between individual frames
that are comparable to, or larger than, the field of view. Series of
frames that were tracked on a moving target can be put in the same
directory if the total track is smaller than 3-5 times the size of a
single-frame field of view.

Moving objects are currently only identified based on the images'
``OBJECT`` keywords. The object name should be as simple as possible,
consisting either of the bodies official number or designation; please
use either a blank or an underscore to separate the designation's year
from the identifier.


Running PP
----------

PP can be run in a fully automated or semi-automated mode, providing
different levels of user interaction.


Fully Automated Mode
~~~~~~~~~~~~~~~~~~~~

In the directory tree example above, PP can be run in different
places, treating the data differently. If you want to run PP only on
data for one field and filter, you can change in that directory and
**run PP locally on all fits files in that directory**::

  cd all_data/field_1/filter_1
  pp_run *fits

If your data are organized as shown in the example above, you can **run
PP from any higher level directory to analyze all underlying
directories in a consecutive way**::
  
  cd all_data
  pp_run all

Passing ``all`` signalizes PP to walk through underlying directories,
starting from the current one. In case you want PP to only run on a
subset of fits files, you can use option `-prefix`, e.g., ::

  pp_run -prefix reduced all 

is the equivalent of using only files that are included in
``reduced*.fits``.



Semi-Automated Mode
~~~~~~~~~~~~~~~~~~~


TBD


PP Diagnostics
--------------

PP generates by default significant amounts of diagnostic information
on each run. These information can be accessed in the individual
directories where the data resides with any web browser, e.g., ::

  firefox all_data/field_2/filter_3/diagnostics.html

If you ran PP with the `all` argument (see above), a file
``summary.html`` will be generated in the root directory (``all_data``),
which provides links to the individual ``index.html`` files.


More information on the diagnostic output is available here:
:ref:`diagnostics`.


Results
-------

PP derives the calibrated photometry for the target that it finds in
the ``OBJECT`` header keyword, as well as one rather bright 'control
star' that is used to check the consistency of the photometric
calibration. Results are written to files
``photometry_<objectname>.dat`` in the respective filter directory.
