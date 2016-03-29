Problems?
=========

Frequently Asked Questions
--------------------------

General Problems
~~~~~~~~~~~~~~~~

**The pipeline does not run and gives me an error similar to `pp_run: Command not found`. What is wrong?**
   There are two things to check: (1) did you properly setup the
   `PHOTPIPEDIR` environment variable (check this by typing ``echo
   $PHOTPIPEDIR`` in a terminal), or (2) the command ``pp_run`` uses a
   symbolic link to ``pp_run.py``; if the former does not work, try
   the latter.

pp_calibrate (Photometric Calibration)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**I keep getting output like `zeropoint for image.fits: Warning: 0 reference stars after source matching for frame image.ldac`. What does it mean?**

   It means that none of the reference stars with measured magnitudes
   in your field of view could be matched with a source in your
   image. As a result, the magnitudes in the photometry output files
   are simply instrumental magnitudes, not calibrated ones. Try using
   a different photometric catalog in ``pp_calibrate``. If your
   field of view is small (<2 arcmin), there might just be no stars
   with known magnitudes in the field, in which case there is not a
   lot that can be done...
   

... This Does Not Solve My Problem...
-------------------------------------

If you encounter problems, e.g., PP stops unexpectedly or crashes,
please don't hesitate to contact me (michael . mommert (at) nau . edu)
and attach the following things to your email:

* the **LOG file** of your PP process (can be found in the
  ``.diagnostics/`` directory, and 

* the **error message** that was printed on the screen.

I will try to get back to you as soon as possible.
