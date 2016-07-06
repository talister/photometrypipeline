Functions
=========

The individual pipeline functions are introduced and explained
below. All functions presented here can be called from the terminal. 

.. function:: pp_run ([-prefix string], [-target string], [-filter string],
              [-fixed_aprad float], images)

   serves as a wrapper for all the individual pipeline processes

   :param -prefix: (optional) the prefix of all science data images if
                   `pp_run` is called using ``images = all``
   :param -target: (optional) the target name to be used in all
                   images, this overwrites the ``OBJECT`` keyword in the
                   FITS headers; note to replace blanks with
                   underscores if the target's name is a designation
   :param -filter: (optional) manual filter name override for the photometric
		   calibration step
   :param -fixed_aprad: (optional) use this fixed aperture radius for
                        all data instead of finding the aperture
                        radius in a curve-of-growth analysis
		   

   :param images: images on which the pipeline is supposed to run,
                  wildcard symbols (``'*'``, ``'?'``) can be used; or,
                  by using ``all``, PP runs on all FITS files in
                  underlying directories (the range of images can be
                  limited by using the `-prefix` option)

   The use of `pp_run` is discussed in the :ref:`quickstart` reference.

   This wrapper should work successfully for most data sets. If the
   analysis fails at some point - or provides inadequate results -
   every single step in the pipeline can be run manually, providing
   the possibility to manually tweak the process parameters.


The following functions describe the individual pipeline processes in
the logical order:


.. function:: pp_prepare ([-ra degrees], [-dec degrees], [-flipx], [-flipy], [-rotate degrees], [-target string], images)

   prepares image files for use in the pipeline

   :param -ra: (optional) manually sets the frame center RA
   :param -dec: (optional) manually sets the frame center declination
   :param -flipx: (optional) forces the image's x-axis to be flipped
                  in wcs coordinates relative to the respective
                  :ref:`telescope_setup` setting
   :param -flipy: (optional) forces the image's y-axis to be flipped
                  in wcs coordinates relative to the respective
                  :ref:`telescope_setup` setting
   :param -rotate: (optional) rotates the image's orientation in the sky
		   (East-of-North) relative to the respective
                   :ref:`telescope_setup` setting
   :param -target: (optional) the target name to be used in all
                   images, this overwrites the ``OBJECT`` keyword in the
                   FITS headers; note to replace blanks with
                   underscores if the target's name is a designation
   :param images: images on which `pp_prepare` is supposed to run

		  
   This function prepares the image data by creating necessary FITS
   header keywords (e.g., the observation midtime ``MIDTIMJD``, the
   pixel scale ``SECPIX``, ...), and by including fake wcs information
   that is required by `SCAMP`_. Existing wcs information are deleted,
   as they might cause confusion with those information generatd by
   `SCAMP`.
	     

.. function:: pp_extract ([-snr float], [-minarea integer],
              [-paramfile path], [-aprad float], [-telescope name],
              [-ignore_saturation], [-quiet], [-write_cat], images)
	      
   wrapper for `Source Extractor`_

   :param -snr: (optional) minimum SNR of sources to be extracted, default: 1.5
   :param -minarea: (optional) minimum number of connected pixels
                    above the SNR threshold for a valid source,
                    default: 3
   :param -paramfile: (optional) manual override for the `Source
                      Extractor` parameter file
   :param -aprad: (optional) aperture photometry aperture radius in
                  pixels; if no aperture radius is given here, the
                  default aperture radius for this
                  telescope/instrument combination is used (see
                  :ref:`telescope_setup` reference)
   :param -telescope: (optional) manual override for the telescope
                      identifier (see :ref:`supported_observatories`)
   :param -ignore_saturation: (optional) using this option will not
                              flag saturated sources; as a result,
                              they are not rejected in the
                              registration and calibration process
   :param -quiet: (optional) suppress output on the screen
   :param images: images to run `pp_extract` on


   `pp_extract` is automatically called by :func:`pp_register` and
   :func:`pp_photometry`. Usually, there is no reason to call this
   function manually.


.. function:: pp_register ([-snr float], [-minarea integer], [-aprad float], [-cat catalogname], images)

   astrometric calibration of the input images using `SCAMP`_ 

   :param -snr: (optional) minimum SNR of sources to be extracted for
                the registration, default: 3
   :param -minarea: (optional) minimum number of connected pixels
                    above the SNR threshold for a valid source,
                    default: :ref:`telescope_setup` setting
   :param -aprad: (optional) aperture photometry aperture radius in
                  pixels; if no aperture radius is given here, the
                  default aperture radius for this
                  telescope/instrument combination is used (see
                  :ref:`telescope_setup` reference)
   :param -cat: (optional) reference catalog override for astrometric
                calibration (a list of supported catalogs is listed
                here: :ref:`supported_catalogs`); if not specific
                catalog is requested, those listed in the
                :ref:`telescope_setup` reference are tried
   :param images: images to run `pp_register` on

   First of all, `pp_register` calls :func:`pp_extract` to identify
   sources in the field of view. ...


.. function:: pp_photometry ([-snr float], [-minarea float], [-aprad float], [-target targetname], [-background_only], [-target_only], images))

   curve-of-growth analysis of the input images and source extraction
   using the derived aperture radius resulting in final instrumental
   magnitudes

   :param -snr: (optional) minimum SNR of sources to be accounted for
                in the analysis, default: 2
   :param -minarea: (optional) minimum number of connected pixels
                    above the SNR threshold for a valid source,
                    default: :ref:`telescope_setup` setting
   :param -aprad: (optional) if this option is used, the
                  curve-of-growth analysis is skipped and instrumental
                  magnitudes are derived with this aperture radius
   :param -target: the target name to be used in all
                   images, this overrides the ``OBJECT`` keyword in the
                   FITS headers; note to replace blanks with
                   underscores if the target's name is a designation
   :param -background_only: only account for background sources in the
                            curve-of-growth analysis
   :param -target_only: only account for the target in the
                        curve-of-growth analysis
   :param image: images to run `pp_photometry` on


   Aperture finding strategies and stuff...


.. function:: pp_calibrate ([-minstars int/float], [-catalog string], [-filter string], images)

   photometric calibration of each input frame in one specific filter
   
   :param -minstars: (optional) minimum number of reference stars used
                     in the photometric calibration; if ``int``, use
                     at least this number of stars; if ``float`` use
                     at least this fraction of the available reference
                     stars; if this option is not used, the default is
                     0.5 (i.e., use at least 50% of all available
                     reference stars)
   :param -catalog: (optional) manual override for the reference
                     catalog; a list of available reference catalogs
                     is available here: :ref:`supported_catalogs`) or
                     using this routine's help function; if this
                     option is not used, the photometric reference
                     catalogs list in the :ref:`telescope_setup` are
                     used
   :param -filter: (optional) manual override for the filter used in
                     the observations; if this option is not used, the
                     filter name is read from the image FITS headers
   :param images: images to run `pp_calibrate` on

   
   Details on the calibration process....


.. function:: pp_distill ([-target string], [-offset float float], 
	      [-fixed_coo float float], images)

   extraction of calibrated photometry for targets

   :param -target: (optional) the target name to be used in all
                   images, this overrides the ``OBJECT`` keyword in the
                   FITS headers; note to replace blanks with
                   underscores if the target's name is a designation
   :param -offset: (optional) position offset to apply on target
                   positions (e.g., Horizons position for moving
                   targets) in arcsec; requires two floats, one for RA
                   and one for Dec
   :param -fixed_coo: (optional) fixed target position (requires two
                   floats, one for RA and one for Dec, coordinates in degrees)
   :param images:  images to run `pp_distill` on

   This function will automatically read the target name from the FITS
   images (or use the manually provided one), pull target positions
   from JPL Horizons, and extract calibrated photometry from the
   database catalogs created with :func:`pp_calibrate` in to a
   ``photometry_<targetname>.dat`` file. In addition to the primary
   target, this function also creates a photometry output file for one
   relatively bright star that is present in the first and the last
   image of the series - this star serves as a control star to check
   the consistency of the derive magnitude zeropoints.

   


.. _Source Extractor: http://www.astromatic.net/software/sextractor
.. _SCAMP: http://www.astromatic.net/software/scamp

