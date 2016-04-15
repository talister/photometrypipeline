"""
Photometry Pipeline Configuation File
2016-03-09, michael.mommert@nau.edu
"""

# Photometry Pipeline 
# Copyright (C) 2016  Michael Mommert, michael.mommert@nau.edu

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.



##### telescope/instrument configurations

# VATT, VATT4k
vatt4k_param = {
    'telescope_instrument' : 'VATT/VATT4k', # telescope/instrument name
    'telescope_keyword'    : 'VATT4k',      # telescope/instrument keyword
    'observatory_code'     : '290',         # MPC observatory code
    'secpix'               : (0.1875, 0.1875), # pixel size (arcsec)
                                               # before binning
    'ext_coeff'            : 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx'                : True,
    'flipy'                : False,
    'rotate'               : 0,

    # instrument-specific FITS header keywords
    'binning'              : ('CCDBIN1', 'CCDBIN2'), # binning in x/y
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec 
    'radec_separator'      : ':',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS|TIME-OBS', # obs date/time
                                                  # keyword; use
                                                  # 'date|time' if
                                                  # separate
    'obsmidtime_jd'        : 'MIDTIMJD', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJECT',  # object name keyword 
    'filter'               : 'FILTER',  # filter keyword
    'filter_translations'  : {'TOP 2 BOT 1': 'V', 'TOP 3 BOT 1': 'R', 
                              'TOP 4 BOT 1': None, 'TOP 5 BOT 1': 'B'},
                             # filtername translation dictionary
    'exptime'              : 'EXPTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword


    # source extractor settings
    'source_minarea'       : 12, # default sextractor source minimum N_pixels
    'aprad_default'        : 5, # default aperture radius in px 
    'aprad_range'          : [2, 10], # [minimum, maximum] aperture radius (px)
    'sex-config-file'      : rootpath+'/setup/vatt4k.sex',
    'mask_file'            : {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file'    : rootpath+'/setup/vatt4k.scamp', 

    # swarp settings
    'copy_keywords'        : ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                              'DATE-OBS,TIME-OBS,RA,DEC,SECPIX,AIRMASS,' +
                              'TEL_KEYW'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file'    : rootpath+'/setup/vatt4k.swarp',  

    # default catalog settings
    'astrometry_catalogs'  : ['URAT-1', '2MASS'], 
    'photometry_catalogs'  : ['SDSS-R9', 'URAT-1', '2MASS'],
}





# CTIO, ANDICAM (CCD)
andicam_param = {
    'telescope_instrument' : 'CTIO/ANDICAM_CCD', # telescope/instrument name
    'telescope_keyword'    : 'ANDICAM',          # telescope/instrument keyword
    'observatory_code'     : '807',         # MPC observatory code
    'secpix'               : (0.185, 0.185), # pixel size (arcsec)
                                               # before binning
    'ext_coeff'            : 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx'                : True, 
    'flipy'                : False, 
    'rotate'               : 0, 

    # instrument-specific FITS header keywords
    'binning'              : ('CCDXBIN', 'CCDYBIN'), # binning in x/y
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec 
    'radec_separator'      : ':',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS|TIME-OBS', # obs date/time
                                                  # keyword; use
                                                  # 'date|time' if
                                                  # separate
    'obsmidtime_jd'        : 'MIDTIMJD', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJECT',  # object name keyword 
    'filter'               : 'CCDFLTID',  # filter keyword
    'filter_translations'  : {'V': 'V', 'R': 'R', 'B': 'B'},
                             # filtername translation dictionary
    'exptime'              : 'EXPTIME', # exposure time keyword (s)
    'airmass'              : 'SECZ', # airmass keyword


    # source extractor settings
    'source_minarea'       : 15, # default sextractor source minimum N_pixels
    'aprad_default'        : 4, # default aperture radius in px 
    'aprad_range'          : [2, 10], # [minimum, maximum] aperture radius (px)
    'sex-config-file'      : rootpath+'/setup/andicam.sex',
    'mask_file'            : {'2,2' : rootpath+'/setup/mask_andicam_2x2.fits'},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file'    : rootpath+'/setup/andicam.scamp', 

    # swarp settings
    'copy_keywords'        : ('OBSERVAT,INSTRUME,CCDFLTID,EXPTIME,OBJECT,' +
                              'DATE-OBS,TIME-OBS,RA,DEC,SECPIX,SECZ,' +
                              'TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file'    : rootpath+'/setup/andicam.swarp',  

    # default catalog settings
    'astrometry_catalogs'  : ['URAT-1', '2MASS'], 
    'photometry_catalogs'  : ['SDSS-R9', 'URAT-1', '2MASS'],
}





# DCT, LMI
dctlmi_param = {
    'telescope_instrument' : 'DCT/LMI', # telescope/instrument name
    'telescope_keyword'    : 'DCTLMI',  # telescope/instrument keyword
    'observatory_code'     : 'G37',         # MPC observatory code
    'secpix'               : (0.12, 0.12 ), # pixel size (arcsec)
                                            # before binning
    'ext_coeff'            : 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx'                : True, 
    'flipy'                : False, 
    'rotate'               : 0, 

    # instrument-specific FITS header keywords
    'binning'              : ('CCDSUM_blank0', 'CCDSUM_blank1'), 
                           # binning in x/y, '_blankN' denotes that both axes
                           # are listed in one keyword, sep. by blanks
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec 
    'radec_separator'      : ':',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS', # obs date/time
                                         # keyword; use
                                         # 'date|time' if
                                         # separate
    'obsmidtime_jd'        : 'MIDTIMJD', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJECT',  # object name keyword 
    'filter'               : 'FILTERS',  # filter keyword
    'filter_translations'  : {'V': 'V', 'R': 'R', 'B': 'B', 'VR': None},
                             # filtername translation dictionary
    'exptime'              : 'EXPTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword


    # source extractor settings
    'source_minarea'       : 9, # default sextractor source minimum N_pixels
    'aprad_default'        : 4, # default aperture radius in px 
    'aprad_range'          : [2, 10], # [minimum, maximum] aperture radius (px)
    'sex-config-file'      : rootpath+'/setup/dctlmi.sex',
    'mask_file'            : {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file'    : rootpath+'/setup/dctlmi.scamp', 

    # swarp settings
    'copy_keywords'        : ('OBSERVAT,INSTRUME,CCDFLTID,EXPTIME,OBJECT,' +
                              'DATE-OBS,RA,DEC,SCALE,AIRMASS,TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file'    : rootpath+'/setup/dctlmi.swarp',  

    # default catalog settings
    'astrometry_catalogs'  : ['URAT-1', '2MASS'], 
    'photometry_catalogs'  : ['SDSS-R9', 'URAT-1', '2MASS'],
}



# Apache Point ARC 3.5m, ARCTIC
arc35arctic_param = {
    'telescope_instrument' : 'ARC35m/ARCTIC', # telescope/instrument name
    'telescope_keyword'    : 'ARC35ARCTIC',   # telescope/instrument keyword
    'observatory_code'     : '705',         # MPC observatory code
    'secpix'               : (0.115, 0.115 ), # pixel size (arcsec)
                                            # before binning
    'ext_coeff'            : 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx'                : True, 
    'flipy'                : False, 
    'rotate'               : 0, 

    # instrument-specific FITS header keywords
    'binning'              : ('CCDBIN1', 'CCDBIN2'), 
                           # binning in x/y, '_blankN' denotes that both axes
                           # are listed in one keyword, sep. by blanks
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec 
    'radec_separator'      : ':',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS', # obs date/time
                                         # keyword; use
                                         # 'date|time' if
                                         # separate
    'obsmidtime_jd'        : 'MIDTIMJD', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJNAME',  # object name keyword 
    'filter'               : 'FILTER',  # filter keyword
    'filter_translations'  : {'SDSS U':'u', 'SDSS G':'g', 'SDSS R': 'r',
                              'SDSS I': 'i', 'SDSS Z': 'z'},
                             # filtername translation dictionary
    'exptime'              : 'EXPTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword


    # source extractor settings
    'source_minarea'       : 10, # default sextractor source minimum N_pixels
    'aprad_default'        : 4, # default aperture radius in px 
    'aprad_range'          : [2, 10], # [minimum, maximum] aperture radius (px)
    'sex-config-file'      : rootpath+'/setup/arc35arctic.sex',
    'mask_file'            : {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file'    : rootpath+'/setup/arc35arctic.scamp', 

    # swarp settings
    'copy_keywords'        : ('OBSERVAT,INSTRUME,FILTER,EXPTIME,OBJNAME,' +
                              'DATE-OBS,RA,DEC,AIRMASS,SECPIX,TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file'    : rootpath+'/setup/arc35arctic.swarp',  

    # default catalog settings
    'astrometry_catalogs'  : ['URAT-1', '2MASS'], 
    'photometry_catalogs'  : ['SDSS-R9', 'URAT-1', '2MASS'],
}


# Apache Point ARC 3.5m, AGILE
arc35agile_param = {
    'telescope_instrument' : 'ARC35m/AGILE', # telescope/instrument name
    'telescope_keyword'    : 'ARC35AGILE',   # telescope/instrument keyword
    'observatory_code'     : '705',         # MPC observatory code
    'secpix'               : (0.13, 0.13 ), # pixel size (arcsec)
                                            # before binning
    'ext_coeff'            : 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx'                : True, 
    'flipy'                : False, 
    'rotate'               : 0, 

    # instrument-specific FITS header keywords
    'binning'              : ('BINX', 'BINY'), 
                           # binning in x/y, '_blankN' denotes that both axes
                           # are listed in one keyword, sep. by blanks
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec 
    'radec_separator'      : ':',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS', # obs date/time
                                         # keyword; use
                                         # 'date|time' if
                                         # separate
    'obsmidtime_jd'        : 'MIDTIMJD', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJNAME',  # object name keyword 
    'filter'               : 'FILTER',  # filter keyword
    'filter_translations'  : {'SDSS u':'u', 'SDSS g':'g', 'SDSS r': 'r',
                              'SDSS i': 'i', 'SDSS z': 'z'},
                             # filtername translation dictionary
    'exptime'              : 'EXPTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword


    # source extractor settings
    'source_minarea'       : 7, # default sextractor source minimum N_pixels
    'aprad_default'        : 4, # default aperture radius in px 
    'aprad_range'          : [2, 10], # [minimum, maximum] aperture radius (px)
    'sex-config-file'      : rootpath+'/setup/arc35agile.sex',
    'mask_file'            : {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file'    : rootpath+'/setup/arc35agile.scamp', 

    # swarp settings
    'copy_keywords'        : ('OBSERVAT,INSTRUME,FILTER,EXPTIME,OBJNAME,' +
                              'DATE-OBS,RA,DEC,AIRMASS,SECPIX,TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file'    : rootpath+'/setup/arc35agile.swarp',  

    # default catalog settings
    'astrometry_catalogs'  : ['URAT-1', '2MASS'], 
    'photometry_catalogs'  : ['SDSS-R9', 'URAT-1', '2MASS'],
}


# Magellan, IMACS
magimacs_param = {
    'telescope_instrument' : 'Magellan/IMACS', # telescope/instrument name
    'telescope_keyword'    : 'MAGIMACS',   # telescope/instrument keyword
    'observatory_code'     : '269',         # MPC observatory code
    'secpix'               : (0.11, 0.11 ), # pixel size (arcsec)
                                            # before binning
    'ext_coeff'            : 0.05,          # typical extinction coefficient


    # image orientation preferences (for each chip)
    'chip_id'              : 'CHIP',        # chip identifier (remove,
                                            # if not existent) 
    # the following keys are dictionaries if 'chip_id' exists, single
    # values otherwise
    'flipx'                : {1:True, 2:True, 3:True, 4:True, 5:True, 6:True,
                              7:True, 8:True}, 
    'flipy'                : {1:False, 2:False, 3:False, 4:False, 5:False,
                              6:False, 7:False, 8:False}, 
    'rotate'               : {1:270, 2:270, 3:270, 4:270, 5:90, 6:90,
                              7:90, 8:90}, 
    'chip_offset_fixed'    : {1:(-0.033, -0.099), 2:(-0.033, -0.033),
                              3:(-0.033, 0.033),  4:(-0.033, 0.099),
                              5:(0.033, -0.033),  6:(0.033, -0.099),
                              7:(0.033, 0.099),   8:(0.033, 0.033)},
                             # chip offset (ra, dec in degress) [optional]

    # instrument-specific FITS header keywords
    'binning'              : ('BINNING_x1', 'BINNING_x2'), 
                           # binning in x/y, '_blankN' denotes that both axes
                           # are listed in one keyword, sep. by blanks
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec 
    'radec_separator'      : ':',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS|UT-TIME', # obs date/time
                                         # keyword; use
                                         # 'date|time' if
                                         # separate
    'obsmidtime_jd'        : 'MIDTIMJD', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJECT',  # object name keyword 
    'filter'               : 'FILTER',  # filter keyword
    'filter_translations'  : {'Sloan_u': 'u', 'Sloan_g': 'g', 'Sloan_r': 'r',
                              'Sloan_i': 'i', 'Sloan_z': 'z'},
                             # filtername translation dictionary
    'exptime'              : 'EXPTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword


    # source extractor settings
    'source_minarea'       : 20, # default sextractor source minimum N_pixels
    'aprad_default'        : 8, # default aperture radius in px 
    'aprad_range'          : [5, 25], # [minimum, maximum] aperture radius (px)
    'sex-config-file'      : rootpath+'/setup/magimacs.sex',
    'mask_file'            : {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file'    : rootpath+'/setup/magimacs.scamp', 

    # swarp settings
    'copy_keywords'        : ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                              'DATE-OBS,UT-TIME,RA,DEC,AIRMASS,' +
                              'SECPIX,TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file'    : rootpath+'/setup/magimacs.swarp',  

    # default catalog settings
    'astrometry_catalogs'  : ['URAT-1', '2MASS'], 
    'photometry_catalogs'  : ['SDSS-R9', 'URAT-1', '2MASS'],
}


# Calar Alto 1.23m, DLR-MKIII
ca123dlrmkiii_param = {
    'telescope_instrument' : 'Calar Alto 1.23m/DLR-MKIII',# telescope/instrument
    'telescope_keyword'    : 'CA123DLRMKIII',   # telescope/instrument keyword
    'observatory_code'     : '493',         # MPC observatory code
    'secpix'               : (0.3132, 0.3132 ), # pixel size (arcsec)
                                            # before binning
    'ext_coeff'            : 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx'                : True, 
    'flipy'                : True, 
    'rotate'               : 0, 

    # instrument-specific FITS header keywords
    'binning'              : ('CCDBINX', 'CCDBINY'), 
                           # binning in x/y, '_blankN' denotes that both axes
                           # are listed in one keyword, sep. by blanks
    'extent'               : ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra'                   : 'RA',  # telescope pointing, RA
    'dec'                  : 'DEC', # telescope pointin, Dec 
    'radec_separator'      : 'XXX',   # RA/Dec hms separator, use 'XXX'
                                    # if already in degrees
    'date_keyword'         : 'DATE-OBS', # obs date/time
                                         # keyword; use
                                         # 'date|time' if
                                         # separate
    'obsmidtime_jd'        : 'MIDTIMJD', # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object'               : 'OBJECT',  # object name keyword 
    'filter'               : 'FILTER',  # filter keyword
    'filter_translations'  : {'V_Johnson': 'V', 'R_Johnson': 'R', 'free': None},
                             # filtername translation dictionary
    'exptime'              : 'EXPTIME', # exposure time keyword (s)
    'airmass'              : 'AIRMASS', # airmass keyword


    # source extractor settings
    'source_minarea'       : 15, # default sextractor source minimum N_pixels
    'aprad_default'        : 4, # default aperture radius in px 
    'aprad_range'          : [2, 15], # [minimum, maximum] aperture radius (px)
    'sex-config-file'      : rootpath+'/setup/ca123dlrmkiii.sex',
    'mask_file'            : {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file'    : rootpath+'/setup/ca123dlrmkiii.scamp', 

    # swarp settings
    'copy_keywords'        : ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                              'DATE-OBS,RA,DEC,AIRMASS,' +
                              'SECPIX,TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file'    : rootpath+'/setup/ca123dlrmkiii.swarp',  

    # default catalog settings
    'astrometry_catalogs'  : ['URAT-1', '2MASS'], 
    'photometry_catalogs'  : ['SDSS-R9', 'URAT-1', '2MASS'],
}


##### access functions for telescope configurations


implemented_telescopes = ['VATT4k', 'ANDICAM', 'DCTLMI', 'ARC35ARCTIC',
                          'ARC35AGILE', 'MAGIMACS']

instrument_identifiers = {'= "Vatt4k"':        'VATT4k',
                          'ANDICAM-CCD':       'ANDICAM',
                          'LMI':               'DCTLMI',
                          'arctic':            'ARC35ARCTIC',
                          'agile':             'ARC35AGILE',
                          'IMACS Long-Camera': 'MAGIMACS',
                          'DLR-MKIII':         'CA123DLRMKIII'}




telescope_parameters = {'VATT4k' :       vatt4k_param, 
                        'ANDICAM':       andicam_param,
                        'DCTLMI':        dctlmi_param,
                        'ARC35ARCTIC':   arc35arctic_param,
                        'ARC35AGILE':    arc35agile_param,
                        'MAGIMACS':      magimacs_param,
                        'CA123DLRMKIII': ca123dlrmkiii_param}
