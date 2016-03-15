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




##### access functions for telescope configurations


implemented_telescopes = ['VATT4k', 'ANDICAM', 'DCTLMI']

instrument_identifiers = {'= "Vatt4k"': 'VATT4k',
                          'ANDICAM-CCD':'ANDICAM', 'LMI': 'DCTLMI'}




telescope_parameters = {'VATT4k' : vatt4k_param, 
                        'ANDICAM': andicam_param,
                        'DCTLMI': dctlmi_param}
