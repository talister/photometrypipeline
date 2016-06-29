#!/usr/bin/env python

""" PP_PREPARE - prepare fits images for photometry pipeline
    v1.0: 2016-02-27, michael.mommert@nau.edu
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


import numpy
import os
import sys
import shutil
import logging
import subprocess
import argparse, shlex
import time
from astropy.io import fits

### pipeline-specific modules
import _pp_conf
from catalog import *
import pp_extract
import diagnostics as diag
import toolbox

# setup logging
logging.basicConfig(filename = _pp_conf.log_filename, 
                    level    = _pp_conf.log_level,
                    format   = _pp_conf.log_formatline, 
                    datefmt  = _pp_conf.log_datefmt)



def prepare(filenames, obsparam, flipx=False, flipy=False, rotate=0,
            man_ra=None, man_dec=None, diagnostics=False,
            display=False):
    """
    prepare wrapper
    output: diagnostic properties
    """

    # start logging
    logging.info('preparing data with parameters: %s' % \
                 (', '.join([('%s: %s' % (var, str(val))) for 
                             var, val in locals().items()])))

    ##### prepare image headers for photometry pipeline

    for filename in filenames:

        if display:
            print 'preparing', filename

        # open image file
        hdulist = fits.open(filename, mode='update')
        header = hdulist[0].header

        # add other headers, if available
        if len(hdulist) > 1:
            for i in range(len(hdulist)):
                try:
                    header += hdulist[i].header
                except:
                    pass

        # account for flips and rotation in telescope configuration
        # if instrument has several chips...
        if 'chip_id' in obsparam:
            chip_id = header[obsparam['chip_id']]
            this_flipx = obsparam['flipx'][chip_id]
            this_flipy = obsparam['flipy'][chip_id]
            this_rotate = obsparam['rotate'][chip_id]
        # if not...
        else:
            this_flipx = obsparam['flipx']
            this_flipy = obsparam['flipy']
            this_rotate = obsparam['rotate']

        if flipx:
            this_flipx = numpy.invert(this_flipx)
        if flipy:
            this_flipy = numpy.invert(this_flipy)
        if rotate > 0:
            this_rotate += rotate

        # read image data
        imdata = hdulist[0].data

        ### add header keywords for Source Extractor
        if 'EPOCH' not in header:
            header['EPOCH'] = (2000, 'PP: required for registration')
        if 'EQUINOX' not in header:
            header['EQUINOX'] = (2000, 'PP: required for registration')

        # read out image binning mode
        if '_' in obsparam['binning'][0]:
            if '_blank' in obsparam['binning'][0]:
                binning_x = float(header[obsparam['binning'][0].\
                                         split('_')[0]].split()[0])
                binning_y = float(header[obsparam['binning'][1].\
                                         split('_')[0]].split()[1])
            if '_x' in obsparam['binning'][0]:
                binning_x = float(header[obsparam['binning'][0].\
                                         split('_')[0]].split('x')[0])
                binning_y = float(header[obsparam['binning'][1].\
                                         split('_')[0]].split('x')[1])
        else:
            binning_x = header[obsparam['binning'][0]]
            binning_y = header[obsparam['binning'][1]]

        # add pixel resolution keyword
        header['SECPIX'] = (\
        numpy.average([obsparam['secpix'][0]*binning_x,
                       obsparam['secpix'][1]*binning_y]), 
        'PP: pixel size in arcsec (after binning)')

        # # set NAXIS1 and NAXIS2 based on data array size
        # # (might be set improperly after cropping)
        # header['NAXIS1'] = (imdata.shape[1], 'PP: based on data array shape')
        # header['NAXIS2'] = (imdata.shape[0], 'PP: based on data array shape')

        # remove keywords that might collide with fake wcs
        for key in header.keys():
            if 'CD' in key and '_' in key:
                header.remove(key)
            elif 'PV' in key and '_' in key:
                header.remove(key)
            elif key in ['CTYPE1', 'CRPIX1', 'CRVAL1', 'CROTA1', 
                       'CFINT1', 'CTYPE2', 'CRPIX2', 'CRVAL2', 
                       'CFINT2',  
                       'LTM1_1', 'LTM2_2', 'WAT0_001', 'LTV1', 
                         'LTV2', 'PIXXMIT', 'PIXOFFST']:
                # removed 'CDELT1', 'CDELT2', 'CRDELT1', 'CRDELT2' from list
                # used by LOWELL31, LOWELL90

                header.remove(key)

        # add header keywords for SCAMP
        header['PHOTFLAG']  = ('F', 'PP: data is not photometric (SCAMP)')
        header['PHOT_K']  = (0.05, 'PP: assumed extinction coefficient')

        # create observation midtime jd 
        if obsparam['date_keyword'].find('|') == -1:
            header['MIDTIMJD'] = \
                (toolbox.dateobs_to_jd(header[obsparam['date_keyword']]) + \
                 float(header[obsparam['exptime']])/2./86400.,
                 'PP: obs midtime')
        else:
            datetime = header[obsparam['date_keyword'].split('|')[0]]+'T'+\
                       header[obsparam['date_keyword'].split('|')[1]]
            datetime = datetime.replace('/', '-')
            header['MIDTIMJD'] = (toolbox.dateobs_to_jd(datetime) + \
                                  float(header[obsparam['exptime']])/2./86400.,
                                  'PP: obs midtime')


        # other keywords
        header['TELINSTR'] = (obsparam['telescope_instrument'],
                              'PP: tel/instr name')
        header['TEL_KEYW'] = (obsparam['telescope_keyword'],
                              'PP: tel/instr keyword') 
        header['FILTER'] = (header[obsparam['filter']], 'PP:copied')
        header['EXPTIME'] = (header[obsparam['exptime']], 'PP: copied')
        if obsparam['airmass'] in header:
            header['AIRMASS'] = (header[obsparam['airmass']], 'PP: copied')
        else:
            header['AIRMASS'] = (1, 'PP: fake airmass')
            
        ##### add fake wcs information that is necessary to run SCAMP

        if obsparam['radec_separator'] == 'XXX':
            ra_deg  = float(header[obsparam['ra']])
            dec_deg = float(header[obsparam['dec']])
        else:
            ra_string  = header[obsparam['ra']].split(
                obsparam['radec_separator'])
            dec_string = header[obsparam['dec']].split(
                obsparam['radec_separator'])
            ra_deg = 15.*(float(ra_string[0])+float(ra_string[1])/60.+\
                          float(ra_string[2])/3600.)
            dec_deg = abs(float(dec_string[0]))+\
                      float(dec_string[1])/60.+float(dec_string[2])/3600.
            if dec_string[0].find('-') > -1:
                dec_deg = -1 * dec_deg

        if man_ra is not None and man_dec is not None:
            ra_deg = float(man_ra)
            dec_deg = float(man_dec)

        #### note to self: put this into code block below
        if obsparam['telescope_instrument'] == 'UKIRTWFCAM':
            ra_deg = float(header['TELRA'])/24.*360. - \
                     float(header['JITTER_X'])/3600.
            dec_deg = float(header['TELDEC']) - \
                      float(header['JITTER_Y'])/3600.

        # apply flips
        xnorm, ynorm = 1, 1
        if this_flipx:
            xnorm = -1
        if this_flipy:
            ynorm = -1

        # check if instrument has a chip offset
        ra_offset, dec_offset = 0, 0
        if (man_ra is None or man_dec is None) and \
           'chip_offset_fixed' in obsparam:
            cid = header[obsparam['chip_id']]
            ra_offset = float(obsparam['chip_offset_fixed'][cid][0])
            dec_offset = float(obsparam['chip_offset_fixed'][cid][1])



        ### create fake header
        #header['WCSDIM'] = (2, 'WCS Dimensionality')
        header['RADECSYS'] = ('FK5', 'PP: fake wcs coordinates')
        header['RADESYS'] = ('FK5', 'PP: fake wcs coordinates')
        header['CTYPE1'] = ('RA---TAN', 'PP: fake Coordinate type')
        header['CTYPE2'] = ('DEC--TAN', 'PP: fake Coordinate type')
        header['CRVAL1'] = (ra_deg+ra_offset, 
                            'PP: fake Coordinate reference value')
        header['CRVAL2'] = (dec_deg+dec_offset, 
                            'PP: fake Coordinate reference value')
        header['CRPIX1'] = (int(float(header[obsparam['extent'][0]])/2), 
                            'PP: fake Coordinate reference pixel')
        header['CRPIX2'] = (int(float(header[obsparam['extent'][1]])/2),
                            'PP: fake Coordinate reference pixel')

        header['CD1_1']  = (xnorm * numpy.cos(this_rotate/180.*numpy.pi) * \
                obsparam['secpix'][0]*binning_x/3600., \
                                             'PP: fake Coordinate matrix')
        header['CD1_2']  = (ynorm * -numpy.sin(this_rotate/180.*numpy.pi) * \
                obsparam['secpix'][1]*binning_y/3600., \
                                             'PP: fake Coordinate matrix')
        header['CD2_1']  = (xnorm * numpy.sin(this_rotate/180.*numpy.pi) * \
                obsparam['secpix'][0]*binning_x/3600., \
                                             'PP: fake Coordinate matrix')
        header['CD2_2']  = (ynorm * numpy.cos(this_rotate/180.*numpy.pi) * \
                obsparam['secpix'][1]*binning_y/3600., \
                                             'PP: fake Coordinate matrix')

        hdulist.flush()
        hdulist.close()


        logging.info('created fake wcs information for image %s' % filename)
    
    ##### create diagnostics 
    if diagnostics:
        diag.create_index(filenames, os.getcwd(), obsparam, display)    

    logging.info('Done! -----------------------------------------------------')

    return None




if __name__ == '__main__':
    
    # command line arguments                                                
    parser = argparse.ArgumentParser(description=('prepare data for '+\
                                                  'photometry pipeline'))
    parser.add_argument('images', help='images to process', nargs='+')
    parser.add_argument("-ra", 
                        help='image center position (RA J2000.0, deg)')
    parser.add_argument("-dec", 
                        help='image center position (Dec J2000.0, deg)')
    parser.add_argument('-flipx', help='flip fake wcs x-axis', 
                        action="store_true")
    parser.add_argument('-flipy', help='flip fake wcs y-axis', 
                        action="store_true")
    parser.add_argument('-rotate', help='rotate fake wcs by angle (deg)', 
                        default=0)
    parser.add_argument("-telescope", help='manual telescope override',
                        default=None)


    args = parser.parse_args()         
    man_ra = args.ra
    if man_ra is not None:
        man_ra = float(man_ra)
    man_dec = args.dec
    if man_dec is not None:
        man_dec = float(man_dec)
    man_flipx = args.flipx
    man_flipy = args.flipy
    man_rotate = float(args.rotate)
    telescope = args.telescope
    filenames = args.images

    ### read telescope information from fits headers
    instruments = []
    for filename in filenames:
        hdulist = fits.open(filename)
        header = hdulist[0].header
        for key in _pp_conf.instrument_keys:
            if key in header:
                instruments.append(header[key])
        hdulist.close()

    if len(instruments) == 0 and telescope is None:
        raise KeyError('cannot identify telescope/instrument; please update' + \
                       '_pp_conf.instrument_keys accordingly')
        
    if telescope is None:
        try:
            telescope = _pp_conf.instrument_identifiers[instruments[0]]
        except KeyError:
            raise KeyError('cannot identify telescope/instrument; please ' + \
                           'update _pp_conf.instrument_identifiers accordingly')

    obsparam = _pp_conf.telescope_parameters[telescope]

    # run prepare wrapper
    preparation = prepare(filenames, obsparam, flipx=man_flipx, flipy=man_flipy,
                          man_ra=man_ra, man_dec=man_dec, rotate=man_rotate,
                          diagnostics=True, display=True)


