#!/usr/bin/env python

""" PP_REGISTER - wcs register frames 
    v1.0: 2015-12-30, michael.mommert@nau.edu
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

# setup logging
logging.basicConfig(filename = _pp_conf.log_filename, 
                    level    = _pp_conf.log_level,
                    format   = _pp_conf.log_formatline, 
                    datefmt  = _pp_conf.log_datefmt)



def register(filenames, telescope, sex_snr, source_minarea, aprad,
             mancat, obsparam, display=False, diagnostics=False,
             dont_run_registration_again=False):
    """
    registration wrapper
    output: diagnostic properties
    """

    # start logging
    logging.info('starting registration with parameters: %s' % \
                 (', '.join([('%s: %s' % (var, str(val))) for 
                             var, val in locals().items()])))
    
    # check if images have been run through pp_prepare
    try:
        midtime_jd = fits.open(filenames[0])[0].header['MIDTIMJD']
    except KeyError:
        raise KeyError(('%s image header incomplete, have the data run ' + 
                        'through pp_prepare?') % filenames[0])
        return None

    ##### run extract routines
    if display:
        print '* extract sources from %d frames' % len(filenames)

    extractparameters = {'sex_snr':sex_snr,
                         'source_minarea':source_minarea, \
                         'aprad':aprad, 'telescope':telescope, \
                         'quiet':False}
    extraction = pp_extract.extract_multiframe(filenames, extractparameters)

    if extraction is None:
        if display:
            print 'ERROR: extraction was not successful'
        logging.error('extraction was not successful')
        return None
    
    ldac_files = [filename[:filename.find('.fit')]+'.ldac' \
                  for filename in filenames]

    ##### run scamp on all image catalogs using different catalogs
    if mancat is not None:
        obsparam['astrometry_catalogs'] = [mancat]

    fileline = " ".join(ldac_files)

    for cat_idx, refcat in enumerate(obsparam['astrometry_catalogs']):

        output = {}

        ### check if sufficient reference stars are available in refcat
        logging.info('check if sufficient reference stars in catalog %s' %
                     refcat)
        hdulist = fits.open(filenames[len(filenames)/2])
        ra = float(hdulist[0].header['CRVAL1'])
        dec = float(hdulist[0].header['CRVAL2'])
        rad = max([float(hdulist[0].header[obsparam['extent'][0]])*
                   obsparam['secpix'][0],
                   float(hdulist[0].header[obsparam['extent'][1]])*
                   obsparam['secpix'][1]])/3600.
        checkrefcat = catalog(refcat)
        n_sources = checkrefcat.download_catalog(ra, dec, rad, 100, 
                                                 save_catalog=False)
        if n_sources < _pp_conf.min_sources_astrometric_catalog:
            logging.info(('Only %d sources in astrometric reference catalog; ' \
                          + 'try other catalog') % n_sources)
            goodfits, badfits = [], filenames
            continue
        else:
            logging.info('%d sources in catalog %s; enough for SCAMP' %
                         (n_sources, refcat))

        # remove existing scamp output
        os.remove('scamp_output.xml') if os.path.exists('scamp_output.xml') \
            else None

        logging.info('run SCAMP on %d image files, match with catalog %s ' % 
                     (len(filenames), refcat))

        # assemble arguments for scamp, run it, and wait for it
        commandline = 'scamp -c '+obsparam['scamp-config-file']+ \
                      ' -ASTREF_CATALOG '+refcat+' '+fileline
        scamp = subprocess.Popen(shlex.split(commandline))
        scamp.wait()

        ##### identify successful and failed WCS registrations based on
        ##### the contrast values provided by SCAMP
        scamp = _pp_conf.read_scamp_output()
        os.rename('scamp_output.xml', 'astrometry_scamp.xml')
        goodfits, badfits = [], []
        fitresults = [] # store scamp outputs 
        for dat in scamp[1]:
            # successful fit
            if ((float(dat[scamp[0]['AS_Contrast']]) < 
                 _pp_conf.scamp_as_contrast_limit)  
                or (float(dat[scamp[0]['XY_Contrast']]) < 
                    _pp_conf.scamp_xy_contrast_limit) 
                or len(dat) == 0):
                filename = dat[scamp[0]['Catalog_Name']]
                for file in os.listdir('.'):
                    if file.find(filename[:filename.find('.ldac')]+'.fit') \
                       > -1:
                        filename = file
                badfits.append(filename)
            # failed fit
            else:
                filename = dat[scamp[0]['Catalog_Name']]
                for file in os.listdir('.'):
                    if file.find(filename[:filename.find('.ldac')]+'.fit')\
                       > -1:
                        filename = file
                goodfits.append(filename)

            fitresults.append(\
                [filename,
                 float(dat[scamp[0]['AS_Contrast']]),
                 float(dat[scamp[0]['XY_Contrast']]),               
                 float(dat[scamp[0]['AstromSigma_Reference']].split()[0]),
                 float(dat[scamp[0]['AstromSigma_Reference']].split()[1]),
                 float(dat[scamp[0]['Chi2_Reference']]), \
                 float(dat[scamp[0]['Chi2_Internal']])])
            
        open('registration_succeeded.lst', 'w').writelines("%s\n" %
                '\n'.join(goodfits))
        if len(goodfits)==0 and len(badfits)==0:
            badfits = filenames
        open('registration_failed.lst', 'w').writelines("%s\n" %
                '\n'.join(badfits))

        output['goodfits']    = goodfits
        output['badfits']     = badfits
        output['fitresults']  = fitresults
        output['catalog']     = refcat


        ##### check registration outcome

        logging.info(' > match succeeded for %d/%d images' % \
                     (len(goodfits), len(filenames)))
        print '\n################################# ' + \
            'REGISTRATION SUMMARY:\n###'
        print '### %d/%d images have been registered successfully' % \
            (len(goodfits), len(filenames))
        print '###\n###############################' + \
            '#######################\n'

        # registration succeeded for most images
        if len(goodfits) > len(badfits):
            break
        # registration failed for most (or all) images
        else:
            logging.info(' > match failed for %d/%d images' % \
                         (len(badfits), len(filenames)))

            ### in case the registration failed, try again!
            # this will make use of the .head files and improves results
            if (cat_idx == len(obsparam['astrometry_catalogs'])-1):
                if not dont_run_registration_again:
                    logging.critical('No match possible with either catalog ' \
                                     + '- try again') 
                    output = register(filenames, telescope, sex_snr,
                                      source_minarea, aprad,
                                      None, obsparam,
                                      display=True, diagnostics=False,
                                      dont_run_registration_again=True)
                    return output
                else:
                    logging.critical('No match possible with either catalog ' \
                                     + '- abort!') 
                    return output
            else:
                logging.info('No match possible with this catalog ' \
                                 + '- try next one in the list') 

                # try next catalog in the list
                continue

    ##### update image headers with wcs solutions where registration
    ##### was successful
    logging.info('update image headers with WCS solutions ') 
    for filename in goodfits:
        # remove fake wcs header keys
        fake_wcs_keys =  ['RADECSYS', 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 
                          'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 
                          'CD2_2']
        hdu = fits.open(filename, mode='update')
        for fake_key in fake_wcs_keys:
            hdu[0].header[fake_key] = ''

        # read new header files
        newhead = open(filename[:filename.find('.fit')]+'.head','r').readlines()

        for line in newhead:
            key = line[:8].strip()
            try:
                value = float(line[10:30].replace('\'', ' ').strip())
            except ValueError:
                value = line[10:30].replace('\'', ' ').strip()
            comment = line[30:].strip()
            if key.find('END') > -1:
                break
            #print key, '|',  value, '|',  comment
            hdu[0].header[key] = (str(value), comment)

        # other header keywords
        hdu[0].header['RADECSYS'] = (hdu[0].header['RADESYS'], 
                                     'copied from RADESYS')
        hdu[0].header['TEL_KEYW'] = (telescope, 'pipeline telescope keyword')
        hdu[0].header['REGCAT'] = (refcat, 'catalog used in WCS registration')
        hdu.flush()
        hdu.close()
    
        # cleaning up (in case the registration succeeded)
        os.remove(filename[:filename.find('.fit')]+'.head')

        
    if len(badfits) == len(filenames):
        if display:
            print 'ERROR: registration failed for all images'
        logging.error('ERROR: registration failed for all images')
        return output
        
    ##### print astrometry output file
    outf = open('best_astrometry.dat', 'w')
    outf.writelines('# filename AS_contrast XY_contrast ' \
                    + 'Chi2_catalog Chi2_int Pos_uncertainty(arcsec)\n')
    for idx, data in enumerate(fitresults):
        outf.writelines('%25.25s %5.2f %5.2f %10.7f %10.7f %7.4f\n' % \
                        (data[0], data[1], data[2], data[5], data[6],
                         numpy.sqrt(data[3]**2+data[4]**2)))
    outf.close()

       
    ### extraction output
    # 
    # -> see pp_extract.py
    #
    ##

    ### output content
    #
    # { 'good_fits'    : list of fits where registration succeeded,
    #   'bad fits'     : list of fits where registration failed,
    #   'fitresults'   : scamp fit results,
    #   'catalog'      : astrometric reference catalog name
    # }
    ###

    ##### create diagnostics 
    if diagnostics:
        diag.add_registration(output, extraction)

    logging.info('Done! -----------------------------------------------------')

    return output




if __name__ == '__main__':
    
    # command line arguments                                                
    parser = argparse.ArgumentParser(description='automated WCS registration')
    parser.add_argument("-snr", help='sextractor SNR threshold', default=3)
    parser.add_argument("-minarea", help='sextractor SNR threshold',
                        default=0)
    parser.add_argument("-aprad", help='aperture radius for photometry (px)', 
                        default=0)
    parser.add_argument("-cat", help='manually select reference catalog', 
                        default=None) 
    parser.add_argument('images', help='images to process', nargs='+')

    args = parser.parse_args()         
    sex_snr = float(args.snr)
    source_minarea = float(args.minarea)
    aprad = float(args.aprad)
    mancat = args.cat
    filenames = args.images


    ### read telescope and filter information from fits headers
    # check that they are the same for all images
    instruments = []
    for filename in filenames:
        hdulist = fits.open(filename)
        header = hdulist[0].header
        for key in _pp_conf.instrument_keys:
            if key in header:
                instruments.append(header[key])

    if len(instruments) == 0:
        raise KeyError('cannot identify telescope/instrument; please update' + \
                       '_pp_conf.instrument_keys accordingly')


    # assign telescope parameters (telescopes.py)
    telescope = _pp_conf.instrument_identifiers[instruments[0]]
    obsparam = _pp_conf.telescope_parameters[telescope]

    # set aperture photometry aperture radius
    if aprad == 0:
        aprad = obsparam['aprad_default']

    # set minarea from obsparam
    if source_minarea == 0:
        source_minarea = obsparam['source_minarea']

    # run registration wrapper
    registration = register(filenames, telescope, sex_snr,
                            source_minarea, aprad, mancat, obsparam,
                            display=True, diagnostics=True)


