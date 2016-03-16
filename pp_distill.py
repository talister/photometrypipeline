#!/usr/bin/env python

""" PP_DISTILL - distill calibrated image databases into one database
                 of select moving or fixed sources
    v1.0: 2016-01-24, michael.mommert@nau.edu
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
import logging
import argparse
import time
import sqlite3
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
from scipy.optimize import minimize
import callhorizons
import scipy.ndimage.interpolation

# pipeline-specific modules
import _pp_conf
from catalog import *
from toolbox import *
import diagnostics as diag

# setup logging
logging.basicConfig(filename = _pp_conf.log_filename, 
                    level    = _pp_conf.log_level,
                    format   = _pp_conf.log_formatline, 
                    datefmt  = _pp_conf.log_datefmt)


def distill(catalogs, man_targetname, offset, display=False,
            diagnostics=False):
    """
    distill wrapper
    """

    # start logging
    logging.info('starting distill with parameters: %s' % \
                 (', '.join([('%s: %s' % (var, str(val))) for 
                             var, val in locals().items()])))

    output = {}

    ### read in database files (if necessary)
    if type(catalogs[0]) == str:
        filenames = catalogs[:]
        catalogs = []
        for filename in filenames:
            filename = filename[:filename.find('.fit')]+'.ldac.db'
            cat = catalog(filename)
            try:
                cat.read_database(filename)
            except IOError:
                logging.error('Cannot find database', filename)
                print 'Cannot find database', filename
                continue
            except sqlite3.OperationalError:
                logging.error('File %s is not a database file' % filename)
                print 'File %s is not a database file' % filename
                continue
            catalogs.append(cat)


    ##### identify target names and types

    objects = [] # one dictionary for each target

    ### select a sufficiently bright star as control star
    # match the first and the last catalog and pick a bright star
    match = catalogs[0].match_with(catalogs[-1],
            match_keys_this_catalog=['ra.deg', 'dec.deg'],
            match_keys_other_catalog=['ra.deg', 'dec.deg'],
            extract_this_catalog=['ra.deg', 'dec.deg', 'FLAGS'],
            extract_other_catalog=['ra.deg', 'dec.deg', 'FLAGS', 'MAG_APER'],
            tolerance=1./3600.)

    if len(match[0][0]) > 0:

        ctlstar_idx = numpy.argsort(match[1][3])[int(0.1*len(match[1][3]))]

        for cat_idx, cat in enumerate(catalogs):
            objects.append({'ident'      : 'control_star',
                            'obsdate.jd' :  cat.obstime,
                            'cat_idx'    :  cat_idx,
                            'ra.deg'     :  match[1][0][ctlstar_idx],
                            'dec.deg'    :  match[1][1][ctlstar_idx]})


    ### check for moving targets
    if display:
        print 'checking JPL Horizons for moving targets'

    obsparam = _pp_conf.telescope_parameters[\
                        catalogs[0].origin\
                                             .split(';')[0].strip()]

    for cat_idx, cat in enumerate(catalogs):
        targetname = cat.obj.replace('_', ' ')
        if man_targetname is not None:
            targetname = man_targetname.replace('_', ' ')
            cat.obj = targetname
        eph = callhorizons.query(targetname)
        eph.set_discreteepochs(cat.obstime[0])
        n = eph.get_ephemerides(obsparam['observatory_code'])
        if n is None or n == 0:
            logging.warning('WARNING: No position from Horizons! '+\
                            'Name (%s) correct?' % cat.obj.replace('_', ' '))
            logging.warning('HORIZONS call: %s' % eph.url)
            print ('WARNING: No position from Horizons! '+\
                   'Name (%s) correct?' % cat.obj.replace('_', ' '))

        else:
            objects.append({'ident': cat.obj,
                            'obsdate.jd': cat.obstime[0],
                            'cat_idx'   : cat_idx,
                            'ra.deg'    : eph[0]['RA']+offset[0]/3600.,
                            'dec.deg'   : eph[0]['DEC']+offset[1]/3600.})
            logging.info('Successfully grabbed Horizons position for %s ' %
                         cat.obj.replace('_', ' '))
            logging.info('HORIZONS call: %s' % eph.url)

    if display:
        print len(objects)/len(catalogs), \
            'targets per frame identified'
    logging.info('%d targets per frame identified' %
                 int(len(objects)/len(catalogs)))




    ##### extract source data for identified targets

    data = []
    targetnames = {}

    # sort objects by catalog idx
    for cat_idx, cat in enumerate(catalogs):

        objects_thiscat = filter(lambda obj:obj['cat_idx']==cat_idx,
                                 objects)

        # create a new catalog
        target_cat = catalog('targetlist:_'+cat.catalogname)
        target_cat.add_fields(['ident', 'ra.deg', 'dec.deg'],
                              [[obj['ident'] for obj in objects_thiscat],
                               [obj['ra.deg'] for obj in objects_thiscat],
                               [obj['dec.deg'] for obj in objects_thiscat]],
                              ['20A', 'D', 'D'])

        # identify magnitudes
        mag_keys = ['MAG_APER', 'MAGERR_APER']
        for key in cat.fields:
            if 'mag' in key:
                mag_keys.append(key)

        # build field lists
        match_keys_other_catalog, extract_other_catalog = [], []        
        for key in ['ra.deg', 'dec.deg', 'XWIN_IMAGE', 'YWIN_IMAGE']:
            if key in cat.fieldnames.keys():
                match_keys_other_catalog.append(cat.fieldnames[key])
                extract_other_catalog.append(cat.fieldnames[key])
            else:
                match_keys_other_catalog.append(key)
                extract_other_catalog.append(key)

        match = target_cat.match_with \
                  (cat,
                   match_keys_this_catalog=('ra.deg', 'dec.deg'),
                   match_keys_other_catalog=match_keys_other_catalog,
                   extract_this_catalog=['ra.deg', 'dec.deg', 'ident'],
                   extract_other_catalog=extract_other_catalog+mag_keys,
                   tolerance=None)


        for i in range(len(match[0][0])):
            # derive calibrated magnitudes, if available
            try:
                cal_mag = match[1][len(extract_other_catalog)+2][i] 
                cal_magerr = match[1][len(extract_other_catalog)+3][i]
            except IndexError:
                # use instrumental magnitudes
                cal_mag = match[1][len(extract_other_catalog)][i]
                cal_magerr = match[1][len(extract_other_catalog)+1][i]

            data.append([match[0][2][i], match[0][0][i], match[0][1][i],
                         match[1][0][i], match[1][1][i], 
                         match[1][len(extract_other_catalog)][i], 
                         match[1][len(extract_other_catalog)+1][i],
                         cal_mag, cal_magerr, 
                         cat.obstime, cat.catalogname,
                         match[1][2][i], match[1][3][i],
                         cat.origin])
            # format: ident, RA_exp, Dec_exp, RA_img, Dec_img,
            #         mag_inst, sigmag_instr, mag_cal, sigmag_cal
            #         obstime, filename, img_x, img_y, origin

            targetnames[match[0][2][i]] = 1

    output['targetnames'] = targetnames

    ##### write results to ASCII file

    for target in targetnames:

        output[target] = []

        if display:
            print 'write photometry results for %s' % target
        outf = open('photometry_%s.dat' % target.replace(' ', '_'), 'w')
        outf.writelines('#                filename     julian_date  ' + \
                        'ast_mag ast_sig          ast_ra      ast_dec    ' + \
                        '[1]   [2]    [3]   [4]    [5]      ZP ZP_sig  ' + \
                        'inst_mag in_sig                 [6] [7]    [8]\n')

        for dat in data:

            # sort measured magnitudes by target
            if dat[0] == target:
                try:
                    filtername = dat[13].split(';')[3]
                except IndexError:
                    filtername = '-'
                try:
                    catalogname = dat[13].split(';')[2]
                except IndexError:
                    catalogname = dat[13].split(';')[1]

                output[target].append(dat)
                outf.write(('%26.26s ' % dat[10].replace(' ', '_')) + 
                           ('%15.7f  ' % dat[9][0]) +
                           ('%7.4f '   % dat[7]) +
                           ('%6.4f '   % dat[8]) +
                           ('%13.8f '  % dat[3]) +
                           ('%+13.8f  '% dat[4]) +
                           ('%5.2f '   % ((dat[1]-dat[3])*3600.)) +
                           ('%5.2f  '  % ((dat[2]-dat[4])*3600.)) +
                           ('%5.2f '   % offset[0]) +
                           ('%5.2f  '  % offset[1]) +
                           ('%5.2f '   % dat[9][1]) +
                           ('%7.4f '   % (dat[7]-dat[5])) +
                           ('%6.4f '   % numpy.sqrt(dat[8]**2-dat[6]**2)) +
                           ('%7.4f '   % dat[5]) +
                           ('%6.4f  '  % dat[6]) +
                           ('%s  '   % catalogname) + 
                           ('%s  '   % filtername) + 
                           ('%s\n'   % dat[13].split(';')[0]))
                
        outf.writelines('#\n# [1]: Horizons_RA - image_RA [arcsec]\n'+
                        '# [2]: Horizons_DDec - image_Dec [arcsec]\n'+
                        '# [3,4]: manual target offsets in RA and DEC ' +
                        '[arcsec]\n'+
                        '# [5]: exposure time (s)\n'+
                        '# [6]: photometric catalog\n' + 
                        '# [7]: photometric band\n' +
                        '# [8]: telescope/instrument\n')
        outf.close()

           
    ### output content
    #
    # { 'targetnames': list of all targets,
    #   '(individual targetname)': [ident, RA_exp, Dec_exp, RA_img, Dec_img,
    #                               mag_inst, sigmag_instr, mag_cal, sigmag_cal
    #                               obstime, filename, img_x, img_y],
    # }
    ###

    ##### create diagnostics 
    if diagnostics:
        if display:
            print 'extracting thumbnail images'
        diag.add_results(output)


    return output



if __name__ == '__main__':

    # command line arguments    
    parser = argparse.ArgumentParser(description='distill sources of interest')
    parser.add_argument('-target', help='target name', default=None)
    parser.add_argument('-offset', help='primary target offset (arcsec)', 
                        nargs=2, default=[0,0])
    parser.add_argument('images', help='images to process', nargs='+')
    args = parser.parse_args()
    man_targetname = args.target
    man_offset = [float(coo) for coo in args.offset]
    filenames = args.images

    # check if input filenames is actually a list
    if len(filenames) == 1:
        if filenames[0].find('.lst') > -1 or filenames[0].find('.list') > -1:
            filenames = [filename[:-1] for filename in open(filenames[0], 'r').\
                         readlines()]

    distillate = distill(filenames, man_targetname, man_offset,
                         display=True, diagnostics=True)


