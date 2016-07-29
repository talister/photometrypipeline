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


def manual_positions(posfile, catalogs):
    """create targets for manually provided positions (using -positions
    option)"""

    print 'target positions manually provided'
    logging.info('target positions manually provided')

    positions = numpy.genfromtxt(posfile, dtype=[('filename', 'S50'), 
                                                 ('MJD', float), 
                                                 ('ra', float), 
                                                 ('dec', float)])
    try:
        assert len(positions) == len(catalogs)
    except AssertionError:
        print (posfile + ' is not complete; has to provide a position ' +
               'for each frame')
        logging.error(posfile+' is not complete; has to provide position ' +
                      'for each frame')
        return []

    objects = []
    for cat_idx, cat in enumerate(catalogs):
        objects.append({'ident'      : 'manual_target',
                        'obsdate.jd' :  cat.obstime,
                        'cat_idx'    :  cat_idx,
                        'ra.deg'     :  positions[cat_idx]['ra'],
                        'dec.deg'    :  positions[cat_idx]['dec']})
        
    return objects


def pick_controlstar(catalogs):
    """match the first and the last catalog and pick a bright star"""

    print 'pick control star'

    match = catalogs[0].match_with(catalogs[-1],
            match_keys_this_catalog=['ra.deg', 'dec.deg'],
            match_keys_other_catalog=['ra.deg', 'dec.deg'],
            extract_this_catalog=['ra.deg', 'dec.deg', 'FLAGS'],
            extract_other_catalog=['ra.deg', 'dec.deg', 'FLAGS', 'MAG_APER'],
            tolerance=1./3600.)

    objects = []
    if len(match[0][0]) > 0:

        ctlstar_idx = numpy.argsort(match[1][3])[int(0.1*len(match[1][3]))]

        for cat_idx, cat in enumerate(catalogs):
            objects.append({'ident'      : 'control_star',
                            'obsdate.jd' :  cat.obstime[0],
                            'cat_idx'    :  cat_idx,
                            'ra.deg'     :  match[1][0][ctlstar_idx],
                            'dec.deg'    :  match[1][1][ctlstar_idx]})
    else:
        print '  no common control star found in first and last frame'
        logging.info('no common control star found in first and last frame')

    return objects


def moving_primary_target(catalogs, man_targetname, offset):

    print 'check JPL Horizons for primary target'

    obsparam = _pp_conf.telescope_parameters[
                        catalogs[0].origin.split(';')[0].strip()]

    objects = []
    for cat_idx, cat in enumerate(catalogs):
        targetname = cat.obj.replace('_', ' ')
        if man_targetname is not None:
            targetname = man_targetname.replace('_', ' ')
            cat.obj = targetname
        eph = callhorizons.query(targetname)
        eph.set_discreteepochs(cat.obstime[0])

        try:
            n = eph.get_ephemerides(obsparam['observatory_code'])
        except ValueError:
            print 'Target (%s) is not an asteroid' % targetname
            logging.warning('Target (%s) is not an asteroid' % targetname)
            n = None
            
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

    return objects


def fixed_target(fixed_coo, catalogs):
    """add fixed target position to object catalog (if different from [0,0])"""

    objects = []
    for cat_idx, cat in enumerate(catalogs):
        objects.append({'ident': 'fixed_target',
                        'obsdate.jd': cat.obstime[0],
                        'cat_idx'   : cat_idx,
                        'ra.deg'    : fixed_coo[0],
                        'dec.deg'   : fixed_coo[1]})
    return objects


def serendipitous_asteroids():
    return []
#     return [{'obsdate.jd': 2457424.974246285, 'cat_idx': 0, 'ident': '2000 SK176', 'dec.deg': -8.41924, 'ra.deg': 209.41220999999996},
# {'obsdate.jd': 2457424.975281852, 'cat_idx': 1, 'ident': '2000 SK176', 'dec.deg': -8.419283831045874, 'ra.deg': 209.41246562833572},
# {'obsdate.jd': 2457424.976279572, 'cat_idx': 2, 'ident': '2000 SK176', 'dec.deg': -8.419326069771497, 'ra.deg': 209.4127165691116},
# {'obsdate.jd': 2457424.97727735, 'cat_idx': 3, 'ident': '2000 SK176', 'dec.deg': -8.419368314306624, 'ra.deg': 209.41297010875692},
# {'obsdate.jd': 2457424.978277558, 'cat_idx': 4, 'ident': '2000 SK176', 'dec.deg': -8.419410659122125, 'ra.deg': 209.41322489138415},
# {'obsdate.jd': 2457424.984377454, 'cat_idx': 5, 'ident': '2000 SK176', 'dec.deg': -8.419668912380356, 'ra.deg': 209.41476360944603},
# {'obsdate.jd': 2457424.98541316, 'cat_idx': 6, 'ident': '2000 SK176', 'dec.deg': -8.419712545266849, 'ra.deg': 209.41502378342423},
# {'obsdate.jd': 2457424.986410602, 'cat_idx': 7, 'ident': '2000 SK176', 'dec.deg': -8.41975465873114, 'ra.deg': 209.41527250678584},
# {'obsdate.jd': 2457424.987408287, 'cat_idx': 8, 'ident': '2000 SK176', 'dec.deg': -8.419797297705333, 'ra.deg': 209.4155227577109},
# {'obsdate.jd': 2457424.988405325, 'cat_idx': 9, 'ident': '2000 SK176', 'dec.deg': -8.4198403180754, 'ra.deg': 209.41577592103883},
# {'obsdate.jd': 2457424.994686216, 'cat_idx': 10, 'ident': '2000 SK176', 'dec.deg': -8.420094996365746, 'ra.deg': 209.41735303572193},
# {'obsdate.jd': 2457424.995725868, 'cat_idx': 11, 'ident': '2000 SK176', 'dec.deg': -8.420139545128215, 'ra.deg': 209.41761250102138},
# {'obsdate.jd': 2457424.996723149, 'cat_idx': 12, 'ident': '2000 SK176', 'dec.deg': -8.420181698807244, 'ra.deg': 209.41786144592314},
# {'obsdate.jd': 2457424.997719352, 'cat_idx': 13, 'ident': '2000 SK176', 'dec.deg': -8.420223650267445, 'ra.deg': 209.41811187690288},
# {'obsdate.jd': 2457424.998718287, 'cat_idx': 14, 'ident': '2000 SK176', 'dec.deg': -8.420265867517708, 'ra.deg': 209.41836555158358},
# {'obsdate.jd': 2457425.006592801, 'cat_idx': 15, 'ident': '2000 SK176', 'dec.deg': -8.420599724452769, 'ra.deg': 209.42034096986737},
# {'obsdate.jd': 2457425.007630706, 'cat_idx': 16, 'ident': '2000 SK176', 'dec.deg': -8.420643565884676, 'ra.deg': 209.42060063257608},
# {'obsdate.jd': 2457425.008629688, 'cat_idx': 17, 'ident': '2000 SK176', 'dec.deg': -8.420684556500422, 'ra.deg': 209.42085172874533},
# {'obsdate.jd': 2457425.009628009, 'cat_idx': 18, 'ident': '2000 SK176', 'dec.deg': -8.420723610203146, 'ra.deg': 209.42110461988887},
# {'obsdate.jd': 2457425.010625266, 'cat_idx': 19, 'ident': '2000 SK176', 'dec.deg': -8.42076, 'ra.deg': 209.42136},
# {'obsdate.jd': 2457424.974246285, 'cat_idx': 0, 'ident': '2001 FV4', 'dec.deg': -8.427449999999999, 'ra.deg': 209.31837999999993},
# {'obsdate.jd': 2457424.975281852, 'cat_idx': 1, 'ident': '2001 FV4', 'dec.deg': -8.427489852516802, 'ra.deg': 209.3184364918058},
# {'obsdate.jd': 2457424.976279572, 'cat_idx': 2, 'ident': '2001 FV4', 'dec.deg': -8.427528248545878, 'ra.deg': 209.3184953944282},
# {'obsdate.jd': 2457424.97727735, 'cat_idx': 3, 'ident': '2001 FV4', 'dec.deg': -8.427566646792878, 'ra.deg': 209.31855682164425},
# {'obsdate.jd': 2457424.978277558, 'cat_idx': 4, 'ident': '2001 FV4', 'dec.deg': -8.427605138547055, 'ra.deg': 209.31861907263678},
# {'obsdate.jd': 2457424.984377454, 'cat_idx': 5, 'ident': '2001 FV4', 'dec.deg': -8.427839885495132, 'ra.deg': 209.31898346314597},
# {'obsdate.jd': 2457424.98541316, 'cat_idx': 6, 'ident': '2001 FV4', 'dec.deg': -8.42787974335857, 'ra.deg': 209.31904856291277},
# {'obsdate.jd': 2457424.986410602, 'cat_idx': 7, 'ident': '2001 FV4', 'dec.deg': -8.427918128665198, 'ra.deg': 209.31911032839363},
# {'obsdate.jd': 2457424.987408287, 'cat_idx': 8, 'ident': '2001 FV4', 'dec.deg': -8.427956523329325, 'ra.deg': 209.31916991789},
# {'obsdate.jd': 2457424.988405325, 'cat_idx': 9, 'ident': '2001 FV4', 'dec.deg': -8.427994893107298, 'ra.deg': 209.3192270918976},
# {'obsdate.jd': 2457424.994686216, 'cat_idx': 10, 'ident': '2001 FV4', 'dec.deg': -8.42823660541511, 'ra.deg': 209.3195995407443},
# {'obsdate.jd': 2457424.995725868, 'cat_idx': 11, 'ident': '2001 FV4', 'dec.deg': -8.428276615126734, 'ra.deg': 209.31966192725892},
# {'obsdate.jd': 2457424.996723149, 'cat_idx': 12, 'ident': '2001 FV4', 'dec.deg': -8.428314994259686, 'ra.deg': 209.3197239494473},
# {'obsdate.jd': 2457424.997719352, 'cat_idx': 13, 'ident': '2001 FV4', 'dec.deg': -8.42835333190346, 'ra.deg': 209.31978522832566},
# {'obsdate.jd': 2457424.998718287, 'cat_idx': 14, 'ident': '2001 FV4', 'dec.deg': -8.428391774664144, 'ra.deg': 209.31984405362797},
# {'obsdate.jd': 2457425.006592801, 'cat_idx': 15, 'ident': '2001 FV4', 'dec.deg': -8.42869481557675, 'ra.deg': 209.32030807478094},
# {'obsdate.jd': 2457425.007630706, 'cat_idx': 16, 'ident': '2001 FV4', 'dec.deg': -8.428734758058589, 'ra.deg': 209.32036730110963},
# {'obsdate.jd': 2457425.008629688, 'cat_idx': 17, 'ident': '2001 FV4', 'dec.deg': -8.428773202647825, 'ra.deg': 209.32042503892413},
# {'obsdate.jd': 2457425.009628009, 'cat_idx': 18, 'ident': '2001 FV4', 'dec.deg': -8.428811621794425, 'ra.deg': 209.32048533294335},
# {'obsdate.jd': 2457425.010625266, 'cat_idx': 19, 'ident': '2001 FV4', 'dec.deg': -8.42885, 'ra.deg': 209.32055},
# {'obsdate.jd': 2457424.974246285, 'cat_idx': 0, 'ident': '2012 TW236', 'dec.deg': -8.47544, 'ra.deg': 209.33318999999997},
# {'obsdate.jd': 2457424.975281852, 'cat_idx': 1, 'ident': '2012 TW236', 'dec.deg': -8.475413056862179, 'ra.deg': 209.33329231184868},
# {'obsdate.jd': 2457424.976279572, 'cat_idx': 2, 'ident': '2012 TW236', 'dec.deg': -8.475391210319247, 'ra.deg': 209.3333927419245},
# {'obsdate.jd': 2457424.97727735, 'cat_idx': 3, 'ident': '2012 TW236', 'dec.deg': -8.475371755719252, 'ra.deg': 209.33349383176127},
# {'obsdate.jd': 2457424.978277558, 'cat_idx': 4, 'ident': '2012 TW236', 'dec.deg': -8.47535302643146, 'ra.deg': 209.33359466448238},
# {'obsdate.jd': 2457424.984377454, 'cat_idx': 5, 'ident': '2012 TW236', 'dec.deg': -8.475215985308424, 'ra.deg': 209.33419351172768},
# {'obsdate.jd': 2457424.98541316, 'cat_idx': 6, 'ident': '2012 TW236', 'dec.deg': -8.475192389733662, 'ra.deg': 209.33429759264206},
# {'obsdate.jd': 2457424.986410602, 'cat_idx': 7, 'ident': '2012 TW236', 'dec.deg': -8.475169447542871, 'ra.deg': 209.33439729547135},
# {'obsdate.jd': 2457424.987408287, 'cat_idx': 8, 'ident': '2012 TW236', 'dec.deg': -8.475145958070772, 'ra.deg': 209.3344969435359},
# {'obsdate.jd': 2457424.988405325, 'cat_idx': 9, 'ident': '2012 TW236', 'dec.deg': -8.475122112403065, 'ra.deg': 209.3345967154096},
# {'obsdate.jd': 2457424.994686216, 'cat_idx': 10, 'ident': '2012 TW236', 'dec.deg': -8.474988318122566, 'ra.deg': 209.33522583895862},
# {'obsdate.jd': 2457424.995725868, 'cat_idx': 11, 'ident': '2012 TW236', 'dec.deg': -8.474963723744242, 'ra.deg': 209.3353268641253},
# {'obsdate.jd': 2457424.996723149, 'cat_idx': 12, 'ident': '2012 TW236', 'dec.deg': -8.474940725539769, 'ra.deg': 209.33542221934076},
# {'obsdate.jd': 2457424.997719352, 'cat_idx': 13, 'ident': '2012 TW236', 'dec.deg': -8.474918030831947, 'ra.deg': 209.3355184673156},
# {'obsdate.jd': 2457424.998718287, 'cat_idx': 14, 'ident': '2012 TW236', 'dec.deg': -8.474895244484607, 'ra.deg': 209.33561759743387},
# {'obsdate.jd': 2457425.006592801, 'cat_idx': 15, 'ident': '2012 TW236', 'dec.deg': -8.474724412766317, 'ra.deg': 209.33640266688045},
# {'obsdate.jd': 2457425.007630706, 'cat_idx': 16, 'ident': '2012 TW236', 'dec.deg': -8.474699650275673, 'ra.deg': 209.33650197247772},
# {'obsdate.jd': 2457425.008629688, 'cat_idx': 17, 'ident': '2012 TW236', 'dec.deg': -8.474675333578814, 'ra.deg': 209.33659836940507},
# {'obsdate.jd': 2457425.009628009, 'cat_idx': 18, 'ident': '2012 TW236', 'dec.deg': -8.474651730121483, 'ra.deg': 209.33669717730356},
# {'obsdate.jd': 2457425.010625266, 'cat_idx': 19, 'ident': '2012 TW236', 'dec.deg': -8.47463, 'ra.deg': 209.3368},
# {'obsdate.jd': 2457424.974246285, 'cat_idx': 0, 'ident': '2015 AU206', 'dec.deg': -8.378939999999998, 'ra.deg': 209.43675},
# {'obsdate.jd': 2457424.975281852, 'cat_idx': 1, 'ident': '2015 AU206', 'dec.deg': -8.37892462942251, 'ra.deg': 209.43681705112343},
# {'obsdate.jd': 2457424.976279572, 'cat_idx': 2, 'ident': '2015 AU206', 'dec.deg': -8.378914478028506, 'ra.deg': 209.43687717613273},
# {'obsdate.jd': 2457424.97727735, 'cat_idx': 3, 'ident': '2015 AU206', 'dec.deg': -8.378906911323863, 'ra.deg': 209.43693478336007},
# {'obsdate.jd': 2457424.978277558, 'cat_idx': 4, 'ident': '2015 AU206', 'dec.deg': -8.378899951062525, 'ra.deg': 209.436991856732},
# {'obsdate.jd': 2457424.984377454, 'cat_idx': 5, 'ident': '2015 AU206', 'dec.deg': -8.378842398337783, 'ra.deg': 209.43735518116355},
# {'obsdate.jd': 2457424.98541316, 'cat_idx': 6, 'ident': '2015 AU206', 'dec.deg': -8.37883148244436, 'ra.deg': 209.43741364476122},
# {'obsdate.jd': 2457424.986410602, 'cat_idx': 7, 'ident': '2015 AU206', 'dec.deg': -8.378819156185022, 'ra.deg': 209.43747087581932},
# {'obsdate.jd': 2457424.987408287, 'cat_idx': 8, 'ident': '2015 AU206', 'dec.deg': -8.3788084317302, 'ra.deg': 209.43753030212176},
# {'obsdate.jd': 2457424.988405325, 'cat_idx': 9, 'ident': '2015 AU206', 'dec.deg': -8.37880089829109, 'ra.deg': 209.43759205903166},
# {'obsdate.jd': 2457424.994686216, 'cat_idx': 10, 'ident': '2015 AU206', 'dec.deg': -8.378743086959227, 'ra.deg': 209.43796876124324},
# {'obsdate.jd': 2457424.995725868, 'cat_idx': 11, 'ident': '2015 AU206', 'dec.deg': -8.378732525914641, 'ra.deg': 209.43803114320397},
# {'obsdate.jd': 2457424.996723149, 'cat_idx': 12, 'ident': '2015 AU206', 'dec.deg': -8.378720261228654, 'ra.deg': 209.43808858779698},
# {'obsdate.jd': 2457424.997719352, 'cat_idx': 13, 'ident': '2015 AU206', 'dec.deg': -8.378709162394792, 'ra.deg': 209.4381448811037},
# {'obsdate.jd': 2457424.998718287, 'cat_idx': 14, 'ident': '2015 AU206', 'dec.deg': -8.37870113500953, 'ra.deg': 209.43820147140488},
# {'obsdate.jd': 2457425.006592801, 'cat_idx': 15, 'ident': '2015 AU206', 'dec.deg': -8.378620534896914, 'ra.deg': 209.43866688117768},
# {'obsdate.jd': 2457425.007630706, 'cat_idx': 16, 'ident': '2015 AU206', 'dec.deg': -8.378613019915377, 'ra.deg': 209.43872700434773},
# {'obsdate.jd': 2457425.008629688, 'cat_idx': 17, 'ident': '2015 AU206', 'dec.deg': -8.378604924638601, 'ra.deg': 209.43878499970927},
# {'obsdate.jd': 2457425.009628009, 'cat_idx': 18, 'ident': '2015 AU206', 'dec.deg': -8.378594427386407, 'ra.deg': 209.4388427734835},
# {'obsdate.jd': 2457425.010625266, 'cat_idx': 19, 'ident': '2015 AU206', 'dec.deg': -8.37858, 'ra.deg': 209.4389},
# {'obsdate.jd': [2457424.974246285,1], 'cat_idx': 0, 'ident': '2016 GW108', 'dec.deg': -8.540019999999998, 'ra.deg': 209.372},
# {'obsdate.jd': [2457424.975281852,1], 'cat_idx': 1, 'ident': '2016 GW108', 'dec.deg': -8.540040496851558, 'ra.deg': 209.37210431487213},
# {'obsdate.jd': [2457424.976279572,1], 'cat_idx': 2, 'ident': '2016 GW108', 'dec.deg': -8.540064901808284, 'ra.deg': 209.37220929270816},
# {'obsdate.jd': [2457424.97727735,1], 'cat_idx': 3, 'ident': '2016 GW108', 'dec.deg': -8.540091893368333, 'ra.deg': 209.37231679777622},
# {'obsdate.jd': [2457424.978277558,1], 'cat_idx': 4, 'ident': '2016 GW108', 'dec.deg': -8.540119575590998, 'ra.deg': 209.37242523884714},
# {'obsdate.jd': [2457424.984377454,1], 'cat_idx': 5, 'ident': '2016 GW108', 'dec.deg': -8.540273294348596, 'ra.deg': 209.373071325477},
# {'obsdate.jd': [2457424.98541316,1], 'cat_idx': 6, 'ident': '2016 GW108', 'dec.deg': -8.540298255671674, 'ra.deg': 209.37318425612656},
# {'obsdate.jd': [2457424.986410602,1], 'cat_idx': 7, 'ident': '2016 GW108', 'dec.deg': -8.540320478879721, 'ra.deg': 209.37329208473298},
# {'obsdate.jd': [2457424.987408287,1], 'cat_idx': 8, 'ident': '2016 GW108', 'dec.deg': -8.540344299750606, 'ra.deg': 209.37339774504747},
# {'obsdate.jd': [2457424.988405325,1], 'cat_idx': 9, 'ident': '2016 GW108', 'dec.deg': -8.540371279271902, 'ra.deg': 209.3735009572041},
# {'obsdate.jd': [2457424.994686216,1], 'cat_idx': 10, 'ident': '2016 GW108', 'dec.deg': -8.540530806554408, 'ra.deg': 209.37416340383075},
# {'obsdate.jd': [2457424.995725868,1], 'cat_idx': 11, 'ident': '2016 GW108', 'dec.deg': -8.540557205848053, 'ra.deg': 209.37427406985654},
# {'obsdate.jd': [2457424.996723149,1], 'cat_idx': 12, 'ident': '2016 GW108', 'dec.deg': -8.54058011577894, 'ra.deg': 209.37438232529405},
# {'obsdate.jd': [2457424.997719352,1], 'cat_idx': 13, 'ident': '2016 GW108', 'dec.deg': -8.540601878363882, 'ra.deg': 209.37448914712667},
# {'obsdate.jd': [2457424.998718287,1], 'cat_idx': 14, 'ident': '2016 GW108', 'dec.deg': -8.540623817985683, 'ra.deg': 209.3745930051677},
# {'obsdate.jd': [2457425.006592801,1], 'cat_idx': 15, 'ident': '2016 GW108', 'dec.deg': -8.540820425917447, 'ra.deg': 209.3754245044099},
# {'obsdate.jd': [2457425.007630706,1], 'cat_idx': 16, 'ident': '2016 GW108', 'dec.deg': -8.540849130321682, 'ra.deg': 209.37553710897066},
# {'obsdate.jd': [2457425.008629688,1], 'cat_idx': 17, 'ident': '2016 GW108', 'dec.deg': -8.540876060425054, 'ra.deg': 209.3756448055592},
# {'obsdate.jd': [2457425.009628009,1], 'cat_idx': 18, 'ident': '2016 GW108', 'dec.deg': -8.54090032819322, 'ra.deg': 209.3757497680407},
# {'obsdate.jd': [2457425.010625266,1], 'cat_idx': 19, 'ident': '2016 GW108', 'dec.deg': -8.54092, 'ra.deg': 209.37585}]

### -------------------


def distill(catalogs, man_targetname, offset, fixed_coo, posfile,
            display=False, diagnostics=False):

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

    ### check for positions file
    if posfile is not None:
        objects += manual_positions(posfile, catalogs)

    ### select a sufficiently bright star as control star
    objects += pick_controlstar(catalogs)

    ### check Horizons for primary target (if a moving target)
    objects += moving_primary_target(catalogs, man_targetname, offset)

    ### add fixed target
    if fixed_coo[0] != 0 and fixed_coo[1] != 0.0:
        objects += fixed_target(fixed_coo, catalogs)

    ### seredipitous asteroids
    objects += serendipitous_asteroids()


    if display:
        print len(objects)/len(catalogs), \
            'potential targets per frame identified:', \
            ", ".join(set([obj['ident'] for obj in objects]))
    logging.info('%d potential targets per frame identified: %s' %
                 (int(len(objects)/len(catalogs)), 
                  ", ".join(set([obj['ident'] for obj in objects]))))



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


        # build field lists for observed catalogs
        match_keys_other_catalog, extract_other_catalog = [], []
        for key in ['ra.deg', 'dec.deg', 'XWIN_IMAGE', 'YWIN_IMAGE', 'FLAGS']:
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
                         cat.origin, match[1][4][i]])
            # format: ident, RA_exp, Dec_exp, RA_img, Dec_img,
            #         mag_inst, sigmag_instr, mag_cal, sigmag_cal
            #         obstime, filename, img_x, img_y, origin, flags

            targetnames[match[0][2][i]] = 1


    output['targetnames'] = targetnames

    ##### write results to ASCII file

    for target in targetnames:

        output[target] = []

        if display:
            print 'write photometry results for %s' % target
        outf = open('photometry_%s.dat' % target.replace(' ', '_'), 'w')
        outf.writelines('#                          filename     julian_date' +
                        'ast_mag ast_sig        ast_ra       ast_dec    ' +
                        '[1]   [2]    [3]   [4]    [5]       ZP ZP_sig ' +
                        'inst_mag in_sig   [6] [7] [8] [9]\n')

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
                outf.write(('%35.35s ' % dat[10].replace(' ', '_')) + 
                           ('%15.7f  ' % dat[9][0]) +
                           ('%8.4f '   % dat[7]) +
                           ('%6.4f '   % dat[8]) +
                           ('%13.8f '  % dat[3]) +
                           ('%+13.8f  '% dat[4]) +
                           ('%5.2f '   % ((dat[1]-dat[3])*3600.)) +
                           ('%5.2f  '  % ((dat[2]-dat[4])*3600.)) +
                           ('%5.2f '   % offset[0]) +
                           ('%5.2f  '  % offset[1]) +
                           ('%5.2f '   % dat[9][1]) +
                           ('%8.4f '   % (dat[7]-dat[5])) +
                           ('%6.4f '   % numpy.sqrt(dat[8]**2-dat[6]**2)) +
                           ('%8.4f '   % dat[5]) +
                           ('%6.4f  '  % dat[6]) +
                           ('%s  '   % catalogname) + 
                           ('%s  '   % filtername) +
                           ('%3d  '  % dat[14]) + 
                           ('%s\n'   % dat[13].split(';')[0]))
                
        outf.writelines('#\n# [1]: Horizons_RA - image_RA [arcsec]\n'+
                        '# [2]: Horizons_DDec - image_Dec [arcsec]\n'+
                        '# [3,4]: manual target offsets in RA and DEC ' +
                        '[arcsec]\n'+
                        '# [5]: exposure time (s)\n'+
                        '# [6]: photometric catalog\n' + 
                        '# [7]: photometric band\n' +
                        '# [8]: Source Extractor flag\n' + 
                        '# [9]: telescope/instrument\n')
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
    parser.add_argument('-fixed_coo', help='target RA and DEC (degrees)', 
                        nargs=2, default=[0,0])
    parser.add_argument('-positions', help='positions file', default=None)
    parser.add_argument('images', help='images to process', nargs='+')
    args = parser.parse_args()
    man_targetname = args.target
    man_offset = [float(coo) for coo in args.offset]
    fixed_coo = [float(coo) for coo in args.fixed_coo]
    posfile = args.positions
    filenames = args.images

    # check if input filenames is actually a list
    if len(filenames) == 1:
        if filenames[0].find('.lst') > -1 or filenames[0].find('.list') > -1:
            filenames = [filename[:-1] for filename in open(filenames[0], 'r').\
                         readlines()]

    distillate = distill(filenames, man_targetname, man_offset, fixed_coo,
                         posfile, display=True, diagnostics=True)


