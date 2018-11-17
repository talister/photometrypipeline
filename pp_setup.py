""" PP_SETUP - photometry pipeline setup file
    v1.0: 2018-11-15, mommermiscience@gmail.com
"""
# Photometry Pipeline
# Copyright (C) 2016-2018  Michael Mommert, mommermiscience@gmail.com

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


"""This configuration file is meant to customize your pipeline
according to your needs. This is not an executable file. Instead, each
pp function will import the classes and use the definitions made in
here - parameters provided to the function through the command line
will override the definitions in this class. Definitions below are
group into different class based on their pp function association.
"""


class Conf():
    """General pipeline configurations"""

    diagnostics = True  # produce diagnostic files and website?


class ConfPrepare(Conf):
    """configuration setup for pp_prepare"""
    pass


class ConfRegister(Conf):
    """configuration setup for pp_register"""
    pass


class ConfPhotometry(Conf):
    """configuration setup for pp_photometry"""
    pass


class ConfCalibrate(Conf):
    """Configuration setup for pp_calibrate"""

    # write photometric calibration raw data into file
    save_caldata = True  # save data on calibration process
    save_caldata_format = 'ascii.basic'  # Table.write formats (ascii, csv...)
    save_caldata_suffix = '_cal.dat'  # file suffix ('.fits' will be clipped)
    save_caldata_usedonly = False  # only output stars used in calibration?

    # add photometric calibration raw data to frame database
    caldata_in_db = True  # add calibration data to database file?


class ConfDistill(Conf):
    """configuration setup for pp_distill"""
    pass


class ConfDiagnostics(Conf):
    """configuration setup for diagnostics"""
    pass


class ConfCombine(Conf):
    """configuration setup for pp_combine"""
    pass


class ConfMPCReport(Conf):
    """configuration setup for pptool_mpcreport"""
    pass


confprepare = ConfPrepare()
confcalibrate = ConfCalibrate()
