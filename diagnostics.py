""" DIAGNOSTICS - diagnostic routines for photometry pipeline
    v1.0: 2016-02-25, mommermiscience@gmail.com
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

import os
import sys
import numpy as np
import logging
import subprocess

from astropy.io import fits
from astropy import wcs
from astropy.visualization import (ZScaleInterval, ImageNormalize,
                                   LogStretch, LinearStretch)
from astropy.time import Time

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pylab as plt
    matplotlib.rcdefaults()  # restore default parameters
except ImportError:
    print('Module matplotlib not found. Please install with: pip install '
          'matplotlib')
    sys.exit()

try:
    from skimage.transform import resize
except ImportError:
    print('Module skimage not found. Please install with: pip install skimage')
    sys.exit()

# pipeline-specific modules
import _pp_conf
import toolbox
from catalog import *

# setup logging
logging.basicConfig(filename=_pp_conf.log_filename,
                    level=_pp_conf.log_level,
                    format=_pp_conf.log_formatline,
                    datefmt=_pp_conf.log_datefmt)


# diagnostics guidelines
#
# - diagnostics.html goes into the data directory
# - all supplementary data and website go into diagroot (see _pp_init.py)
# - if there are sub-directories with data, create diagnostics.html in
#   each directory with data; summary.html links to other directories
###

class Diagnostics_Html():

    from pp_setup import confdiagnostics as conf

    def create_website(self, filename, content='',
                       content_dir='.diagnostics'):
        """
        create empty website for diagnostics output
        """
        html = ("<!DOCTYPE html PUBLIC '-//W3C//DTD HTML 4.01//EN'>\n"
                "<HTML>\n"
                "<HEAD>\n"
                "  <TITLE>Photometry Pipeline - Diagnostics</TITLE>\n"
                "  <LINK rel=\"stylesheet\" href=\"{:s}"
                "diagnostics_stylesheet.css\">\n"
                "</HEAD>\n"
                "<BODY>\n"
                "{:s}\n"
                "</BODY>\n"
                "</HTML>\n").format(os.getenv('PHOTPIPEDIR'), content)

        outf = open(filename, 'w')
        outf.writelines(html)
        outf.close()

    def append_website(self, filename, content, insert_at='</BODY>',
                       replace_from='X?!do not replace anything!?X'):
        """
        append content to an existing website:
        insert content before line that contains `insert_at`
        replace lines between `replace_from` and `insert_at` (by
        default, nothing is replaced)
        """
        # read existing code
        existing_html = open(filename, 'r').readlines()

        # insert content into existing html
        outf = open(filename, 'w')
        delete = False
        for line in existing_html:
            if replace_from in line:
                delete = True
                continue
            if insert_at in line:
                outf.writelines(content)
                delete = False
            if delete:
                continue
            outf.writelines(line)
        outf.close()


class Prepare_Diagnostics(Diagnostics_Html):
    """diagnostics run as part of pp_prepare"""

    function_tag = "<!-- pp_prepare -->"

    def frame_table(self, filenames, obsparam):

        logging.info('create data summary table')

        if self.conf.individual_frame_pages:
            self.frame_pages(filenames, obsparam)

        # create frame information table
        html = "<P><TABLE CLASS=\"gridtable\">\n"
        html += ("<TR><TH>Idx</TH>"
                 "<TH>Filename</TH>"
                 "<TH>Observation Midtime (UT)</TH>"
                 "<TH>Object Name</TH>"
                 "<TH>Airmass</TH>"
                 "<TH>Exptime (s)</TH>"
                 "<TH>Pixel Size (\")"
                 "<TH>Binning</TH>"
                 "<TH>FoV (')</TH></TR>\n")

        for idx, filename in enumerate(filenames):
            hdulist = fits.open(filename, ignore_missing_end=True)
            header = hdulist[0].header
            binning = toolbox.get_binning(header, obsparam)
            try:
                objectname = header[obsparam['object']]
            except KeyError:
                objectname = 'Unknown Target'

            if self.conf.individual_frame_pages:
                framename = "<A HREF=\"{:s}\">{:s}</A>".format(
                    '.diagnostics/'+filename+'.html', filename)
                self.frame_preview(filename)

                # update frame page
                framehtml = ("<!-- Quickview -->\n"
                             "<A HREF=\"#quickview\" "
                             "ONCLICK=\"toggledisplay"
                             "('quickview');\"><H2>Quickview Image</H2>"
                             "</A>\n"
                             "<IMG ID=\"quickview\" SRC=\"{:s}\" "
                             "STYLE=\"display: none\"\>\n\n").format(
                    filename+'.png')
                self.append_website(
                    '.diagnostics/{:s}.html'.format(filename),
                    framehtml, replace_from='<!-- Quickview -->')
            else:
                framename = filename

            html += ("<TR><TD>{:d}</TD>"
                     "<TD>{:s}</TD>"
                     "<TD>{:s}</TD>"
                     "<TD>{:s}</TD>"
                     "<TD>{:4.2f}</TD>"
                     "<TD>{:.1f}</TD>"
                     "<TD>{:.2f} x {:.2f}</TD>"
                     "<TD>{:d} x {:d}</TD>"
                     "<TD>{:.1f} x {:.1f}</TD>\n"
                     "</TR>\n").format(
                         idx+1, framename,
                         Time(header["MIDTIMJD"], format='jd').iso,
                         str(objectname),
                         float(header[obsparam['airmass']]),
                         float(header[obsparam['exptime']]),
                         obsparam['secpix'][0],
                         obsparam['secpix'][1],
                         binning[0], binning[1],
                         float(header[obsparam['extent'][0]]) *
                         obsparam['secpix'][0]*binning[0]/60.,
                         float(header[obsparam['extent'][1]]) *
                         obsparam['secpix'][1]*binning[1]/60.)

        html += '</TABLE>\n'

        return html

    def frame_preview(self, filename):
        """create preview image for one frame"""

        logging.info('create image preview for file {:s}'.format(
            filename))

        hdulist = fits.open(filename, ignore_missing_end=True)

        # create frame image
        imgdat = hdulist[0].data
        # clip extreme values
        # imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
        #                  np.percentile(imgdat, 99))

        # normalize imgdat to pixel values 0 < px < 1
        imgdat = (imgdat - np.min(imgdat)) / np.max(imgdat)
        # resize image larger than 1000px on one side
        imgdat = resize(imgdat,
                        (min(imgdat.shape[0], 1000),
                         min(imgdat.shape[1], 1000)))

        norm = ImageNormalize(
            imgdat, interval=ZScaleInterval(),
            stretch={'linear': LinearStretch(),
                     'log': LogStretch()}[self.conf.image_stretch])

        plt.figure(figsize=(5, 5))

        img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
        # remove axes
        plt.axis('off')
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)

        framefilename = '.diagnostics/' + filename + '.png'
        plt.savefig(framefilename, format='png', bbox_inches='tight',
                    pad_inches=0, dpi=200)
        logging.info('image preview for file {:s} written to {:s}'.format(
            filename, os.path.abspath('.diagnostics/' + filename + '.png')))

        plt.close()
        hdulist.close()

    def frame_pages(self, filenames, obsparam):

        logging.info('setting up individual frame diagnostics report pages')

        for filename in filenames:
            header = fits.open(filename)[0].header
            html = ("<script>\n"
                    "  function toggledisplay(elementID)\n"
                    "  {\n"
                    "  (function(style) {\n"
                    "  style.display = style.display === 'none' ? '' :"
                    "'none';\n"
                    "  })(document.getElementById(elementID).style);\n"
                    "  }\n"
                    "</script>\n\n")
            html += ("<H1>{:s} Diagnostics</H1>"
                     "<P><TABLE CLASS=\"gridtable\">\n"
                     "<TR><TH>Telescope/Instrument</TH><TD>{:s} ({:s})</TD>"
                     "</TR>\n"
                     "<TR><TH>Target/Field Identifier</TH><TD>{:s}</TD>"
                     "</TR>\n"
                     "<TR><TH>RA</TH><TD>{:s}</TD></TR>\n"
                     "<TR><TH>Dec</TH><TD>{:s}</TD></TR>\n"
                     "<TR><TH>Exposure Time (s)</TH><TD>{:s}</TD></TR>\n"
                     "<TR><TH>Observation Midtime</TH><TD>{:s}</TD></TR>\n"
                     "</TABLE>\n\n").format(
                         filename,
                         obsparam['telescope_instrument'],
                         obsparam['telescope_keyword'],
                         header[obsparam['object']],
                         str(header[obsparam['ra']]),
                         str(header[obsparam['dec']]),
                         str(header[obsparam['exptime']]),
                         str(Time(header['MIDTIMJD'], format='jd').iso),
            )

            self.create_website('.diagnostics/{:s}.html'.format(filename),
                                html)
            logging.info(('diagnostics report page for file {:s} '
                          'written to {:s}').format(
                              filename,
                              '.diagnostics/{:s}.html'.format(filename)))

    def add_index(self, filenames, directory, obsparam):
        """
        create index.html
        diagnostic root website
        """

        logging.info('create frame table')

        # create header information
        refheader = fits.open(filenames[0],
                              ignore_missing_end=True)[0].header
        raw_filtername = refheader[obsparam['filter']]
        translated_filtername = obsparam['filter_translations'][
            refheader[obsparam['filter']]]

        html = ("{:s}\n<H1>Photometry Pipeline Diagnostic Output</H1>\n"
                "<TABLE CLASS=\"gridtable\">\n"
                "  <TR><TH>Data Directory</TH><TD>{:s}</TD></TR>\n"
                "  <TR><TH>Telescope/Instrument</TH><TD>{:s}</TD></TR>\n"
                "  <TR><TH>Number of Frames</TH><TD>{:d}</TD></TR>\n"
                "  <TR><TH>Raw Filter Identifier</TH><TD>{:s}</TD></TR>\n"
                "  <TR><TH>Translated Filter Identifier</TH>"
                "<TD>{:s}</TD></TR>\n"
                "  <TR><TH>Log File</TH>"
                "      <TD><A HREF=\"{:s}\">available here</A></TD></TR>"
                "</TABLE>\n").format(
                    self.function_tag,
                    directory,
                    obsparam['telescope_instrument'],
                    len(filenames),
                    raw_filtername,
                    translated_filtername,
                    ('.diagnostics/' +
                     _pp_conf.log_filename.split('.diagnostics/')[1]))

        html += "<H3>Data Summary</H3>\n"
        html += self.frame_table(filenames, obsparam)

        self.create_website(_pp_conf.index_filename, html)


# registration results website

class Registration_Diagnostics(Diagnostics_Html):

    function_tag = "<!-- pp_register -->"

    def registration_table(self, data, extraction_data, obsparam):

        logging.info('creating image registration overview table')

        html = ("<TABLE CLASS=\"gridtable\">\n<TR>\n"
                "<TH>Filename</TH><TH>C<SUB>AS</SUB></TH>"
                "<TH>C<SUB>XY</SUB></TH>"
                "<TH>&sigma;<SUB>RA</SUB> (arcsec)</TH>"
                "<TH>&sigma;<SUB>DEC</SUB> (arcsec)</TH>"
                "<TH>&chi;<SUP>2</SUP><SUB>Reference</SUB></TH>"
                "<TH>&chi;<SUP>2</SUP><SUB>Internal</SUB></TH>\n</TR>\n")

        for dat in data['fitresults']:
            framefilename = '.diagnostics/{:s}.html'.format(dat[0])
            filename = '<A HREF=\"{:s}\">{:s}</A>'.format(
                framefilename, dat[0])

            html += ("<TR><TD>{:s}</TD>"
                     + "<TD>{:4.1f}</TD><TD>{:4.1f}</TD>"
                     + "<TD>{:5.3f}</TD><TD>{:5.3f}</TD>"
                     + "<TD>{:e}</TD><TD>{:e}</TD>\n</TR>\n").format(
                         filename, dat[1], dat[2], dat[3],
                         dat[4], dat[5], dat[6])
        html += "</TABLE>\n"
        html += ("<P CLASS=\"caption\"><STRONG>Legend</STRONG>: "
                 "C<SUB>AS</SUB>: position "
                 "angle/scale contrast (values >{:.1f} are ok); ").format(
                     _pp_conf.scamp_as_contrast_limit)
        html += ("C<SUB>XY</SUB>: xy-shift contrast "
                 "(values >{:.1f} are ok); ").format(
            _pp_conf.scamp_xy_contrast_limit)
        html += ("&sigma;<SUB>RA</SUB> and &sigma;<SUB>DEC</SUB> "
                 "refer to the internal astrometric uncertainties as "
                 "provided by SCAMP; &chi;<SUP>2</SUP><SUB>Reference</SUB> "
                 "and &chi;<SUP>2</SUP><SUB>Internal</SUB> refer to the "
                 "&chi;<SUP>2</SUP> statistics based on the reference "
                 "catalog and the respective frame as provided by SCAMP."
                 "</P>\n")

        return html

    def registration_maps(self, data, extraction_data, obsparam):

        logging.info('create registration maps with reference stars')

        # load reference catalog
        refcat = catalog(data['catalog'])
        for filename in os.listdir('.'):
            if data['catalog'] in filename and '.cat' in filename:
                refcat.read_ldac(filename)
                break

        # create frame images
        for dat in extraction_data:
            framefilename = '.diagnostics/{:s}_astrometry.png'.format(
                dat['fits_filename'])
            imgdat = fits.open(dat['fits_filename'],
                               ignore_missing_end=True)[0].data
            resize_factor = min(1., 1000./np.max(imgdat.shape))

            header = fits.open(dat['fits_filename'],
                               ignore_missing_end=True)[0].header

            # turn relevant header keys into floats
            # astropy.io.fits bug
            for key, val in list(header.items()):
                if 'CD1_' in key or 'CD2_' in key or \
                   'CRVAL' in key or 'CRPIX' in key or \
                   'EQUINOX' in key:
                    header[key] = float(val)

            plt.figure(figsize=(5, 5))
            # create fake image to ensure image dimensions and margins
            img = plt.imshow(np.ones((1000, 1000))*np.nan, origin='lower')

            # remove axes
            plt.axis('off')
            img.axes.get_xaxis().set_visible(False)
            img.axes.get_yaxis().set_visible(False)

            # plot reference sources
            if refcat.shape[0] > 0:
                try:
                    w = wcs.WCS(header)
                    world_coo = np.array(list(zip(refcat['ra_deg'],
                                                  refcat['dec_deg'])))
                    img_coo = w.wcs_world2pix(world_coo, True)
                    img_coo = [c for c in img_coo
                               if (c[0] > 0 and c[1] > 0 and
                                   c[0] < header[obsparam['extent'][0]] and
                                   c[1] < header[obsparam['extent'][1]])]
                    plt.scatter([c[0]*resize_factor for c in img_coo],
                                [c[1]*resize_factor for c in img_coo],
                                s=5, marker='o', edgecolors='red',
                                linewidth=0.3, facecolor='none')
                except astropy.wcs._wcs.InvalidTransformError:
                    logging.error('could not plot reference sources due to '
                                  'astropy.wcs._wcs.InvalidTransformError; '
                                  'most likely unknown distortion '
                                  'parameters.')

            plt.savefig(framefilename, bbox_inches='tight',
                        pad_inches=0, dpi=200, transparent=True)
            logging.info(('registration map image file for image {:s} '
                          'written to {:s}').format(
                              filename, os.path.abspath(framefilename)))

            plt.close()

    def add_registration(self, data, extraction_data):
        """
        add registration results to website
        """
        obsparam = extraction_data[0]['parameters']['obsparam']

        # update index.html
        html = self.function_tag+'\n'
        html += ('<H2>Registration</H2>\n'
                 '<P>Registration based on {:s} catalog: ').format(
            data['catalog'])
        if len(data['badfits']) == 0:
            html += ('<STRONG><FONT COLOR="GREEN">All frames registered '
                     'successfully</FONT></STRONG></P>')
        else:
            html += ('<STRONG><FONT COLOR="RED">{:d} files could not be '
                     'registered</FONT></STRONG></P>').format(
                         len(data['badfits']))

        if self.conf.show_registration_table:
            html += self.registration_table(data, extraction_data, obsparam)

        if (self.conf.individual_frame_pages and
                self.conf.show_registration_star_map):
            self.registration_maps(data, extraction_data, obsparam)

            for framedata in data['fitresults']:
                # update frame page
                filename = framedata[0]
                if filename in data['goodfits']:
                    resultstring = ('<FONT COLOR="GREEN">Registered '
                                    'successfully</FONT>')
                else:
                    resultstring = ('<FONT COLOR="RED">Registration '
                                    'faild</FONT>')

                framehtml = (
                    "<!-- Registration -->\n"
                    "<A HREF=\"#registration\" "
                    "ONCLICK=\"toggledisplay('registration');\">"
                    "<H2>Astrometric Registration</H2></A>\n"
                    "<DIV ID=\"registration\" STYLE=\"display: none\">\n"
                    "<TABLE CLASS=\"gridtable\">\n<TR>\n"
                    "<TH>Filename</TH><TH>C<SUB>AS</SUB></TH>"
                    "<TH>C<SUB>XY</SUB></TH>"
                    "<TH>&sigma;<SUB>RA</SUB> (arcsec)</TH>"
                    "<TH>&sigma;<SUB>DEC</SUB> (arcsec)</TH>"
                    "<TH>&chi;<SUP>2</SUP><SUB>Reference</SUB></TH>"
                    "<TH>&chi;<SUP>2</SUP><SUB>Internal</SUB></TH>\n</TR>\n"
                    "<TR><TD>{:s}</TD>"
                    "<TD>{:4.1f}</TD><TD>{:4.1f}</TD>"
                    "<TD>{:5.3f}</TD><TD>{:5.3f}</TD>"
                    "<TD>{:e}</TD><TD>{:e}</TD>\n</TR>\n"
                    "</TABLE>\n"
                    "<STRONG>{:s}</STRONG>"
                    "<DIV CLASS=\"parent_image\">\n"
                    "  <IMG CLASS=\"back_image\" SRC=\"{:s}\" />\n"
                    "  <IMG CLASS=\"front_image\" SRC=\"{:s}\" />\n"
                    "</DIV>\n</DIV>\n\n").format(
                        filename, framedata[1], framedata[2], framedata[3],
                        framedata[4], framedata[5], framedata[6],
                        resultstring,
                        filename+'.png', filename+"_astrometry.png")
                self.append_website(
                    '.diagnostics/{:s}.html'.format(filename),
                    framehtml, replace_from='<!-- Registration -->')

        self.append_website(_pp_conf.index_filename, html,
                            replace_from=self.function_tag)


class Photometry_Diagnostics(Diagnostics_Html):

    function_tag = "<!-- pp_photometry -->"

    def curve_of_growth_plot(self, data):
        parameters = data['parameters']
        growth_filename = '.diagnostics/curve_of_growth.png'

        f, (ax1, ax2) = plt.subplots(2, sharex=True)

        ax1.set_xlim([min(parameters['aprad']), max(parameters['aprad'])])
        ax1.set_ylabel('Fractional Combined Flux')
        if not parameters['target_only']:
            ax1.plot(parameters['aprad'], data['background_flux'][0],
                     color='black', linewidth=1,
                     label='background sources')
            ax1.fill_between(parameters['aprad'],
                             (data['background_flux'][0] -
                              data['background_flux'][1]),
                             (data['background_flux'][0] +
                              data['background_flux'][1]),
                             color='black', alpha=0.2)
        if not parameters['background_only']:
            ax1.plot(parameters['aprad'], data['target_flux'][0],
                     color='red', linewidth=1,
                     label='target')
            ax1.fill_between(parameters['aprad'],
                             (data['target_flux'][0] -
                              data['target_flux'][1]),
                             (data['target_flux'][0] +
                              data['target_flux'][1]),
                             color='red', alpha=0.2)
        ax1.set_ylim([0, ax1.get_ylim()[1]])
        ax1.plot([data['optimum_aprad'], data['optimum_aprad']],
                 [ax1.get_ylim()[0], ax1.get_ylim()[1]],
                 linewidth=2, color='blue')
        ax1.plot([plt.xlim()[0], plt.xlim()[1]],
                 [data['fluxlimit_aprad'], data['fluxlimit_aprad']],
                 color='black', linestyle='--')
        ax1.grid()
        ax1.legend(loc=4)

        ax2.set_ylim([-0.1, 1.1])
        ax2.set_ylabel('SNR')
        if not parameters['target_only']:
            ax2.errorbar(parameters['aprad'], data['background_snr'],
                         color='black', linewidth=1)
        if not parameters['background_only']:
            ax2.errorbar(parameters['aprad'], data['target_snr'],
                         color='red', linewidth=1)
        ax2.plot([data['optimum_aprad'], data['optimum_aprad']],
                 [plt.ylim()[0], plt.ylim()[1]],
                 linewidth=2, color='blue')
        ax2.grid()
        ax2.set_xlabel('Aperture Radius (px)')
        plt.savefig(growth_filename, format='png')
        plt.close()
        data['growth_filename'] = growth_filename

    def fwhm_vs_time_plot(self, extraction, data):
        fwhm_filename = '.diagnostics/fwhm.png'

        frame_midtimes = np.array([frame['time'] for frame in extraction])
        fwhm = [np.median(frame['catalog_data']['FWHM_IMAGE'])
                for frame in extraction]
        fwhm_sig = [np.std(frame['catalog_data']['FWHM_IMAGE'])
                    for frame in extraction]

        plt.title('Median PSF FWHM per Frame')
        plt.xlabel('Minutes after {:s} UT'.format(
            Time(frame_midtimes.min(), format='jd',
                 out_subfmt='date_hm').iso))
        plt.ylabel('Point Source FWHM (px)')
        plt.scatter((frame_midtimes-frame_midtimes.min())*1440,
                    fwhm, marker='o',
                    color='black')
        xrange = [plt.xlim()[0], plt.xlim()[1]]
        plt.plot(xrange, [data['optimum_aprad']*2, data['optimum_aprad']*2],
                 color='blue')
        plt.xlim(xrange)
        plt.ylim([0, max([data['optimum_aprad']*2+1, max(fwhm)])])

        plt.grid()
        plt.savefig(fwhm_filename, format='png')
        plt.close()
        data['fwhm_filename'] = fwhm_filename

    def add_photometry(self, data, extraction):
        """
        add photometry results to website
        """
        # create curve-of-growth plot
        self.curve_of_growth_plot(data)

        # create fwhm vs time plot
        self.fwhm_vs_time_plot(extraction, data)

        # update index.html
        html = self.function_tag+'\n'
        html += ("<H2>Instrumental Photometry</H2>\n"
                 "<TABLE CLASS=\"gridtable\">\n"
                 "<TR><TH>Photometry Method</TH><TD>{:s}</TD></TR>\n"
                 "<TR><TH>Source Extractor MINAREA (px)</TH>"
                 "<TD>{:.1f}</TD></TR>\n"
                 "<TR><TH>Source Extractor Detection Threshold (&sigma;)"
                 "</TH><TD>{:.1f}</TD></TR>\n").format(
                     {'APER': 'Aperture Photometry'}[_pp_conf.photmode],
                     extraction[0]['parameters']['source_minarea'],
                     extraction[0]['parameters']['sex_snr'])

        if _pp_conf.photmode == 'APER':

            if data['n_target'] > 0 and data['n_bkg'] > 0:
                apsrc = ("{:d} target detections and {:d} "
                         "background detections").format(
                    data['n_target'], data['n_bkg'])
            elif data['n_target'] == 0 and data['n_bkg'] > 0:
                apsrc = "{:d} frames with background detections".format(
                    data['n_bkg'])
            elif data['n_bkg'] == 0 and data['n_target'] > 0:
                apsrc = "{:d} frames with target detections".format(
                    data['n_target'])
            else:
                apsrc = "manually defined"

            html += ("<TR><TH>Aperture Radius (px)</TH>"
                     "<TD>{:.2f}</TD></TR>\n"
                     "<TR><TH>Aperture Radius Basis</TH>"
                     "<TD>{:s}</TD></TR>\n"
                     "<TR><TH>Aperture Radius Strategy</TH>"
                     "<TD>{:s}</TD></TR>\n").format(
                         data['optimum_aprad'],
                         apsrc,
                         data['aprad_strategy']
            )

        html += "</TABLE>\n"

        html += "<P><IMG SRC=\"{:s}\">\n".format(data['growth_filename'])
        html += "<IMG SRC=\"{:s}\">\n".format(data['fwhm_filename'])

        self.append_website(_pp_conf.index_filename, html,
                            replace_from=self.function_tag)


class Calibration_Diagnostics(Diagnostics_Html):

    function_tag = "<!-- pp_calibrate -->"

    def zeropoint_overview_plot(self, data):
        """produce a plot of magnitude zeropoint as a function of time"""

        logging.info('create zeropoint overview plot')

        times = np.array([dat['obstime'][0] for dat in data['zeropoints']])
        zp = [dat['zp'] for dat in data['zeropoints']]
        zperr = [dat['zp_sig'] for dat in data['zeropoints']]

        plt.plot()
        plt.errorbar((times-times.min())*1440, zp, yerr=zperr, linestyle='',
                     color='blue', marker='s', capsize=3)
        plt.xlabel('Minutes after {:s} UT'.format(
            Time(times.min(), format='jd',
                 out_subfmt='date_hm').iso))
        plt.ylabel(
            '{:s}-Band Magnitude Zeropoints (mag)'.format(
                data['filtername']))
        plt.ylim([plt.ylim()[1], plt.ylim()[0]])
        plt.grid()
        plt.savefig('.diagnostics/zeropoints.png', format='png')
        logging.info('zeropoint overview plot written to {:s}'.format(
            os.path.abspath('.diagnostics/zeropoints.png')))
        plt.close()
        data['zpplot'] = 'zeropoints.png'

    def phot_calibration_plot(self, data, idx):
        """produce curve-of-growth plot for each frame"""

        f, (ax1, ax3) = plt.subplots(2)
        plt.subplots_adjust(hspace=0.3)

        ax1.set_title('%s: %s-band from %s' %
                      (data['catalogs'][idx].catalogname,
                       data['filtername'],
                       data['ref_cat'].catalogname))
        ax1.set_xlabel('Number of Reference Stars')
        ax1.set_ylabel('Magnitude Zeropoint', fontdict={'color': 'red'})

        zp_idx = data['zeropoints'][idx]['zp_idx']
        clipping_steps = data['zeropoints'][idx]['clipping_steps']

        x = [len(clipping_steps[i][3]) for i
             in range(len(clipping_steps))]

        ax1.errorbar(x, [clipping_steps[i][0] for i
                         in range(len(clipping_steps))],
                     yerr=[clipping_steps[i][1] for i
                           in range(len(clipping_steps))], color='red')
        ax1.set_ylim(ax1.get_ylim()[::-1])  # reverse y axis
        ax1.plot([len(clipping_steps[zp_idx][3]),
                  len(clipping_steps[zp_idx][3])],
                 ax1.get_ylim(), color='black')

        ax1.grid(linestyle='--')

        ax2 = ax1.twinx()
        ax2.plot(x, [clipping_steps[i][2] for i
                     in range(len(clipping_steps))],
                 color='blue')
        ax2.set_ylabel(r'reduced $\chi^2$', fontdict={'color': 'blue'})
        ax2.set_yscale('log')

        # residual plot
        ax3.set_xlabel('Reference Star Magnitude')
        ax3.set_ylabel('Calibration-Reference (mag)')

        match = data['zeropoints'][idx]['match']
        x = match[0][0][clipping_steps[zp_idx][3]]
        residuals = (match[1][0][clipping_steps[zp_idx][3]]
                     + clipping_steps[zp_idx][0]
                     - match[0][0][clipping_steps[zp_idx][3]])
        residuals_sig = np.sqrt(match[1][1][clipping_steps[zp_idx][3]]**2
                                + clipping_steps[zp_idx][1]**2)

        ax3.errorbar(x, residuals, yerr=residuals_sig, color='black',
                     marker='o', linestyle='')
        x_range = ax3.get_xlim()
        ax3.plot(x_range, [0, 0], color='black', linestyle='--')
        ax3.set_xlim(x_range)
        ax3.set_ylim(ax3.get_ylim()[::-1])  # reverse y axis

        ax3.grid(linestyle='--')

        plotfilename = '.diagnostics/{:s}_photcal.png'.format(
            data['catalogs'][idx].catalogname)
        plt.savefig(plotfilename, format='png')
        data['zeropoints'][idx]['plotfilename'] = plotfilename

    def calibration_raw_data_tables(self, dat):

        html = "<TD><TABLE CLASS=\"gridtable\">\n<TR>\n"
        html += ("<TH>Idx</TH><TH>Source Name</TH><TH>RA</TH><TH>Dec</TH>"
                 "<TH>Catalog (mag)</TH>"
                 "<TH>Instrumental (mag)</TH><TH>Calibrated (mag)</TH>"
                 "<TH>Residual (mag</TH>\n</TR>\n")
        for i, idx in enumerate(dat['zp_usedstars']):
            name = str(dat['match'][0][2][idx])
            if isinstance(name, bytes):
                name = name.decode('utf8')
            html += ("<TR><TD>{:d}</TD><TD>{:s}</TD><TD>{:12.8f}</TD>"
                     "<TD>{:12.8f}</TD><TD>{:.3f}+-{:.3f}</TD>"
                     "<TD>{:.3f}+-{:.3f}</TD>"
                     "<TD>{:.3f}+-{:.3f}</TD><TD>{:.3f}</TD>"
                     "</TR>").format(
                         i+1, name,
                         dat['match'][0][3][idx],
                         dat['match'][0][4][idx],
                         dat['match'][0][0][idx],
                         dat['match'][0][1][idx],
                         dat['match'][1][0][idx],
                         dat['match'][1][1][idx],
                         dat['zp']+dat['match'][1][0][idx],
                         np.sqrt(dat['zp_sig']**2 +
                                 dat['match'][1][1][idx]**2),
                         (dat['zp']+dat['match'][1][0][idx]) -
                         dat['match'][0][0][idx])
        html += "</TABLE>\n"

        return html

    def calibration_star_maps(self, dat):
        """create thumbnail images with calibration stars marked"""

        fits_filename = (dat['filename'][:dat['filename'].find('.ldac')]
                         + '.fits')
        imgdat = fits.open(fits_filename,
                           ignore_missing_end=True)[0].data
        resize_factor = min(1., 1000./np.max(imgdat.shape))

        # normalize imgdat to pixel values 0 < px < 1
        if np.min(imgdat) < 0:
            imgdat = imgdat + np.min(imgdat)
        if np.max(imgdat) > 1:
            imgdat = imgdat / np.max(imgdat)

        imgdat = resize(imgdat,
                        (min(imgdat.shape[0], 1000),
                         min(imgdat.shape[1], 1000)))
        header = fits.open(fits_filename,
                           ignore_missing_end=True)[0].header

        norm = ImageNormalize(
            imgdat, interval=ZScaleInterval(),
            stretch={'linear': LinearStretch(),
                     'log': LogStretch()}[self.conf.image_stretch])

        # turn relevant header keys into floats
        for key, val in list(header.items()):
            if 'CD1_' in key or 'CD2_' in key or \
               'CRVAL' in key or 'CRPIX' in key or \
               'EQUINOX' in key:
                header[key] = float(val)

        plt.figure(figsize=(4, 4))
        img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')

        # remove axes
        plt.axis('off')
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)

        # plot reference sources
        if len(dat['match'][0][3]) > 0 and len(dat['match'][0][4]) > 0:
            try:
                w = wcs.WCS(header)
                world_coo = [[dat['match'][0][3][idx],
                              dat['match'][0][4][idx]]
                             for idx in dat['zp_usedstars']]
                img_coo = w.wcs_world2pix(world_coo, True)
                plt.scatter([c[0]*resize_factor for c in img_coo],
                            [c[1]*resize_factor for c in img_coo],
                            s=10, marker='o', edgecolors='red',
                            linewidth=0.3,
                            facecolor='none')
                for i in range(len(dat['zp_usedstars'])):
                    plt.annotate(str(i+1),
                                 xy=((img_coo[i][0]*resize_factor)+15,
                                     img_coo[i][1]*resize_factor),
                                 color='red',
                                 horizontalalignment='left',
                                 verticalalignment='center')
            except astropy.wcs._wcs.InvalidTransformError:
                logging.error('could not plot reference sources due to '
                              'astropy.wcs._wcs.InvalidTransformError; '
                              'most likely unknown distortion '
                              'parameters.')

        catframe = ('.diagnostics/{:s}.'
                    'fits_reference_stars.png').format(
            dat['filename'][:dat['filename'].find('.ldac')])
        plt.savefig(catframe, format='png', bbox_inches='tight',
                    pad_inches=0, dpi=200)
        plt.close()

    def add_calibration(self, data, instrumental=False):
        """
        wrapper to add calibration results to diagnostics website
        """

        html = self.function_tag+'\n'
        html += "<H2>Photometric Calibration</H2>\n"

        if not instrumental:
            # create zeropoint overview plot
            self.zeropoint_overview_plot(data)

            # main diagnostics website content
            html += ("<TABLE CLASS=\"gridtable\">\n"
                     "<TR><TH>Reference Catalog</TH><TD>{:s}</TD></TR>\n"
                     "<TR><TH>Reference Catalog History</TH>"
                     "<TD>{:s}</TD></TR>\n"
                     "<TR><TH>Target Filter</TH><TD>{:s}</TD></TR>\n"
                     "</TABLE>\n").format(
                         data['ref_cat'].catalogname,
                         data['ref_cat'].history,
                         data['filtername'])

            # build overview table
            html += ("<P><TABLE CLASS=\"gridtable\">\n<TR>\n"
                     "<TH>Filename</TH><TH>Zeropoint (mag)</TH>"
                     "<TH>&sigma; (mag)</TH>"
                     "<TH>N<SUP>*</SUP><SUB>used</SUB></TH>"
                     "<TH>N<SUP>*</SUP><SUB>matched</SUB></TH>\n</TR>\n")
            for idx, dat in enumerate(data['zeropoints']):

                # update frame pages
                if self.conf.individual_frame_pages:
                    framename = "<A HREF=\"{:s}\">{:s}</A>".format(
                        '.diagnostics/'+dat['filename'][:-4]+'fits'+'.html',
                        dat['filename'][:-4]+'fits')

                    framehtml = ("<!-- Calibration -->\n"
                                 "<A HREF=\"#calibration_overview\" "
                                 "ONCLICK=\"toggledisplay"
                                 "('calibration_overview');\">"
                                 "<H2>Photometric Calibration</H2></A>\n"
                                 "<DIV ID=\"calibration_overview\" "
                                 "STYLE=\"display: none\"\>\n")

                    framehtml += ("<P><TABLE CLASS=\"gridtable\">\n"
                                  "<TR><TH>Reference Catalog</TH>"
                                  "<TD>{:s}</TD></TR>\n"
                                  "<TR><TH>Reference Catalog History</TH>"
                                  "<TD>{:s}</TD></TR>\n"
                                  "<TR><TH>Target Filter</TH>"
                                  "<TD>{:s}</TD></TR>\n"
                                  "<TR><TH>Zeropoint (mag)</TH>"
                                  "<TD>{:7.4f}+-{:.4f}</TD></TR>\n"
                                  "<TR><TH>N<SUP>*</SUP><SUB>used</SUB>"
                                  "</TH>"
                                  "<TD>{:d}</TD></TR>\n"
                                  "<TH>N<SUP>*</SUP><SUB>matched</SUB></TH>"
                                  "<TD>{:d}</TD></TR>\n").format(
                                      data['ref_cat'].catalogname,
                                      data['ref_cat'].history,
                                      data['filtername'],
                                      dat['zp'], dat['zp_sig'],
                                      dat['zp_nstars'],
                                      len(dat['match'][0][0]))
                    framehtml += "</TABLE></P>\n"

                    # frame calibration data
                    catframe = ('.diagnostics/{:s}.'
                                'fits_reference_stars.png').format(
                                    dat['filename'][:dat['filename'].find(
                                        '.ldac')])

                    # build individual curve of growth plots
                    if self.conf.show_phot_calibration_plots:
                        self.phot_calibration_plot(data, idx)
                        framehtml += (
                            "<A HREF=\"#calibration_plot\" "
                            "ONCLICK=\"toggledisplay"
                            "('calibration_plot');\">"
                            "<H3>Calibration Analysis</H3></A>\n"
                            "<DIV ID=\"calibration_plot\">\n"
                            "<P><IMG SRC={:s} \></DIV>\n").format(
                            dat['plotfilename'].split('.diagnostics/')[1])

                    # build individual catalog maps
                    if self.conf.show_calibration_star_map:
                        self.calibration_star_maps(dat)
                        framehtml += (
                            "<A HREF=\"#calibration_starmap\" "
                            "ONCLICK=\"toggledisplay"
                            "('calibration_starmap');\">"
                            "<H3>Calibration Map</H3></A>\n"
                            "<DIV ID=\"calibration_starmap\" "
                            "STYLE=\"display: none\"\>\n"
                            "<P><IMG SRC={:s} \></DIV>\n").format(
                            catframe.split('.diagnostics/')[1])

                    # build individual catalog data table websites
                    if self.conf.show_calibration_star_table:
                        framehtml += (
                            "<A HREF=\"#calibration_table\" "
                            "ONCLICK=\"toggledisplay"
                            "('calibration_table');\">"
                            "<H3>Calibration Data Table</H3></A>\n"
                            "<DIV ID=\"calibration_table\" "
                            "STYLE=\"display: none\"\>\n"
                            "<P>{:s}</DIV>\n").format(
                            self.calibration_raw_data_tables(dat))

                    framehtml += "</DIV>\n\n"

                    self.append_website(
                        '.diagnostics/{:s}.html'.format(
                            dat['filename'][:-4]+'fits'),
                        framehtml, replace_from='<!-- Calibration -->')

                else:
                    framename = dat['filename'][:-4]+'fits'

                html += ("<TR><TD>{:s}</TD>"
                         "<TD>{:7.4f}</TD><TD>{:7.4f}</TD><TD>{:d}</TD>"
                         + "<TD>{:d}</TD>\n</TR>").format(
                             framename, dat['zp'],
                             dat['zp_sig'], dat['zp_nstars'],
                             len(dat['match'][0][0]))
            html += "</TABLE></P>\n"

            html += "<P><IMG SRC=\"{:s}\" ALT=\"Zeropoints\">\n".format(
                '.diagnostics/'+data['zpplot'])
        else:
            html += ("Instrumental magnitudes are reported "
                     "(filter used: {:s})\n").format(
                         str(data['filtername']))

        self.append_website(_pp_conf.index_filename, html,
                            replace_from=self.function_tag)


class Distill_Diagnostics(Diagnostics_Html):

    function_tag = "<!-- pp_distill -->"

    def lightcurve_plots(self, data):
        data['lightcurveplots'] = {}
        for target in data['targetnames']:

            logging.info('create lightcurve plot for {:s}'.format(target))

            midtimes = np.array([dat[9][0] for dat in data[target]])

            plt.plot()
            plt.title(target.replace('_', ' '))
            plt.xlabel('Minutes after {:s} UT'.format(
                Time(midtimes.min(), format='jd',
                     out_subfmt='date_hm').iso))
            plt.ylabel('Magnitude')
            plt.errorbar((midtimes-midtimes.min())*1440,
                         [dat[7] for dat in data[target]],
                         yerr=[dat[8] for dat in data[target]],
                         linestyle='', color='red',
                         marker='o', capsize=3)
            plt.ylim([plt.ylim()[1], plt.ylim()[0]])
            plt.xticklabels = [Time(t, format='jd').iso
                               for t in plt.xticks()[0]]
            plt.grid()
            plt.savefig('.diagnostics/{:s}.png'.format(
                        target.translate(_pp_conf.target2filename)),
                        format='png')
            logging.info('lightcurve plot for {:s} written to {:s}'.format(
                target, os.path.abspath('.diagnostics/{:s}.png'.format(
                    target.translate(_pp_conf.target2filename)))))
            plt.close()
            data['lightcurveplots'][target] = ('.diagnostics/'
                                               '{:s}.png').format(
                target.translate(_pp_conf.target2filename))

    def thumbnail_images(self, data):
        data['thumbnailplots'] = {}

        boxsize = 300  # thumbnail boxsize
        for target in data['targetnames']:

            if sys.version_info < (3, 0):
                target = str(target)

            data['thumbnailplots'][target] = []
            for dat in data[target]:
                for fitsfilename in ['.fits', '.fit']:
                    fitsfilename = (dat[10][:dat[10].find('.ldac')] +
                                    fitsfilename)
                    if os.path.isfile(fitsfilename):
                        break
                hdulist = fits.open(fitsfilename, ignore_missing_end=True)

                logging.info('create thumbnail image for {:s}/{:s}'.format(
                    target, fitsfilename))

                # turn relevant header keywords into floats
                # should be fixed in astropy.wcs
                for key, val in list(hdulist[0].header.items()):
                    if 'CD1' in key or 'CD2' in key or \
                       'CRVAL' in key or 'CRPIX' in key or \
                       'EQUINOX' in key:
                        hdulist[0].header[key] = float(val)

                w = wcs.WCS(hdulist[0].header)
                obj_x, obj_y = dat[11], dat[12]
                image_coords = w.wcs_world2pix(np.array([[dat[1], dat[2]]]),
                                               True)
                exp_x, exp_y = image_coords[0][0], image_coords[0][1]

                # create margin around image allowing for any cropping
                composite = np.zeros((hdulist[0].data.shape[0]+2*boxsize,
                                      hdulist[0].data.shape[1]+2*boxsize))

                composite[boxsize:boxsize +
                          hdulist[0].data.shape[0],
                          boxsize:boxsize +
                          hdulist[0].data.shape[1]] = hdulist[0].data

                # extract thumbnail data accordingly
                thumbdata = composite[int(boxsize+obj_y-boxsize/2):
                                      int(boxsize+obj_y+boxsize/2),
                                      int(boxsize+obj_x-boxsize/2):
                                      int(boxsize+obj_x+boxsize/2)]

                # run statistics over center of the frame around the target
                if thumbdata.shape[0] > 0 and thumbdata.shape[1] > 0:
                    norm = ImageNormalize(
                        thumbdata, interval=ZScaleInterval(),
                        stretch={'linear': LinearStretch(),
                                 'log': LogStretch()}[
                                     self.conf.image_stretch])

                    # extract aperture radius
                    if _pp_conf.photmode == 'APER':
                        aprad = float(hdulist[0].header['APRAD'])

                    # create plot
                    fig = plt.figure()
                    img = plt.imshow(thumbdata, cmap='gray',
                                     origin='lower', norm=norm)
                    # remove axes
                    plt.axis('off')
                    img.axes.get_xaxis().set_visible(False)
                    img.axes.get_yaxis().set_visible(False)

                    plt.annotate('{:s}\n{:5.3f}+-{:5.3f} mag'.format(
                        fitsfilename, dat[7], dat[8]), (3, 10),
                        color='white')

                    # place aperture
                    if _pp_conf.photmode == 'APER':
                        targetpos = plt.Circle((boxsize/2, boxsize/2),
                                               aprad, ec='red', fc='none',
                                               linewidth=1)
                    else:
                        targetpos = plt.Rectangle(
                            (boxsize/2-7, boxsize/2-7),
                            15, 15, ec='red', fc='none',
                            linewidth=1)
                    plt.gca().add_patch(targetpos)

                    # place expected position (if within thumbnail)
                    if (abs(exp_x-obj_x) <= boxsize/2 and
                            abs(exp_y-obj_y) <= boxsize/2):
                        plt.scatter(exp_x-obj_x+boxsize/2,
                                    exp_y-obj_y+boxsize/2,
                                    marker='+', s=100, color='green')

                    thumbfilename = ('.diagnostics/' +
                                     target.translate(
                                         _pp_conf.target2filename) + '_' +
                                     fitsfilename[:fitsfilename.
                                                  find('.fit')] +
                                     '_thumb.png')
                    plt.savefig(thumbfilename, format='png',
                                bbox_inches='tight',
                                pad_inches=0)
                    plt.close()
                    hdulist.close()
                    data['thumbnailplots'][target].append((fitsfilename,
                                                           thumbfilename))
                else:
                    logging.warning('cannot produce thumbnail image ' +
                                    'for {:s} in frame {:s}').format(
                                        target, dat[10])
                    continue

    def target_animations(self, data):
        data['gifs'] = {}

        for target in data['targetnames']:
            gif_filename = '{:s}.gif'.format(
                target.translate(_pp_conf.target2filename))
            logging.info('converting images to gif: {:s}'.format(
                gif_filename))
            root = os.getcwd()
            os.chdir(_pp_conf.diagroot)
            try:
                convert = subprocess.Popen(
                    ['convert', '-delay', '50',
                     ('{:s}*thumb.png'.format(target.translate(
                         _pp_conf.target2filename))),
                     '-loop', '0',
                     ('%s' % gif_filename)])

                convert.wait()
            except:
                logging.warning('could not produce gif animation for '
                                + 'target {:s}'.format(target))
            data['gifs'][target] = '.diagnostics/' + gif_filename
            os.chdir(root)

    def add_frame_report(self, data):
        data['resultswebsites'] = {}
        for target in data['targetnames']:

            if sys.version_info < (3, 0):
                target = str(target)

            html = "<H2>{:s} - Photometric Results</H2>\n".format(
                target)
            html += "<P><IMG SRC=\"{:s}\">\n".format(
                data['lightcurveplots'][target].split(
                    '.diagnostics/')[1])
            html += "<IMG SRC=\"{:s}\">\n".format(
                data['gifs'][target].split('.diagnostics/')[1])

            # create summary table
            html += "<TABLE CLASS=\"gridtable\">\n<TR>\n"
            html += ("<TH>Filename</TH><TH>Julian Date</TH>"
                     "<TH>Target (mag)</TH>"
                     "<TH>sigma (mag)</TH><TH>Target RA (deg)</TH>"
                     "<TH>Target Dec (deg)</TH><TH>RA Offset (\")</TH>"
                     "<TH>Dec Offset (\")</TH>\n</TR>\n")
            for dat in data[target]:
                html += ("<TR><TD><A HREF=\"#{:s}\">{:s}</A></TD>"
                         "<TD>{:15.7f}</TD><TD>{:7.4f}</TD>"
                         "<TD>{:6.4f}</TD><TD>{:13.8f}</TD>"
                         "<TD>{:+13.8f}</TD><TD>{:5.2f}</TD>"
                         "<TD>{:5.2f}</TD>\n</TR>\n").format(
                             dat[10], dat[10], dat[9][0],
                             dat[7], dat[8], dat[3], dat[4],
                             ((dat[1]-dat[3])*3600.),
                             ((dat[2]-dat[4])*3600.))
            html += "</TABLE>\n"

            # plot individual thumbnails
            html += "<H3>Thumbnails</H3>\n"
            for idx, plts in enumerate(data['thumbnailplots'][target]):
                html += ("<P>{:s}<IMG ID=\"{:s}\" "
                         "SRC=\"{:s}\">\n").format(
                             plts[0],
                             data[target][idx][10],
                             plts[1].split('.diagnostics/')[1])
            filename = ('.diagnostics/' +
                        target.translate(
                            _pp_conf.target2filename) +
                        '_' + 'results.html')
            self.create_website(filename, html)
            data['resultswebsites'][target] = filename

    def add_results(self, data, imagestretch='linear'):
        """
        add results to website
        """

        self.lightcurve_plots(data)

        self.thumbnail_images(data)

        self.target_animations(data)

        self.add_frame_report(data)

        html = self.function_tag+'\n'
        html += "<H2>Photometry Results</H2>\n"
        html += ("<P>photometric data obtained for {:d} "
                 "object(s): \n").format(
                     len(data['targetnames']))
        for target in data['targetnames']:
            print(target)
            html += "<BR><A HREF=\"{:s}\">{:s}</A>\n".format(
                data['resultswebsites'][target], target)
        for target in data['targetnames']:
            html += "<P><IMG SRC=\"{:s}\">\n".format(
                data['lightcurveplots'][target])
            html += "<IMG SRC=\"{:s}\">\n".format(data['gifs'][target])

        self.append_website(_pp_conf.index_filename, html,
                            replace_from=self.function_tag)


def abort(where):
    """
    use this function to add information to index.html that the
    pipeline crashed and where
    """
    html = ("<P><FONT COLOR=\"RED\">Pipeline crashed "
            + "unexpectedly in module %s; refere to <A HREF=\"%s\">log</A> "
            + "for additional information</FONT>\n") % (
                _pp_conf.log_filename, where)

    append_website(_pp_conf.index_filename, html)

    return None


preparation = Prepare_Diagnostics()
registration = Registration_Diagnostics()
photometry = Photometry_Diagnostics()
calibration = Calibration_Diagnostics()
distill = Distill_Diagnostics()
