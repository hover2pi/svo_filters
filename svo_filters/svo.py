#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A Python wrapper for the SVO Filter Profile Service
"""
from glob import glob
import inspect
import os
from pkg_resources import resource_filename
import warnings
import itertools

import astropy.table as at
import astropy.io.votable as vo
import astropy.units as q
from astropy.utils.exceptions import AstropyWarning
from astroquery.svo_fps import SvoFps
from bokeh.plotting import figure, show
import bokeh.palettes as bpal
import numpy as np


warnings.simplefilter('ignore', category=AstropyWarning)
EXTINCTION = {'PS1.g': 3.384, 'PS1.r': 2.483, 'PS1.i': 1.838, 'PS1.z': 1.414,
              'PS1.y': 1.126, 'SDSS.u': 4.0, 'SDSS.g': 3.384, 'SDSS.r': 2.483,
              'SDSS.i': 1.838, 'SDSS.z': 1.414, '2MASS.J': 0.650,
              '2MASS.H': 0.327, '2MASS.Ks': 0.161}

SYSTEMATICS = {}


class Filter:
    """
    Creates a Filter object to store a photometric filter profile
    and metadata

    Attributes
    ----------
    path: str
        The absolute filepath for the bandpass data, an ASCII file with
        a wavelength column in Angstroms and a response column of values
        ranging from 0 to 1
    refs: list, str
        The references for the bandpass data
    rsr: np.ndarray
        The wavelength and relative spectral response (RSR) arrays
    Band: str
        The band name
    CalibrationReference: str
        The paper detailing the calibration
    FWHM: float
        The FWHM for the filter
    Facility: str
        The telescope facility
    FilterProfileService: str
        The SVO source
    MagSys: str
        The magnitude system
    PhotCalID: str
        The calibration standard
    PhotSystem: str
        The photometric system
    ProfileReference: str
        The SVO reference
    WavelengthCen: float
        The center wavelength
    WavelengthEff: float
        The effective wavelength
    WavelengthMax: float
        The maximum wavelength
    WavelengthMean: float
        The mean wavelength
    WavelengthMin: float
        The minimum wavelength
    WavelengthPeak: float
        The peak wavelength
    WavelengthPhot: float
        The photon distribution based effective wavelength
    WavelengthPivot: float
        The wavelength pivot
    WavelengthUCD: str
        The SVO wavelength unit
    WavelengthUnit: str
        The wavelength unit
    WidthEff: float
        The effective width
    ZeroPoint: float
        The value of the zero point flux
    ZeroPointType: str
        The system of the zero point
    ZeroPointUnit: str
        The units of the zero point
    filterID: str
        The SVO filter ID

    """
    def __init__(self, band, filter_directory=None, wave_units=q.um, flux_units=q.erg/q.s/q.cm**2/q.AA, **kwargs):
        """
        Loads the bandpass data into the Filter object

        Parameters
        ----------
        band: str
            The bandpass filename (e.g. 2MASS.J)
        filter_directory: str
            The directory containing the filter files
        wave_units: str, astropy.units.core.PrefixUnit  (optional)
            The wavelength units
        flux_units: str, astropy.units.core.PrefixUnit  (optional)
            The zeropoint flux units
        """
        if filter_directory is None:
            filter_directory = resource_filename('svo_filters', 'data/filters/')

        # Check if TopHat
        if band.lower().replace('-', '').replace(' ', '') == 'tophat':

            # check kwargs for limits
            wave_min = kwargs.get('wave_min')
            wave_max = kwargs.get('wave_max')
            filepath = ''

            if wave_min is None or wave_max is None:
                raise ValueError("Please provide **{'wave_min', 'wave_max'} to create top hat filter.")

            else:
                # Load the filter
                n_pix = kwargs.get('n_pixels', 100)
                self.load_TopHat(wave_min, wave_max, n_pix)

        else:

            # List of all bands on file
            bands = filters(filter_directory)

            # Read file if the band is in the filter directory
            if band in bands:

                # Get the file
                filepath = glob(os.path.join(filter_directory, band + '*'))[0]

                # Get the first line to determine format
                with open(filepath) as f:
                    top = f.readline()

                # Read in XML file
                if top.startswith('<?xml'):

                    self.load_xml(filepath)

                # Read in txt file
                elif filepath.endswith('.txt'):

                    self.load_txt(filepath)

                else:

                    raise TypeError("File must be XML or ascii format.")

            # Otherwise try a Web query or throw an error
            else:
              
                err = """No filters match {}\n\nFILTERS ON FILE: {}\n\nA full list of available filters from the\nSVO Filter Profile Service can be found at\nhttp: //svo2.cab.inta-csic.es/theory/fps3/\n\nTry again with the desired filter as '<facility>/<instrument>.<filter_name>', e.g. '2MASS/2MASS.J'""".format(band, ', '.join(bands))

                # Try a web query
                if '/' in band:
                    try:
                        self.load_web(band)
                    except:
                        raise IndexError(err)
                
                else:
                    
                    # Or throw an error
                    raise IndexError(err)

        # Set the wavelength and throughput
        self._wave_units = q.AA
        self._wave = np.array([self.raw[0]]) * self.wave_units
        self._throughput = np.array([self.raw[1]])

        # Set n_bins and pixels_per_bin
        self.n_bins = 1
        self.pixels_per_bin = self.raw.shape[-1]

        # Rename some values and apply units
        self.wave_min = self.WavelengthMin * self.wave_units
        self.wave_max = self.WavelengthMax * self.wave_units
        self.wave_eff = self.WavelengthEff * self.wave_units
        self.wave_center = self.WavelengthCen * self.wave_units
        self.wave_mean = self.WavelengthMean * self.wave_units
        self.wave_peak = self.WavelengthPeak * self.wave_units
        self.thru_peak = self.raw[1].max()
        self.wave_phot = self.WavelengthPhot * self.wave_units
        self.wave_pivot = self.WavelengthPivot * self.wave_units
        self.width_eff = self.WidthEff * self.wave_units
        self.fwhm = self.FWHM * self.wave_units
        self.hm_x1 *= self.wave_units
        self.hm_x2 *= self.wave_units
        self.zp = self.ZeroPoint * q.Unit(self.ZeroPointUnit)

        # Delete redundant attributes
        del self.WavelengthMin, self.WavelengthMax, self.WavelengthEff
        del self.WavelengthCen, self.WavelengthMean, self.WavelengthPeak
        del self.WavelengthPhot, self.WavelengthPivot, self.WidthEff, self.FWHM
        del self.ZeroPointUnit, self.ZeroPoint
        try:
            del self.WavelengthUnit
        except AttributeError:
            pass

        # Set the wavelength units
        if wave_units is not None:
            self.wave_units = wave_units

        # Set zeropoint flux units
        if flux_units is not None:
            self._flux_units = self.zp.unit
            self.flux_units = flux_units

        # Get references
        self.refs = []
        try:
            if isinstance(self.CalibrationReference, str):
                self.refs = [self.CalibrationReference.split('=')[-1]]
        except:
            self.CalibrationReference = None

        # Set a base name
        self.name = self.filterID.split('/')[-1]

        # Try to get the extinction vector R from Green et al. (2018)
        self.ext_vector = EXTINCTION.get(self.name, 0)

        # Set the systematic uncertainty (default 2 percent)
        self.systematics = SYSTEMATICS.get(self.name, 0.02)

        # Bin
        if kwargs:
            bwargs = {k: v for k, v in kwargs.items() if k in
                      inspect.signature(self.bin).parameters.keys()}
            self.bin(**bwargs)

    def apply(self, spectrum, plot=False):
        """
        Apply the filter to the given [W, F], or [W, F, E] spectrum

        Parameters
        ----------
        spectrum: array-like
            The wavelength [um] and flux of the spectrum
            to apply the filter to
        plot: bool
            Plot the original and filtered spectrum

        Returns
        -------
        np.ndarray
            The filtered spectrum and error

        """
        # Convert to filter units if possible
        f_units = 1.
        if hasattr(spectrum[0], 'unit'):
            spectrum[0] = spectrum[0].to(self.wave_units)
        if hasattr(spectrum[1], 'unit'):
            spectrum[1] = spectrum[1].to(self.flux_units)
            f_units = self.flux_units
        if len(spectrum) >= 3 and hasattr(spectrum[2], 'unit'):
            spectrum[2] = spectrum[2].to(self.flux_units)

        # Make into iterable arrays
        wav, flx, *err = [np.asarray(i) for i in spectrum]

        # Check for error array
        if len(err) == 0:
            err = np.ones_like(flx) * np.nan
            unc = False
        else:
            err = err[0]
            unc = True

        # Make flux 2D
        if len(flx.shape) == 1:
            flx = np.expand_dims(flx, axis=0)
            err = np.expand_dims(err, axis=0)

        # Make throughput 3D
        rsr = np.copy(self.rsr)

        # Make empty filtered arrays
        filtered_flx = np.zeros((rsr.shape[0], flx.shape[0], rsr.shape[2]))
        filtered_err = np.zeros_like(filtered_flx)

        # Rebin the input spectra to the filter wavelength array
        # and apply the RSR curve to the spectrum
        for i, bn in enumerate(rsr):
            for j, (f, e) in enumerate(zip(flx, err)):
                filtered_flx[i][j] = np.interp(bn[0], wav, f, left=np.nan, right=np.nan) * bn[1]
                filtered_err[i][j] = np.interp(bn[0], wav, e, left=np.nan, right=np.nan) * bn[1]

        # Propagate the filter systematic uncertainties
        if unc:
            filtered_err += filtered_flx*self.systematics

        if plot:

            # Make the figure
            COLORS = color_gen('Category10')
            xlab = 'Wavelength [{}]'.format(self.wave_units)
            ylab = 'Flux Density [{}]'.format(self.flux_units)
            fig = figure(title=self.filterID, x_axis_label=xlab, y_axis_label=ylab)

            # Plot the unfiltered spectrum
            fig.line(wav, flx[0], legend='Input spectrum', color='black')
 
            # Plot the uncertainties
            if unc:
                band_x = np.append(wav, wav[::-1])
                band_y = np.append(flx - err, (flx + err)[::-1])
                fig.patch(band_x, band_y, color='black', fill_alpha=0.1, line_alpha=0)

            # Plot each spectrum bin
            for wav, bn, bne in zip(self.wave, filtered_flx, filtered_err):
                color = next(COLORS)
                fig.line(wav, bn[0], color=color)

                # Plot the uncertainties
                if unc:
                    band_x = np.append(wav, wav[::-1])
                    band_y = np.append(bn[0] - bne[0], (bn[0] + bne[0])[::-1])
                    fig.patch(band_x, band_y, color=color, fill_alpha=0.1, line_alpha=0)

            show(fig)

        return filtered_flx.squeeze()*f_units, filtered_err.squeeze()*f_units

    def bin(self, n_bins=1, pixels_per_bin=None, wave_min=None, wave_max=None):
        """
        Break the filter up into bins and apply a throughput to each bin,
        useful for G141, G102, and other grisms

        Parameters
        ----------
        n_bins: int
            The number of bins to dice the throughput curve into
        pixels_per_bin: int (optional)
            The number of channels per bin, which will be used
            to calculate n_bins
        wave_min: astropy.units.quantity (optional)
            The minimum wavelength to use
        wave_max: astropy.units.quantity (optional)
            The maximum wavelength to use
        """
        # Get wavelength limits
        if wave_min is not None:
            self.wave_min = wave_min
        if wave_max is not None:
            self.wave_max = wave_max

        # Trim the wavelength by the given min and max
        raw_wave = self.raw[0]
        whr = np.logical_and(raw_wave * q.AA >= self.wave_min, raw_wave * q.AA <= self.wave_max)
        self.wave = (raw_wave[whr] * q.AA).to(self.wave_units)
        self.throughput = self.raw[1][whr]
        print('Bandpass trimmed to {} - {}'.format(self.wave_min, self.wave_max))

        # Calculate the number of bins and channels
        pts = len(self.wave)
        if isinstance(pixels_per_bin, int):
            self.pixels_per_bin = pixels_per_bin
            self.n_bins = int(pts/self.pixels_per_bin)
        elif isinstance(n_bins, int):
            self.n_bins = n_bins
            self.pixels_per_bin = int(pts/self.n_bins)
        else:
            raise ValueError("Please specify 'n_bins' OR 'pixels_per_bin' as integers.")

        print('{} bins of {} pixels each.'.format(self.n_bins,
                                                  self.pixels_per_bin))

        # Trim throughput edges so that there are an integer number of bins
        new_len = self.n_bins * self.pixels_per_bin
        start = (pts - new_len) // 2
        self.wave = self.wave[start:new_len+start].reshape(self.n_bins, self.pixels_per_bin)
        self.throughput = self.throughput[start:new_len+start].reshape(self.n_bins, self.pixels_per_bin)

    @property
    def flux_units(self):
        """A getter for the flux units"""
        return self._flux_units

    @flux_units.setter
    def flux_units(self, units):
        """
        A setter for the flux units

        Parameters
        ----------
        units: str, astropy.units.core.PrefixUnit
            The desired units of the zeropoint flux density
        """
        # Check that the units are valid
        dtypes = (q.core.PrefixUnit, q.quantity.Quantity, q.core.CompositeUnit)
        if not isinstance(units, dtypes):
            raise ValueError(units, "units not understood.")

        # Check that the units changed
        if units != self.flux_units:

            # Convert to new units
            sfd = q.spectral_density(self.wave_eff)
            self.zp = self.zp.to(units, equivalencies=sfd)

            # Store new units
            self._flux_units = units

    def info(self, fetch=False):
        """
        Print a table of info about the current filter
        """
        # Get the info from the class
        tp = (int, bytes, bool, str, float, tuple, list, np.ndarray)
        info = [[k, str(v)] for k, v in vars(self).items() if isinstance(v, tp) and k not in ['rsr', 'raw'] and not k.startswith('_')]

        # Make the table
        table = at.Table(np.asarray(info).reshape(len(info), 2), names=['Attributes', 'Values'])

        # Sort and print
        table.sort('Attributes')

        if fetch:
            return table
        else:
            table.pprint(max_width=-1, max_lines=-1, align=['>', '<'])

    def load_TopHat(self, wave_min, wave_max, pixels_per_bin=100):
        """
        Loads a top hat filter given wavelength min and max values

        Parameters
        ----------
        wave_min: astropy.units.quantity (optional)
            The minimum wavelength to use
        wave_max: astropy.units.quantity (optional)
            The maximum wavelength to use
        n_pixels: int
            The number of pixels for the filter
        """
        # Get min, max, effective wavelengths and width
        self.pixels_per_bin = pixels_per_bin
        self.n_bins = 1
        self._wave_units = q.AA
        wave_min = wave_min.to(self.wave_units)
        wave_max = wave_max.to(self.wave_units)

        # Create the RSR curve
        self._wave = np.linspace(wave_min, wave_max, pixels_per_bin)
        self._throughput = np.ones_like(self.wave)
        self.raw = np.array([self.wave.value, self.throughput])

        # Calculate the effective wavelength
        wave_eff = ((wave_min + wave_max) / 2.).value
        width = (wave_max - wave_min).value

        # Add the attributes
        self.path = ''
        self.refs = ''
        self.Band = 'Top Hat'
        self.CalibrationReference = ''
        self.FWHM = width
        self.hm_x1 = wave_min.value
        self.hm_x2 = wave_max.value
        self.Facility = '-'
        self.FilterProfileService = '-'
        self.MagSys = '-'
        self.PhotCalID = ''
        self.PhotSystem = ''
        self.ProfileReference = ''
        self.WavelengthMin = wave_min.value
        self.WavelengthMax = wave_max.value
        self.WavelengthCen = wave_eff
        self.WavelengthEff = wave_eff
        self.WavelengthMean = wave_eff
        self.WavelengthPeak = wave_eff
        self.WavelengthPhot = wave_eff
        self.WavelengthPivot = wave_eff
        self.WavelengthUCD = ''
        self.WidthEff = width
        self.ZeroPoint = 0
        self.ZeroPointType = ''
        self.ZeroPointUnit = 'Jy'
        self.filterID = 'Top Hat'

    def load_txt(self, filepath):
        """Load the filter from a txt file

        Parameters
        ----------
        file: str
            The filepath
        """
        # Get the transmission data from the file
        self.raw = np.genfromtxt(filepath, unpack=True)

        # Convert to Angstroms if microns
        if self.raw[0][-1] < 100:
            self.raw[0] = self.raw[0] * 10000

        # Load it into the object
        self.filterID = os.path.splitext(os.path.basename(filepath))[0]
        self.load_raw()

    def load_raw(self, wave_units='AA'):
        """
        The raw data to calculate the filter properties from

        Parameters
        ----------
        data: sequence
            The wavelength and throughput for the filter
        """
        funit = q.erg / q.s / q.cm**2 / q.AA
        self.WavelengthUnit = wave_units
        self.ZeroPointUnit = str(funit)
        x, f = self.raw

        # Rebin Vega to filter
        vega_file = resource_filename('svo_filters', 'data/spectra/vega.txt')
        vega_data = np.genfromtxt(vega_file, unpack=True)[: 2]
        vega_data[0] *= 10000
        vega = rebin_spec(vega_data, x)

        # Calculate the filter's properties
        self.ZeroPoint = np.trapz(f * x * vega, x=x) / np.trapz(f * x, x=x)
        self.WavelengthPeak = x[np.argmax(f)]
        self.WavelengthMin = x[np.where(f > f.max() / 100.)][0]
        self.WavelengthMax = x[np.where(f > f.max() / 100.)][-1]
        self.WavelengthEff = np.trapz(f * x**2 * vega, x=x) / np.trapz(f * x * vega, x=x)
        self.WavelengthMean = np.trapz(f * x, x=x) / np.trapz(f, x=x)
        self.WidthEff = np.trapz(f, x=x) / f.max()
        self.WavelengthPivot = np.sqrt(np.trapz(f, x=x) / np.trapz(f / x**2, x=x))
        self.WavelengthPhot = np.trapz(f * vega * x**3, x=x) / np.trapz(f * vega * x**2, x=x)

        # Half max stuff
        halfmax = f.max() / 2.
        self.hm_x1 = x[f > halfmax][0]
        self.hm_x2 = x[f > halfmax][-1]
        self.FWHM = self.hm_x2 - self.hm_x1
        self.WavelengthCen = (self.hm_x1 + self.hm_x2) / 2.

        # Add missing attributes
        self.path = ''
        self.pixels_per_bin = self.raw.shape[-1]
        self.n_bins = 1

    def load_web(self, filt):
        """Load the filter from a Web query

        Parameters
        ----------
        filt: str
            The '<facility>/<instrument>.<filter_name>' of the filter, e.g. '2MASS/2MASS.J'
        """
        # Get the transmission data for the filter
        data = SvoFps.get_transmission_data(filt)

        # Make into arrays
        self.raw = np.array([np.array(data['Wavelength']), np.array(data['Transmission'])])

        # Load it into the object
        self.filterID = filt
        self.load_raw(wave_units=str(data['Wavelength'].unit))

    def load_xml(self, filepath):
        """Load the filter from a txt file

        Parameters
        ----------
        filepath: str
            The filepath for the filter
        """
        # Parse the XML file
        vot = vo.parse_single_table(filepath)
        self.raw = np.array([list(i) for i in vot.array]).T

        # Parse the filter metadata
        for p in [str(p).split() for p in vot.params]:

            # Extract the key/value pairs
            key = p[1].split('"')[1]
            val = p[-1].split('"')[1]

            # Do some formatting
            flt1 = p[2].split('"')[1] == 'float'
            flt2 = p[3].split('"')[1] == 'float'
            if flt1 or flt2:
                val = float(val)

            else:
                val = val.replace('b&apos;', '').replace('&apos', '').replace('&amp;', '&').strip(';')

            # Set the attribute
            if key != 'Description':
                setattr(self, key, val)

        halfmax = self.raw[1].max() / 2.
        self.hm_x1 = self.raw[0][self.raw[1] > halfmax][0]
        self.hm_x2 = self.raw[0][self.raw[1] > halfmax][-1]

        # Create some attributes
        self.path = filepath
        self.pixels_per_bin = self.raw.shape[-1]
        self.n_bins = 1

    def overlap(self, spectrum):
        """Tests for overlap of this filter with a spectrum

        Example of full overlap:

            |---------- spectrum ----------|
               |------ self ------|

        Examples of partial overlap: :

            |---------- self ----------|
               |------ spectrum ------|

            |---- spectrum ----|
               |----- self -----|

            |---- self ----|
               |---- spectrum ----|

        Examples of no overlap: :

            |---- spectrum ----|  |---- other ----|

            |---- other ----|  |---- spectrum ----|

        Parameters
        ----------
        spectrum: sequence
            The [W, F] spectrum with astropy units

        Returns
        -------
        ans : {'full', 'partial', 'none'}
            Overlap status.
        """
        swave = self.wave[np.where(self.throughput != 0)]
        s1, s2 = swave.min(), swave.max()

        owave = spectrum[0]
        o1, o2 = owave.min(), owave.max()

        if (s1 >= o1 and s2 <= o2):
            ans = 'full'

        elif (s2 < o1) or (o2 < s1):
            ans = 'none'

        else:
            ans = 'partial'

        return ans

    def plot(self, details=False, fig=None, draw=True):
        """
        Plot the filter

        Parameters
        ----------
        details: bool
            Plot the filter details
        fig: bokeh.plotting.figure (optional)
            A figure to plot on
        draw: bool
            Draw the figure, else return it

        Returns
        -------
        bokeh.plotting.figure
            The filter figure
        """
        COLORS = color_gen('Category10')

        # Make the figure
        if fig is None:
            xlab = 'Wavelength [{}]'.format(self.wave_units)
            ylab = 'Throughput'
            title = self.filterID
            fig = figure(title=title, x_axis_label=xlab, y_axis_label=ylab)

        # Plot the raw curve
        fig.line((self.raw[0] * q.AA).to(self.wave_units), self.raw[1], alpha=0.1, line_width=8, color='black')

        # Plot each bin
        for x, y in self.rsr:
            fig.line(x, y, color=next(COLORS), line_width=2)

        # Plot details
        if details:

            dcolor = 'red'

            # Min and Max
            fig.line([self.wave_min] * 2, [0, self.thru_peak], color=dcolor, line_dash='dashed', legend_label='wave_min')
            fig.line([self.wave_max] * 2, [0, self.thru_peak], color=dcolor, line_dash='dashed', legend_label='wave_max')

            # FWHM
            fig.line([self.hm_x1, self.hm_x2], [self.thru_peak / 2.] * 2, color=dcolor, line_dash='dotted', legend_label='fwhm')

            # Max throughput
            fig.circle([self.wave_peak.value], [self.thru_peak], fill_color=dcolor, line_color=dcolor, size=8, legend_label='max_thru')

            # Effective wavelength
            fig.line([self.wave_eff] * 2, [0, self.thru_peak], color=dcolor, line_dash='solid', legend_label='wave_eff')

            # Click policy
            fig.legend.click_policy = 'hide'

        if draw:
            show(fig)
        else:
            return fig

    @property
    def rsr(self):
        """A getter for the relative spectral response (rsr) curve"""
        arr = np.array([self.wave.value, self.throughput]).swapaxes(0, 1)

        return arr

    @property
    def throughput(self):
        """A getter for the throughput"""
        return self._throughput

    @throughput.setter
    def throughput(self, points):
        """A setter for the throughput

        Parameters
        ----------
        throughput: sequence
            The array of throughput points
        """
        # Test shape
        if not points.shape == self.wave.shape:
            raise ValueError("Throughput and wavelength must be same shape.")

        self._throughput = points

    @property
    def wave(self):
        """A getter for the wavelength"""
        return self._wave

    @wave.setter
    def wave(self, wavelength):
        """A setter for the wavelength

        Parameters
        ----------
        wavelength: astropy.units.quantity.Quantity
            The array with units
        """
        # Test units
        if not isinstance(wavelength, q.quantity.Quantity):
            raise ValueError("Wavelength must be in length units.")

        self._wave = wavelength
        self.wave_units = wavelength.unit

    @property
    def wave_units(self):
        """A getter for the wavelength units"""
        return self._wave_units

    @wave_units.setter
    def wave_units(self, units):
        """
        A setter for the wavelength units

        Parameters
        ----------
        units: str, astropy.units.core.PrefixUnit
            The wavelength units
        """
        # Make sure it's length units
        if not units.is_equivalent(q.m):
            raise ValueError(units, ": New wavelength units must be a length.")

        # Update the units
        self._wave_units = units

        # Update all the wavelength values
        self._wave = self.wave.to(self.wave_units).round(5)
        self.wave_min = self.wave_min.to(self.wave_units).round(5)
        self.wave_max = self.wave_max.to(self.wave_units).round(5)
        self.wave_eff = self.wave_eff.to(self.wave_units).round(5)
        self.wave_center = self.wave_center.to(self.wave_units).round(5)
        self.wave_mean = self.wave_mean.to(self.wave_units).round(5)
        self.wave_peak = self.wave_peak.to(self.wave_units).round(5)
        self.wave_phot = self.wave_phot.to(self.wave_units).round(5)
        self.wave_pivot = self.wave_pivot.to(self.wave_units).round(5)
        self.width_eff = self.width_eff.to(self.wave_units).round(5)
        self.fwhm = self.fwhm.to(self.wave_units).round(5)
        self.hm_x1 = self.hm_x1.to(self.wave_units).round(5)
        self.hm_x2 = self.hm_x2.to(self.wave_units).round(5)


def color_gen(colormap='viridis', key=None, n=15):
    """Color generator for Bokeh plots

    Parameters
    ----------
    colormap: str, sequence
        The name of the color map

    Returns
    -------
    generator
        A generator for the color palette
    """
    if colormap in dir(bpal):
        palette = getattr(bpal, colormap)

        if isinstance(palette, dict):
            if key is None:
                key = list(palette.keys())[0]
            palette = palette[key]

        elif callable(palette):
            palette = palette(n)

        else:
            raise TypeError("pallette must be a bokeh palette name or a sequence of color hex values.")

    elif isinstance(colormap, (list, tuple)):
        palette = colormap

    else:
        raise TypeError("pallette must be a bokeh palette name or a sequence of color hex values.")

    yield from itertools.cycle(palette)


def filters(filter_directory=None):
    """
    Get a list of the available filters

    Parameters
    ----------
    filter_directory: str
        The directory containing the filter relative spectral response curves
    update: bool
        Check the filter directory for new filters and generate pickle of table
    fmt: str
        The format for the returned table

    Returns
    -------
    list
        The list of band names
    """
    if filter_directory is None:
        filter_directory = resource_filename('svo_filters', 'data/filters/')

    # Get list of files from dir
    files = glob(os.path.join(filter_directory, '*'))

    # Grab the basename as the filter name
    bands = [b.replace(filter_directory, '').replace('.txt', '').replace('.xml', '') for b in files]

    return bands


def rebin_spec(spec, wavnew, oversamp=100, plot=False):
    """
    Rebin a spectrum to a new wavelength array while preserving
    the total flux

    Parameters
    ----------
    spec: array-like
        The wavelength and flux to be binned
    wavenew: array-like
        The new wavelength array

    Returns
    -------
    np.ndarray
        The rebinned flux

    """
    wave, flux = spec
    nlam = len(wave)
    x0 = np.arange(nlam, dtype=float)
    x0int = np.arange((nlam-1.) * oversamp + 1., dtype=float)/oversamp
    w0int = np.interp(x0int, x0, wave)
    spec0int = np.interp(w0int, wave, flux)/oversamp

    # Set up the bin edges for down-binning
    maxdiffw1 = np.diff(wavnew).max()
    w1bins = np.concatenate(([wavnew[0]-maxdiffw1],
                             .5*(wavnew[1::]+wavnew[0: -1]),
                             [wavnew[-1]+maxdiffw1]))

    # Bin down the interpolated spectrum:
    w1bins = np.sort(w1bins)
    nbins = len(w1bins)-1
    specnew = np.zeros(nbins)
    inds2 = [[w0int.searchsorted(w1bins[ii], side='left'),
              w0int.searchsorted(w1bins[ii+1], side='left')]
             for ii in range(nbins)]

    for ii in range(nbins):
        specnew[ii] = np.sum(spec0int[inds2[ii][0]: inds2[ii][1]])

    return specnew
