#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A Python wrapper for the SVO Filter Profile Service
"""
from __future__ import print_function
from astropy.utils.exceptions import AstropyWarning
from glob import glob
import astropy.table as at
import astropy.io.votable as vo
import astropy.units as q
import matplotlib.pyplot as plt
import warnings
import pkg_resources
import numpy as np
import urllib
import os

warnings.simplefilter('ignore', category=AstropyWarning)
WL_KEYS = ['FWHM', 'WavelengthCen', 'WavelengthEff', 'WavelengthMax',
           'WavelengthMean', 'WavelengthMin', 'WavelengthPeak', 
           'WavelengthPhot', 'WavelengthPivot', 'WidthEff']

class Filter(object):
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
    FWHM:float
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
    def __init__(self, band, filter_directory=pkg_resources.resource_filename('svo_filters', 'data/filters/'), 
                 wl_units=q.um, **kwargs):
        """
        Loads the bandpass data into the Filter object
        
        Parameters
        ----------
        band: str
            The bandpass filename (e.g. 2MASS.J)
        filter_directory: str
            The directory containing the filter files
        wl_units: str, astropy.units.core.PrefixUnit, astropy.units.core.CompositeUnit (optional)
            The wavelength units
        """
        # Get list of filters
        filters = filter_list(filter_directory)
        filepath = filters['path']+band
        
        # If the filter is missing, ask what to do
        if filepath not in filters['files']:
            
            print('Current filters:',
                  ', '.join(filters['bands']),
                  '\n')
        
            print('No filters match',filepath)
            dl = input('Would you like me to download it? [y/n] ')
            
            if dl.lower()=='y':
                
                # Prompt for new filter
                print('\nA full list of available filters from the\n'\
                      'SVO Filter Profile Service can be found at\n'\
                      'http://svo2.cab.inta-csic.es/theory/fps3/\n')
                band = input('Enter the band name to retrieve (e.g. 2MASS/2MASS.J): ')
                
                # Download the XML (VOTable) file
                baseURL = 'http://svo2.cab.inta-csic.es/svo/theory/fps/fps.php?ID='
                filepath = filter_directory+os.path.basename(band)
                _ = urllib.request.urlretrieve(baseURL+band, filepath)
                
                # Print the new filepath
                print('Band stored as',filepath)
            
            else:
                return
            
        # Try to read filter info
        try:
            
            # Parse the XML file
            vot = vo.parse_single_table(filepath)
            self.rsr = np.array([list(i) for i in vot.array]).T
            
            # Parse the filter metadata
            for p in [str(p).split() for p in vot.params]:
                
                # Extract the key/value pairs
                key = p[1].split('"')[1]
                val = p[-1].split('"')[1]
                
                # Do some formatting
                if p[2].split('"')[1]=='float'\
                or p[3].split('"')[1]=='float':
                    val = float(val)
                
                else:
                    val = val.replace('b&apos;','')\
                             .replace('&apos','')\
                             .replace('&amp;','&')\
                             .strip(';')
                
                # Set the attribute
                if key!='Description':
                    setattr(self, key, val)
                    
            # Set wavelength units
            self.WavelengthUnit = q.Unit(self.WavelengthUnit)
            for key in WL_KEYS:
                setattr(self, key, getattr(self, key)*self.WavelengthUnit)
                
            # Create some attributes
            self.path = filepath
            self.n_channels = len(self.rsr[0])
            self.n_bins = 1
            self.raw = self.rsr
            
            # Get the bin centers
            w_cen = np.nanmean(self.rsr[0])
            f_cen = np.nanmean(self.rsr[1])
            self.centers = np.asarray([[w_cen],[f_cen]])
            
            try:
                self.refs = [self.CalibrationReference.split('=')[-1]]
            except:
                self.refs = []
                
            # Bin from the get-go
            if kwargs:
                self.bin(**kwargs)
                
            # Set the wavelength units
            if wl_units:
                self.set_units(wl_units)
            
        # If empty, delete XML file
        except TypeError:
            
            print('No filter named',band)
            # if os.path.isfile(filepath):
            #     os.remove(filepath)
                
            return
            
    def apply(self, spectrum, plot=False):
        """
        Apply the filter to the given spectrum
        
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
            The filtered spectrum
            
        """
        # Make into iterable arrays
        wav, flx = [np.asarray(i) for i in spectrum]
        
        # Make flux 2D
        if len(flx.shape)==1:
            flx = np.expand_dims(flx, axis=0)
        
        # Make throughput 3D
        rsr = np.copy(self.rsr)
        if len(rsr.shape)==2:
            rsr = np.expand_dims(rsr, axis=0)
        
        # Make empty filtered array
        filtered = np.zeros((rsr.shape[0],flx.shape[0],rsr.shape[2]))
        
        # Rebin the input spectra to the filter wavelength array
        # and apply the RSR curve to the spectrum
        for i,bn in enumerate(rsr):
            for j,f in enumerate(flx):
                filtered[i][j] = np.interp(bn[0], wav, f)*bn[1]
                
        if plot:
            plt.loglog(wav, flx[0])
            for n,bn in enumerate(rsr):
                plt.loglog(bn[0], filtered[n][0])
        
        del rsr, wav, flx
        
        return filtered.squeeze()
        
    def bin(self, n_bins='', n_channels='', bin_throughput=''):
        """
        Break the filter up into bins and apply a throughput to each bin,
        useful for G141, G102, and other grisms
        
        Parameters
        ----------
        n_bins: int
            The number of bins to dice the throughput curve into
        n_cahnnels: int (optional)
            The number of channels per bin, which will be used to calculate n_bins
        bin_throughput: array-like (optional)
            The throughput for each bin (top hat by default)
            must be of length n_channels
        """
        # Calculate the number of bins and channels
        rsr = len(self.raw[0])
        if n_channels and isinstance(n_channels,int):
            self.n_channels = int(n_channels)
            self.n_bins = int(rsr/self.n_channels)
        elif n_bins and isinstance(n_bins,int):
            self.n_bins = int(n_bins)
            self.n_channels = int(rsr/self.n_bins)
        elif not n_bins and not n_channels \
        and isinstance(bin_throughput, (list,tuple,np.ndarray)):
            pass
        else:
            print('Please specify n_bins or n_channels as integers.')
            return
            
        print('{} bins of {} channels each.'.format(self.n_bins,self.n_channels))
        
        # Trim throughput edges so that there are an integer number of bins
        new_len = self.n_bins*self.n_channels
        start = (rsr-new_len)//2
        self.rsr = np.copy(self.raw[:,start:new_len+start])
        
        # Reshape the throughput array
        self.rsr = self.rsr.reshape(2,self.n_bins,self.n_channels)
        self.rsr = self.rsr.swapaxes(0,1)
        
        # Get the bin centers
        w_cen = np.nanmean(self.rsr[:,0,:], axis=1)
        f_cen = np.nanmean(self.rsr[:,1,:], axis=1)
        self.centers = np.asarray([w_cen,f_cen])
        
        # Get the bin throughput function
        if not isinstance(bin_throughput, (list,tuple,np.ndarray)):
            bin_throughput = np.ones(self.n_channels)
            
        # Make sure the shape is right
        if len(bin_throughput)==self.n_channels:
            
            # Save the attribute
            self.bin_throughput = np.asarray(bin_throughput)
            
            # Apply the bin throughput
            self.rsr[:,1] *= self.bin_throughput
            
        else:
            print('bin_throughput must be an array of length',self.n_channels)
            print('Using top hat throughput for each bin.')
                
    def plot(self):
        """
        Plot the filter
        """
        # If the filter is binned, plot each with bin centers
        try:
            for x,y in self.rsr:
                plt.plot(x, y)
            plt.plot(*self.centers, ls='None', marker='.', c='k')
            plt.plot(self.raw[0], self.raw[1], lw=6, alpha=0.1, zorder=0)
            
        # Otherwise just plot curve
        except:
            plt.plot(*self.rsr)
            
        plt.xlabel('Wavelength [{}]'.format(str(self.WavelengthUnit)))
        plt.ylabel('Throughput')
        
    def info(self):
        """
        Print a table of info about the current filter
        """
        # Get the info from the class 
        tp = (int, bytes, bool, str, float, tuple, list, np.ndarray)
        exclude = ['rsr', 'bin_throughput', 'raw']
        info = [[k,str(v)] for k,v in vars(self).items() if isinstance(v, tp)
                and k not in exclude]
                
        # Make the table
        table = at.Table(np.asarray(info).reshape(len(info),2),
                 names=['Attributes','Values'])
        
        # Sort and print
        table.sort('Attributes')
        table.pprint(max_width=-1, align=['>','<'])
        
    def set_units(self, wl_units=q.um):
        """
        Set the wavelength and flux units
        
        Parameters
        ----------
        wl_units: str, astropy.units.core.PrefixUnit, astropy.units.core.CompositeUnit
            The wavelength units
        """
        # Set wavelength units
        old_unit = self.WavelengthUnit
        self.WavelengthUnit = q.Unit(wl_units)
        for key in WL_KEYS:
            setattr(self, key, getattr(self, key).to(self.WavelengthUnit))
            
        # Update the rsr curve
        const = (old_unit/self.WavelengthUnit).decompose()._scale 
        self.raw[0] *= const
        self.rsr[:,0] *= const
        self.centers[0] *= const
        
def filter_list(filter_directory=pkg_resources.resource_filename('svo_filters', 'data/filters/')):
    """
    Get a list of the available filters
    
    Parameters
    ----------
    filter_directory: str
        The directory containing the filter relative spectral response curves
    
    Returns
    -------
    list
        The list of band names
    """
    files = glob(filter_directory+'*')
    bands = [os.path.basename(b) for b in files]
    
    return {'files':files, 'bands':bands, 'path':filter_directory}
