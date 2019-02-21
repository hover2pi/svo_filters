"""Some tests to make sure the filters work as intended"""
import unittest

import astropy.units as q
import astropy.table as at
from bokeh.plotting.figure import Figure
import numpy as np

from svo_filters import svo

 
class TestFilter(unittest.TestCase):
    """Tests for Filter class"""
    def setUp(self):
        pass

    def test_info(self):
        """Test that the info attr works"""
        filt = svo.Filter('2MASS.J')
        self.assertTrue(filt.info() == None)
        self.assertTrue(isinstance(filt.info(fetch=True), at.Table))
        self.assertRaises(ValueError, setattr, filt, 'wave', np.arange(10))
        self.assertRaises(ValueError, setattr, filt, 'throughput', np.arange(1))

    def test_centers(self):
        """Test that the centers are returned"""
        filt = svo.Filter('2MASS.J')
        self.assertEqual(filt.centers.shape, (2, 1))

    def test_no_filter(self):
        """Test if exception is raised for bogus filter"""
        self.assertRaises(IOError, svo.Filter, 'BAD_FILTER')

    def test_filter_xml(self):
        """Test if Filter object is created properly"""
        filt = svo.Filter('2MASS.J')

        self.assertTrue(isinstance(filt, svo.Filter))

    def test_filter_txt(self):
        """Test if Filter object is created properly"""
        filt = svo.Filter('NIRISS.GR700XD.1')

        self.assertTrue(isinstance(filt, svo.Filter))

    def test_filter_tophat(self):
        """Test if Filter object is created properly"""
        # Good tophat
        filt = svo.Filter('tophat', wave_min=0.8*q.um, wave_max=1.2*q.um)
        self.assertTrue(isinstance(filt, svo.Filter))

        # Bad tophat
        self.assertRaises(ValueError, svo.Filter, 'tophat')

    def test_filter_bin(self):
        """Test that the binning works"""
        # Test n_bins
        filt = svo.Filter('WFC3_IR.G141', n_bins=10)
        self.assertTrue(filt.wave.shape[0] == 10)

        # Test pixels_per_bin
        filt = svo.Filter('WFC3_IR.G141', pixels_per_bin=50)
        self.assertTrue(filt.wave.shape[1] == 50)

        # Test neither throws an error
        self.assertRaises(ValueError, svo.Filter, 'WFC3_IR.G141', n_bins=None,
                          pixels_per_bin=None)

    def test_filter_units_bad(self):
        """Test that changing the wave units to non-length raises exception"""
        filt = svo.Filter('WFC3_IR.G141')

        # Fun syntax to test attribute setting
        self.assertRaises(ValueError, setattr, filt, 'wave_units', q.second)
        self.assertRaises(ValueError, setattr, filt, 'flux_units', q.second)

    def test_spectrum(self):
        """Test that the filter gets applied to a spectrum properly"""
        filt = svo.Filter('2MASS.J')
        spec = [np.linspace(0.9, 2, 1000)*q.um,
                np.ones(1000)*q.erg/q.s/q.cm**2/q.AA]

        # Apply the filter to the spectrum
        filtered = filt.apply(spec)
        self.assertEqual(filtered.shape, filt.wave.squeeze().shape)

        # Check overlap
        self.assertEqual(svo.Filter('2MASS.J').overlap(spec), 'full')
        self.assertEqual(svo.Filter('2MASS.Ks').overlap(spec), 'partial')
        self.assertEqual(svo.Filter('WISE.W4').overlap(spec), 'none')

    def test_filter_apply_binned(self):
        """Test that the filter gets applied to a spectrum properly"""
        filt = svo.Filter('2MASS.J', n_bins=4)
        spec = [np.linspace(0.9, 2, 1000)*q.um,
                np.ones(1000)*q.erg/q.s/q.cm**2/q.AA]

        # Apply the filter to the spectrum
        filtered = filt.apply(spec, plot=True)
        self.assertEqual(filtered.shape, filt.wave.squeeze().shape)

    def test_filter_apply_binned_err(self):
        """Test that the filter gets applied to a spectrum with errors properly"""
        filt = svo.Filter('2MASS.J', n_bins=4)
        spec = [np.linspace(0.9, 2, 1000)*q.um,
                np.ones(1000)*q.erg/q.s/q.cm**2/q.AA,
                np.ones(1000)*0.05*q.erg/q.s/q.cm**2/q.AA]

        # Apply the filter to the spectrum
        filtered, err_filtered = filt.apply(spec)
        self.assertEqual(err_filtered.shape, filt.wave.squeeze().shape)

    def test_plot(self):
        """Test that the plots are produced properly"""
        filt = svo.Filter('2MASS.J')
        filt.plot()
        plt = filt.plot(draw=False)

        self.assertTrue(type(plt) == Figure)

class TestFilterList(unittest.TestCase):
    """Tests for filter function"""
    def setUp(self):
        self.filts = svo.filters(data=True)
        self.filts_dict = svo.filters(data=True, fmt='dict')

    def test_dir(self):
        """Test new directory"""
        self.assertRaises(IOError, svo.filters, filter_directory='foobar')

    def test_update(self):
        """Test that the update works"""
        self.assertRaises(IOError, svo.filters, update=True)

    def test_filter_dtypes(self):
        """Test if a table of filters is returned"""
        self.assertTrue(isinstance(self.filts, at.Table))
        self.assertTrue(isinstance(self.filts_dict, dict))

    def test_filter_len(self):
        """Test that the table isn't empty"""
        self.assertTrue(len(self.filts) > 0)
