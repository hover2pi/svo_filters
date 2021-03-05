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

    def test_no_filter(self):
        """Test if exception is raised for bogus filter"""
        self.assertRaises(IndexError, svo.Filter, 'BAD_FILTER')

    def test_filter_web(self):
        """Test if Filter object is created from web query"""
        filt = svo.Filter('Generic/Johnson.B')

        self.assertTrue(isinstance(filt, svo.Filter))

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
        filtered, err_filtered = filt.apply(spec)
        self.assertEqual(filtered.shape, filt.wave.squeeze().shape)

        # Check that the propagated errors are all nans since none were given
        self.assertTrue(np.all([np.isnan(i) for i in err_filtered]))

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
        filtered, err_filtered = filt.apply(spec, plot=True)
        self.assertEqual(filtered.shape, filt.wave.squeeze().shape)

    def test_filter_apply_binned_err(self):
        """Test that the filter gets applied to a spectrum with errors properly"""
        filt = svo.Filter('2MASS.J', n_bins=4)
        funit = q.erg/q.s/q.cm**2/q.AA
        spec = [np.linspace(0.9, 2, 1000)*q.um,
                np.ones(1000)*funit,
                np.ones(1000)*0.05*funit]

        # Apply the filter to the spectrum
        filtered, err_filtered = filt.apply(spec)

        # Check that the units are still there
        self.assertEqual(filtered.unit,funit)
        self.assertEqual(err_filtered.shape, filt.wave.squeeze().shape)

    def test_filter_apply_no_units(self):
        """Test that the apply method works with and without units"""
        filt = svo.Filter('2MASS.J')
        spec = [np.linspace(0.9, 2, 1000),
                np.ones(1000)]

        # Apply the filter to the spectrum
        filtered, err_filtered = filt.apply(spec)
        self.assertFalse(hasattr(filtered, 'unit'))
        self.assertTrue(np.all([np.isnan(i) for i in err_filtered]))

    def test_plot(self):
        """Test that the plots are produced properly"""
        filt = svo.Filter('2MASS.J')
        filt.plot(details=True)
        plt = filt.plot(draw=False)

        self.assertTrue(type(plt) == Figure)

class TestFilterList(unittest.TestCase):
    """Tests for filter function"""
    def setUp(self):
        pass

    def test_default(self):
        """Test default directory"""
        filts = svo.filters()
        self.assertTrue(len(filts) > 0)
