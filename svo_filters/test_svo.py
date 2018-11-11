"""Some tests to make sure the filters work as intended"""
import unittest

import astropy.units as q
import astropy.table as at
import numpy as np

from svo_filters import svo

 
class TestFilter(unittest.TestCase):
    """Tests for Filter class"""
    def setUp(self):
        pass

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
        filt = svo.Filter('tophat', wave_min=0.8*q.um, wave_max=1.2*q.um)

        self.assertTrue(isinstance(filt, svo.Filter))

    def test_filter_bin(self):
        """Test that the binning works"""
        filt = svo.Filter('WFC3_IR.G141', n_bins=10)

        self.assertTrue(filt.wave.shape[0] == 10)

    def test_filter_units_bad(self):
        """Test that changing the wave units to non-length raises exception"""
        filt = svo.Filter('WFC3_IR.G141')

        # Fun syntax to test attribute setting
        self.assertRaises(ValueError, setattr, filt, 'wave_units', q.second)

    def test_filter_apply(self):
        """Test that the filter gets applied to a spectrum properly"""
        filt = svo.Filter('2MASS.J')
        spec = [np.linspace(0.9, 2, 1000)*q.um,
                np.ones(1000)*q.erg/q.s/q.cm**2/q.AA]

        # Apply the filter to the spectrum
        filtered = filt.apply(spec)
        self.assertEqual(filtered.shape, filt.wave.squeeze().shape)

    def test_filter_apply_binned(self):
        """Test that the filter gets applied to a spectrum properly"""
        filt = svo.Filter('2MASS.J', n_bins=4)
        spec = [np.linspace(0.9, 2, 1000)*q.um,
                np.ones(1000)*q.erg/q.s/q.cm**2/q.AA]

        # Apply the filter to the spectrum
        filtered = filt.apply(spec)
        self.assertEqual(filtered.shape, filt.wave.squeeze().shape)


class TestFilterList(unittest.TestCase):
    """Tests for filter function"""
    def setUp(self):
        self. filts = svo.filters(data=True)

    def test_filter_table(self):
        """Test if a table of filters is returned"""
        self.assertTrue(isinstance(self.filts, at.Table))

    def test_filter_len(self):
        """Test that the table isn't empty"""
        self.assertTrue(len(self.filts) > 0)
