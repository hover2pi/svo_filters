"""Some tests to make sure the filters work as intended"""
import unittest

import astropy.units as q
import astropy.table as at

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
