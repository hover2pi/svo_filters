"""Some tests to make sure the filters work as intended"""
import unittest

import astropy.units as q

from svo_filters import svo

 
class TestFilters(unittest.TestCase):
    """Tests for Filter class"""
    def setUp(self):
        pass

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
