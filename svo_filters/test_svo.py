"""Some tests to make sure the filters work as intended"""
import unittest
import astropy.units as q
import svo

 
class TestFilters(unittest.TestCase):
    """Tests for Filter class"""
    def setUp(self):
        pass

    def test_filter_class(self):
        """Test if Filter object is created properly"""
        filt = svo.Filter('2MASS.J')

        is_filter = isinstance(filt, svo.Filter)

        self.assertTrue(is_filter)

    def test_tophat(self):
        """Test if Filter object is created properly"""
        filt = svo.Filter('tophat', wl_min=0.8*q.um, wl_max=1.2*q.um)

        is_filter = isinstance(filt, svo.Filter)

        self.assertTrue(is_filter)

if __name__ == '__main__':
    unittest.main()