"""Some tests to make sure the filters work as intended"""

import astropy.units as q
from . import Filter

def test_filter_class():
    """Test if Filter object is created properly"""
    filt = Filter('2MASS.J')
    
    is_filter = isinstance(filt, Filter)
    
    assert is_filter
    
def test_tophat():
    """Test if Filter object is created properly"""
    filt = Filter('tophat', wl_min=0.8*q.um, wl_max=1.2*q.um)
    
    is_filter = isinstance(filt, Filter)
    
    assert is_filter