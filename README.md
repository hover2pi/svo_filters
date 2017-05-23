# Introduction

The [Spanish Virtual Observatory (SVO) Filter Profile Service](http://svo2.cab.inta-csic.es/theory/fps3/) is a great resource for homogenized photometric filter curves and metadata. With `svo_filters`, I tried to create a lightweight and flexible package to incorporate these filters into Python applications.

Install in the usual fashion with

```
cd svo_filters
python setup.py install
```

Then import like so


```python
from svo_filters import svo
```

The actual filters are stored locally as XML files and can be viewed with


```python
svo.filter_list()
```




    {'bands': ['2MASS.H',
      '2MASS.J',
      '2MASS.Ks',
      'IRAC.I1',
      'IRAC.I2',
      'IRAC.I3',
      'IRAC.I4',
      'Kepler.K',
      'WFC3_IR.G102',
      'WFC3_IR.G141',
      'WFC3_UVIS2.F814W',
      'WFC3_UVIS2.F850LP',
      'WISE.W1',
      'WISE.W2',
      'WISE.W3',
      'WISE.W4'],
     'files': ['/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/2MASS.H',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/2MASS.J',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/2MASS.Ks',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/IRAC.I1',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/IRAC.I2',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/IRAC.I3',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/IRAC.I4',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/Kepler.K',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/WFC3_IR.G102',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/WFC3_IR.G141',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/WFC3_UVIS2.F814W',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/WFC3_UVIS2.F850LP',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/WISE.W1',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/WISE.W2',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/WISE.W3',
      '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/WISE.W4'],
     'path': '/Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/'}



# Load a Photometric Filter

To load a filter object, just pass the band name to the `Filter` class. We can then view the metadata and see a plot!


```python
H_band = svo.Filter('2MASS.H')
H_band.info()
H_band.plot()
```

         Attributes                                              Values                                        
    -------------------- --------------------------------------------------------------------------------------
                    Band H                                                                                     
    CalibrationReference http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2003AJ....126.1090C&db_key=AST
                    FWHM 0.260964753837 um                                                                     
                Facility 2MASS                                                                                 
    FilterProfileService ivo://svo/fps                                                                         
                  MagSys Vega                                                                                  
               PhotCalID 2MASS/2MASS.H/Vega                                                                    
              PhotSystem 2MASS                                                                                 
        ProfileReference http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html                    
           WavelengthCen 1.6487192828100004 um                                                                 
                     ... ...                                                                                   
           WavelengthUCD em.wl                                                                                 
                WidthEff 0.250940236716 um                                                                     
               ZeroPoint 1024.0                                                                                
           ZeroPointType Pogson                                                                                
           ZeroPointUnit Jy                                                                                    
                 centers [[ 1.63794816]
     [ 0.53618449]]                                                        
                filterID 2MASS/2MASS.H                                                                         
                  n_bins 1                                                                                     
              n_channels 58                                                                                    
                    path /Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/2MASS.H          
                    refs ['AST']                                                                               
    Length = 28 rows



![png](svo_demo_files/svo_demo_7_1.png)


The actual relative spectral response curve is stored as an array of the wavelength and throughput.


```python
H_band.rsr.shape
```




    (2, 58)



# Load a Grism

Filters can also be binned arbitrarily, for use with grisms. We can pass integers to the `n_bins` or `n_channels` arguments to specify the number of wavelength bins or channels per bin, respectively.


```python
G141 = svo.Filter('WFC3_IR.G141', n_bins=15)
G141.info()
G141.plot()
```

    15 bins of 634 channels each.
         Attributes                                                                                                                                                                                                                                                                           Values                                                                                                                                                                                                                                                                     
    -------------------- ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                Comments />                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                    FWHM 0.57288181298114 um                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                Facility HST                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
    FilterProfileService ivo://svo/fps                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
              Instrument WFC3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
                  MagSys Vega                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
               PhotCalID HST/WFC3_IR.G141/Vega                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
              PhotSystem WFC3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
        ProfileReference http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
           WavelengthCen 1.3890218043528004 um                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
                     ... ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
           WavelengthUCD em.wl                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
                WidthEff 0.52172163640941 um                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
               ZeroPoint 1330.8028480785                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
           ZeroPointType Pogson                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
           ZeroPointUnit Jy                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                 centers [[  1.01926517e+00   1.05965090e+00   1.10163665e+00   1.14528608e+00
        1.19066489e+00   1.23784184e+00   1.28688800e+00   1.33787763e+00
        1.39088738e+00   1.44599760e+00   1.50329149e+00   1.56285512e+00
        1.62477911e+00   1.68915677e+00   1.75608516e+00]
     [  1.72295573e-03   3.03102843e-02   2.26018965e-01   3.25663447e-01
        3.65870863e-01   4.06875938e-01   4.32371706e-01   4.55397606e-01
        4.56173390e-01   4.64371055e-01   4.66175199e-01   4.43126798e-01
        4.16419178e-01   1.63810879e-01   1.06626963e-02]]
                filterID HST/WFC3_IR.G141                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                  n_bins 15                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
              n_channels 634                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                    path /Users/jfilippazzo/Documents/Modules/svo_filters/svo_filters/filters/WFC3_IR.G141                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                    refs []                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
    Length = 28 rows



![png](svo_demo_files/svo_demo_12_1.png)


This can get even fancier by passing an arbitrary array of length `n_channels` to the `bin_throughput` argument in order to change the shape of the bin throughput.


```python
import numpy as np
G141 = svo.Filter('WFC3_IR.G141', n_bins=15)

# Let's just do a top-hat
throughput = np.ones(G141.n_channels)
throughput[:100] = 0
throughput[-100:] = 0
G141.bin(bin_throughput=throughput)
G141.plot()
```

    15 bins of 634 channels each.
    15 bins of 634 channels each.



![png](svo_demo_files/svo_demo_14_1.png)


Filters are fun!

## Licensed

This project is Copyright (c) Joe Filippazzo and licensed under the terms of the BSD 3-Clause license. See the licenses folder for more information.

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge
