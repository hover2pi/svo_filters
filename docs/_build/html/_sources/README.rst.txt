svo_filters
===========

The `Spanish Virtual Observatory (SVO) Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps3/>`_ is a great resource for homogenized photometric filter curves and metadata. With `svo_filters`, I tried to create a lightweight and flexible package to incorporate these filters into Python applications.

Installation
------------

To install **svo_filters**, do::

    pip install svo_filters

Alternatively, you can clone from Github with::

    git clone https://github.com/hover2pi/svo_filters.git
    python svo_filters/setup.py install

Load a Photometric Filter
----------------------------

The actual filters are stored locally as XML files and can be viewed with::

    from svo_filters import svo
    svo.filters()

To create a filter object, pass a bandpass name to the :py:class:`svo.Filter()` class::

    H_band = svo.Filter('2MASS.H')

You can see someinformation about the filter with::

    H_band.info()

And you can plot the bandpass like so::

    H_band.plot()

.. image:: ../svo_filters/data/plots/H.png
    :align: center
    :height: 300px
    :alt: Filter bandpass

Load a Grism
------------

Filters can also be binned arbitrarily, for use with grisms. We can pass integers to the ``n_bins`` or ``pixels_per_bin`` arguments to specify the number of wavelength bins or pixels per bin, respectively::

    G141 = svo.Filter('WFC3_IR.G141', n_bins=15)

.. image:: ../svo_filters/data/plots/G141.png
    :align: center
    :height: 300px
    :alt: Filter bandpass

Apply a Filter to a Spectrum
----------------------------

Filters can be applied to a spectrum by passing a sequence of [W, F] or [W, F, E] with astropy units to the :py:meth:`~svo.Filter.apply` method::

    filtered = G141.apply(spec, plot=True)

.. image:: ../svo_filters/data/plots/filtered.png
    :align: center
    :height: 300px
    :alt: Filter bandpass

Contents
========

.. toctree::
   :maxdepth: 2

   svo

