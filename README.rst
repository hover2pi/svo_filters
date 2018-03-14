Introduction
------------

The [Spanish Virtual Observatory (SVO) Filter Profile Service](http://svo2.cab.inta-csic.es/theory/fps3/) is a great resource for homogenized photometric filter curves and metadata. With `svo_filters`, I tried to create a lightweight and flexible package to incorporate these filters into Python applications.

Install in the usual fashion with

.. code-block:: python
    pip install svo_filters

Then import like so


.. code-block:: python
    from svo_filters import svo

The actual filters are stored locally as XML files and can be viewed with


.. code-block:: python
    svo.filters()


Filters are fun!
