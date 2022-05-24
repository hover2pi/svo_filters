#!/usr/bin/env python
import os
from setuptools import setup, find_packages

REQUIRES = ['astropy',
            'astroquery',
            'bokeh',
            'ipython',
            'matplotlib',
            'numpy',
            'numpydoc',
            'pysynphot',
            'pytest',
            'pyyaml']

DEPENDENCY_LINKS = [
    'git+https://github.com/astropy/astroquery.git@ccc96185beeff86f3a12a31a00a801afcebe1dbe']

FILES = []
for root, _, files in os.walk("svo_filters"):
    FILES += [os.path.join(root.replace("svo_filters/", ""), fname)
              for fname in files if not fname.endswith(".py") and not fname.endswith(".pyc")]

setup(
    name='svo_filters',
    version='0.4.2',
    description='A Python wrapper for the SVO Filter Profile Service',
    packages=find_packages(
        ".",
        exclude=["*.tests"]),
    package_data={
        'svo_filters': FILES},
    install_requires=REQUIRES,
    dependency_links=DEPENDENCY_LINKS,
    author='Joe Filippazzo',
    author_email='jfilippazzo@stsci.edu',
    license='MIT',
    url='https://github.com/hover2pi/svo_filters',
    long_description='',
    zip_safe=True,
    use_2to3=False)
