#!/bin/bash

echo "Creating a Python $PYTHON_VERSION environment"
conda create -n svo python=$PYTHON_VERSION || exit 1
source activate svo

echo "Installing packages..."
conda install numpy
conda install -c astropy astropy-helpers
conda install -c conda-forge extension-helpers
conda install -c conda-forge pytest-astropy-header
conda install flake8 beautifulsoup4 lxml astropy
git clone https://github.com/astropy/astroquery.git
cd astroquery
pip install -e .
pip install pytest pytest-cov coveralls