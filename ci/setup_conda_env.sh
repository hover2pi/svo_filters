#!/bin/bash

echo "Creating a Python $PYTHON_VERSION environment"
conda create -n svo python=$PYTHON_VERSION || exit 1
source activate svo

echo "Installing packages..."
conda install numpy
conda install -c astropy astropy-helpers
conda install -c conda-forge extension-helpers
conda install flake8 beautifulsoup4 lxml astropy
git clone https://github.com/astropy/astroquery.git
cd astroquery
python setup.py install
pip install pytest pytest-cov coveralls