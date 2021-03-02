#!/bin/bash

echo "Creating a Python $PYTHON_VERSION environment"
conda create -n svo python=$PYTHON_VERSION || exit 1
source activate svo

echo "Installing packages..."
conda install flake8 beautifulsoup4 lxml numpy astropy astropy_helpers
conda install -c conda-forge extension-helpers
git clone https://github.com/astropy/astroquery.git
cd astroquery
python setup.py install
pip install pytest pytest-cov coveralls