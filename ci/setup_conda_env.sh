#!/bin/bash

echo "Creating a Python $PYTHON_VERSION environment"
conda create -n svo python=$PYTHON_VERSION || exit 1
source activate svo

echo "Installing packages..."
conda install flake8 beautifulsoup4 lxml numpy astropy