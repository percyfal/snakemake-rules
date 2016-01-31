#!/bin/bash

set -e
set -x

git fetch origin master

MINICONDA="Miniconda3-$MINICONDA_VERSION-Linux-x86_64"
MINICONDA_URL="http://repo.continuum.io/miniconda/$MINICONDA.sh"

wget $MINICONDA_URL
bash $MINICONDA.sh -b -p $HOME/miniconda
rm -rf $MINICONDA.sh

python -V

conda config --set always_yes yes
conda update -q conda
conda info -a
conda config --add channels bokeh
conda config --add channels percyfal

conda build 
