#!/bin/bash
set -e
set -x

echo $RECIPE_DIR
BLD_DIR=`pwd`
SRC_DIR=$RECIPE_DIR/..

pushd $SRC_DIR

$PYTHON setup.py install --record=record.txt

# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
