#!/usr/bin/env bash
set -x
set -e

# Builds and uploads packages for linux-32, linux-64, and macos, for py2 and py3

. $HOME/miniconda3/etc/profile.d/conda.sh

if ! command -v conda; then
    echo "conda executble not in PATH: $PATH"
    exit
fi
if ! command -v anaconda; then
    echo "anaconda executble not in PATH: $PATH"
    exit
fi

function build() {
    local NAME=$1
    echo "Building $NAME";

    CHANNELS="-c vladsaveliev -c bioconda -c defaults -c conda-forge"
    for PY in 3.6 ; do
        PACKAGE_PATH=$(conda build $NAME $CHANNELS --output --py $PY | tail -n1)
        if [ -f $PACKAGE_PATH ] ; then
            echo "$PACKAGE_PATH exists, skipping"
        else
            echo "Building $PACKAGE_PATH"
            conda build $NAME $CHANNELS --py $PY
            anaconda upload $PACKAGE_PATH
            BASEDIR=$(dirname $(dirname $PACKAGE_PATH))
            FILENAME=$(basename $PACKAGE_PATH)
            echo "Converting packages into $BASEDIR"
            for PLATFORM in osx-64 linux-64 ; do
                if [ -f $BASEDIR/$PLATFORM/$FILENAME ] ; then
                    echo "$BASEDIR/$PLATFORM/$FILENAME exists, skipping"
                else
                    conda convert -p $PLATFORM $PACKAGE_PATH -o $BASEDIR -f
                    anaconda upload $BASEDIR/$PLATFORM/$FILENAME
                fi
            done
        fi
    done
}

if [ -z "$1" ]; then
    echo "Specify folder to build"
    exit 1
fi
build $1

set +x
set +e
