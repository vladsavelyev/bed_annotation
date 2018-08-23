#!/usr/bin/env bash

VERSION=$1
if [ -z $VERSION ] ; then
    echo 'Usage: $1 VERSION'  >&2
    exit 1
fi

cat $VERSION > VERSION.txt
git add VERSION.txt
git commit -m 'Bump version $VERSION'
git tag $VERSION
git push --tags