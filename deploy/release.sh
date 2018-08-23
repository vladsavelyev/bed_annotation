#!/usr/bin/env bash

VERSION=$1
if [ -z $VERSION ] ; then
    echo "Usage: $1 VERSION" >&2
    exit 1
fi

echo $VERSION > VERSION.txt
git add VERSION.txt
git commit -m "Release $VERSION"
git tag $VERSION
git push --tags
