#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

clean_openfoam .

rm logRun errRun
