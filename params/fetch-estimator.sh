#!/bin/sh
# The lattice estimator used for the MLWE security figures.
set -e
REV=f18533a19433f6fb1d9fb396006f462adc6b8ad3
test -d lattice-estimator || git clone https://github.com/malb/lattice-estimator
cd lattice-estimator
git checkout -q "$REV"
git --no-pager log -1 --format='malb/lattice-estimator at %H (%ci)'
