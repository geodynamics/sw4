#!/usr/bin/env bash

set -eu -o pipefail

sudo apt-get update

sudo apt-get install proj-bin libhdf5-mpi-dev libfftw3-dev libfftw3-mpi-dev liblapack-dev python3-h5py curl
