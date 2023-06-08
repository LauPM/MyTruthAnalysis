#!/bin/bash
# run with $ source config.sh $ to start using larsoft

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
source localProducts_larsoft__e20_prof/setup

# kx509
# voms-proxy-init -noregen -rfc -voms dune

setup larsoft v09_58_02  -q e20:prof  
#setup larbatch v01_57_01

mrbslp #set local products

# mrbsetenv #local
# ups active

#BUILD#
# mrb i -j16 #NUNCA MAS#
# make install (-j16?) #en la carpeta build que has cambiado