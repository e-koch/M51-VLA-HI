#!/usr/bin/env python
#
# Run this script OUTSIDE CASA in an environment that has astropy,
# spectral-cube, scipy, and numpy.
#
# Reproduces all derived products (convolved cubes, noise maps, masks,
# moment maps) for the M51 VLA HI data.
#

##############################################################################
# Load routines, change directory, initialise handlers
##############################################################################

import os
import sys
import importlib

# Pipeline directory: location of the cloned phangs_imaging_scripts repo
pipedir = '/home/erickoch/phangs_imaging_scripts/'

# Master key pointing to all project keys
key_file = '/home/erickoch/M51_HI/keys/master_key.txt'

os.chdir(pipedir)

from phangsPipeline import phangsLogger as pl
importlib.reload(pl)
pl.setup_logger(level='DEBUG', logfile=None)

from phangsPipeline import handlerKeys as kh
from phangsPipeline import handlerDerived as der
importlib.reload(kh)
importlib.reload(der)

this_kh = kh.KeyHandler(master_key=key_file)
this_der = der.DerivedHandler(key_handler=this_kh)

this_kh.make_missing_directories(imaging=True, derived=True,
                                 postprocess=True, release=True)

##############################################################################
# Set targets, configs, and line products
##############################################################################

# Target: m51 only (the sub-parts m51_1..m51_9 are used only during imaging)
this_der.set_targets(only=['m51'])

# All interferometric configurations present in the data
this_der.set_interf_configs(only=['D', 'C+D', 'B+C+D', 'A+B+C+D'])

# Feathered configurations (GBT-combined), if single-dish data are available
# this_der.set_feather_configs(only=['D+tp', 'C+D+tp', 'B+C+D+tp', 'A+B+C+D+tp'])

# HI line products at each channel width; OH lines excluded from derived pipeline
this_der.set_line_products(only=['hi21cm_5kms', 'hi21cm_10kms'])
this_der.set_no_cont_products(True)

##############################################################################
# Choose which steps to run
##############################################################################

do_convolve   = True   # smooth to angular resolutions in derived_key (6",10",30",60")
do_noise      = True   # fit noise model in signal-free channels
do_strictmask = True   # high-S/N (4σ seed / 2σ wing) 3-D mask
do_broadmask  = True   # union of strict masks across all linked configs
do_moments    = True   # moment maps + error maps for both mask types
do_secondary  = True   # secondary moments that depend on first-round moments (mom1wprior)

##############################################################################
# Execute the pipeline step by step
##############################################################################

# Convolve cubes to the angular resolutions listed in derived_key.txt
# (6as, 10as, 30as, 60as) for each config/product combination.
if do_convolve:
    this_der.loop_derive_products(
        do_convolve=True, do_noise=False, do_strictmask=False,
        do_broadmask=False, do_moments=False, do_secondary=False)

# Estimate per-channel noise from signal-free regions; produces *_noise.fits.
if do_noise:
    this_der.loop_derive_products(
        do_convolve=False, do_noise=True, do_strictmask=False,
        do_broadmask=False, do_moments=False, do_secondary=False)

# Build strict masks: 4σ seed, ≥2 connected channels, 2σ wings.
# Parameters are set in derived_key.txt (strictmask_kw).
if do_strictmask:
    this_der.loop_derive_products(
        do_convolve=False, do_noise=False, do_strictmask=True,
        do_broadmask=False, do_moments=False, do_secondary=False)

# Build broad masks by combining strict masks across all configs listed
# in mask_configs in derived_key.txt (['A+B+C+D','B+C+D','C+D','D']).
if do_broadmask:
    this_der.loop_derive_products(
        do_convolve=False, do_noise=False, do_strictmask=False,
        do_broadmask=True, do_moments=False, do_secondary=False)

# Compute moment maps defined in moment_key.txt:
#   strict: mom0, mom1, mom2, ew  (+ error maps)
#   broad:  mom0, tpeak, tpeak12p5, mom1
if do_moments:
    this_der.loop_derive_products(
        do_convolve=False, do_noise=False, do_strictmask=False,
        do_broadmask=False, do_moments=True, do_secondary=False)

# Secondary moments that depend on first-round results (e.g. mom1wprior
# uses the broad mom1 as a velocity prior).
if do_secondary:
    this_der.loop_derive_products(
        do_convolve=False, do_noise=False, do_strictmask=False,
        do_broadmask=False, do_moments=False, do_secondary=True)
