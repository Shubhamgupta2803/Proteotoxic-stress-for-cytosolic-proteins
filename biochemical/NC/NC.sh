#!/bin/bash

PYTHON="/Library/Developer/CommandLineTools/usr/bin/python3"
  $PYTHON NC.py \
    --indir 90/ted/ \
    --outdir 90/NC/ 

PYTHON="/Library/Developer/CommandLineTools/usr/bin/python3"
  $PYTHON NC.py \
    --indir 80/ted/ \
    --outdir 80/NC/ 

PYTHON="/Library/Developer/CommandLineTools/usr/bin/python3"
  $PYTHON NC.py \
    --indir 70/ted/ \
    --outdir 70/NC/ 

PYTHON="/Library/Developer/CommandLineTools/usr/bin/python3"
  $PYTHON NC.py \
    --indir 60/ted/ \
    --outdir 60/NC/ 

PYTHON="/Library/Developer/CommandLineTools/usr/bin/python3"
  $PYTHON NC.py \
    --indir 50/ted/ \
    --outdir 50/NC/ 
   