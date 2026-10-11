#!/bin/bash
set -e  # Stop immediately if any command fails

PYTHON="/Library/Developer/CommandLineTools/usr/bin/python3"
BASE="/Users/shubham/labwork/deseq_old"


    # -----------------------------
    # Step 1: Main count pipeline
    # -----------------------------
    $PYTHON script1.py \
        --with_chx configs/config_with_chx.json \
        --without_chx configs/config_without_chx.json \
        --output_dir "/users/shubham/labwork/deseq_old/90" \
        --k_fraction 0.9

    echo "----- script1.py finished -----"


    # -----------------------------
    # Step 1: Main count pipeline
    # -----------------------------
    $PYTHON script1.py \
        --with_chx configs/config_with_chx.json \
        --without_chx configs/config_without_chx.json \
        --output_dir "/users/shubham/labwork/deseq_old/80" \
        --k_fraction 0.8

    echo "----- script1.py finished -----"

    # -----------------------------
    # Step 1: Main count pipeline
    # -----------------------------
    $PYTHON script1.py \
        --with_chx configs/config_with_chx.json \
        --without_chx configs/config_without_chx.json \
        --output_dir "/users/shubham/labwork/deseq_old/70" \
        --k_fraction 0.7

    echo "----- script1.py finished -----"

    # -----------------------------
    # Step 1: Main count pipeline
    # -----------------------------
    $PYTHON script1.py \
        --with_chx configs/config_with_chx.json \
        --without_chx configs/config_without_chx.json \
        --output_dir "/users/shubham/labwork/deseq_old/60" \
        --k_fraction 0.6

    echo "----- script1.py finished -----"

    # -----------------------------
    # Step 1: Main count pipeline
    # -----------------------------
    $PYTHON script1.py \
        --with_chx configs/config_with_chx.json \
        --without_chx configs/config_without_chx.json \
        --output_dir "/users/shubham/labwork/deseq_old/50" \
        --k_fraction 0.5

    echo "----- script1.py finished -----"