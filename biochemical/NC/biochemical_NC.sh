#!/bin/bash

for DIR in 90 80 70 60 50
do
    case $DIR in
        90) COMP="I_non_>=15,T_non_<15;I_non_>=20,T_non_<20" ;;
        80) COMP="I_non_>=30,T_non_<30;I_non_>=35,T_non_<35" ;;
        70) COMP="I_non_>=40,T_non_<40;I_non_>=45,T_non_<45" ;;
        60) COMP="I_non_>=50,T_non_<50;I_non_>=55,T_non_<55" ;;
        50) COMP="I_non_>=60,T_non_<60;I_non_>=65,T_non_<65" ;;
    esac

    KF=$(echo "$DIR / 100" | bc -l)

    echo "=========================================="
    echo "Running biochemical plots for L_${DIR}"

    /Library/Developer/CommandLineTools/usr/bin/python3 biochemical_NC_plot.py \
        --features "$DIR/NC/biochemical" \
        --output "$DIR/plots/biochemical_NC" \
        --k_fraction "$KF" \
        --compare "$COMP"

    echo "Finished L_${DIR}"
done