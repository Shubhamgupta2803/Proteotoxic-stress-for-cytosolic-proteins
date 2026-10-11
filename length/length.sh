
#!/bin/bash

PYTHON="/Library/Developer/CommandLineTools/usr/bin/python3"

for KF in 0.9 0.8 0.7 0.6 0.5
do
  DIR=$(echo "$KF * 100" | bc | cut -d'.' -f1)

  # -----------------------------
  # Set comparisons per K
  # -----------------------------
  case $KF in
    0.9)
      COMP="I_mito_>=15,I_mito_<15;I_non_>=15,I_non_<15"
      ;;
    0.8)
      COMP="I_mito_>=30,I_mito_<30;I_non_>=30,I_non_<30"
      ;;
    0.7)
      COMP="I_mito_>=40,I_mito_<40;I_non_>=40,I_non_<40"
      ;;
    0.6)
      COMP="I_mito_>=50,I_mito_<50;I_non_>=50,I_non_<50"
      ;;
    0.5)
      COMP="I_mito_>=60,I_mito_<60;I_non_>=60,I_non_<60"
      ;;
  esac

  echo "========== Running analysis pipeline for $KF =========="
  echo "Comparisons: $COMP"

  $PYTHON length.py \
    --annotated ${DIR}/ted \
    --output ${DIR}/plots/length \
    --k_fraction $KF \
    --compare "$COMP"

  echo "----- plotting finished for $KF -----"
done