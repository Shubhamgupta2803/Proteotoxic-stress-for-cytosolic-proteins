for DIR in 90 80 70 60 50
do

  case $DIR in
    90) COMP="I_non_>=15,T_non_<15;I_non_>=20,T_non_<20" ;;
    80) COMP="I_non_>=30,T_non_<30;I_non_>=35,T_non_<35" ;;
    70) COMP="I_non_>=40,T_non_<40;I_non_>=45,T_non_<45" ;;
    60) COMP="I_non_>=50,T_non_<50;I_non_>=55,T_non_<55" ;;
    50) COMP="I_non_>=60,T_non_<60;I_non_>=65,T_non_<65" ;;
  esac

  echo "Running for DIR=$DIR"

  /Library/Developer/CommandLineTools/usr/bin/python3 biochemical_plot.py \
    --features $DIR/biochemical \
    --output $DIR/plots/biochemical/ \
    --k_fraction 0.$DIR \
    --compare "$COMP"

done