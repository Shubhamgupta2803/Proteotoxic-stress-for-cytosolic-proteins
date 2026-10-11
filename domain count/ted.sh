 echo "=================================================="
 echo "Running analysis pipeline for k_fraction = 0.9"
/Library/Developer/CommandLineTools/usr/bin/python3 ted.py \
  --indir 90/final_merged_output_localization.csv \
  --comparisons ">=15" "<15:interactome" ">=20" \
  --outdir 90/ted/ 
 echo "=========pipeline completed for 0.9==============="


 echo "=================================================="
 echo "Running analysis pipeline for k_fraction = 0.8"
/Library/Developer/CommandLineTools/usr/bin/python3 ted.py \
  --indir 80/final_merged_output_localization.csv \
  --comparisons ">=30" "<30:interactome" ">=35" \
  --outdir 80/ted/ 
 echo "=========pipeline completed for 0.8==============="

 echo "=================================================="
 echo "Running analysis pipeline for k_fraction = 0.7"
/Library/Developer/CommandLineTools/usr/bin/python3 ted.py \
  --indir 70/final_merged_output_localization.csv \
  --comparisons ">=40" "<40:interactome" ">=45" \
  --outdir 70/ted/ 
 echo "=========pipeline completed for 0.7==============="

 echo "=================================================="
 echo "Running analysis pipeline for k_fraction = 0.6"
/Library/Developer/CommandLineTools/usr/bin/python3 ted.py \
  --indir 60/final_merged_output_localization.csv \
  --comparisons ">=50" "<50:interactome" ">=55" \
  --outdir 60/ted/ 
 echo "=========pipeline completed for 0.6==============="

 echo "=================================================="
 echo "Running analysis pipeline for k_fraction = 0.5"
/Library/Developer/CommandLineTools/usr/bin/python3 ted.py \
  --indir 50/final_merged_output_localization.csv \
  --comparisons ">=60" "<60:interactome" ">=65"  \
  --outdir 50/ted/ 
 echo "=========pipeline completed for 0.5==============="
