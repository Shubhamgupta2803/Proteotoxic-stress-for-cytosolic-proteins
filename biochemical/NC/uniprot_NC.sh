echo "========== Running analysis pipeline for 0.9 =========="


/Library/Developer/CommandLineTools/usr/bin/python3 uniprot_NC.py \
    --indir 90/NC \
    --output 90/NC/biochemical/ \
    --iupred_path /Users/shubham/Downloads/iupred2a/iupred2a.py

echo "----- uniprot finished -----"
echo "-----Biochemical properties extraction of proteins with CHX for 90% count is completed-----"

echo "========== Running analysis pipeline for 0.8 =========="


/Library/Developer/CommandLineTools/usr/bin/python3 uniprot_NC.py \
    --indir 80/NC \
    --output 80/biochemical/ \
    --iupred_path /Users/shubham/Downloads/iupred2a/iupred2a.py

echo "----- uniprot finished -----"
echo "-----Biochemical properties extraction of proteins with CHX for 80% count is completed-----"

echo "========== Running analysis pipeline for 0.7 =========="


/Library/Developer/CommandLineTools/usr/bin/python3 uniprot_NC.py \
    --indir 70/NC \
    --output 70/NC/biochemical/ \
    --iupred_path /Users/shubham/Downloads/iupred2a/iupred2a.py

echo "----- uniprot finished -----"
echo "-----Biochemical properties extraction of proteins with CHX for 70% count is completed-----"

echo "========== Running analysis pipeline for 0.6 =========="


/Library/Developer/CommandLineTools/usr/bin/python3 uniprot_NC.py \
    --indir 60/NC \
    --output 60/NC/biochemical/ \
    --iupred_path /Users/shubham/Downloads/iupred2a/iupred2a.py

echo "----- uniprot finished -----"
echo "-----Biochemical properties extraction of proteins with CHX for 60% count is completed-----"

echo "========== Running analysis pipeline for 0.5 =========="


/Library/Developer/CommandLineTools/usr/bin/python3 uniprot_NC.py \
    --indir 50/NC \
    --output 50/NC/biochemical/ \
    --iupred_path /Users/shubham/Downloads/iupred2a/iupred2a.py

echo "----- uniprot finished -----"
echo "-----Biochemical properties extraction of proteins with CHX for 50% count is completed-----"