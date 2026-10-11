echo "========== Running analysis pipeline for 0.9 =========="


/Library/Developer/CommandLineTools/usr/bin/python3 uniprot3.py \
    --indir 90/ted \
    --output 90/biochemical/ \
    --iupred_path /Users/shubham/Downloads/iupred2a/iupred2a.py

echo "----- uniprot finished -----"
echo "-----Biochemical properties extraction of proteins with CHX for 90% count is completed-----"

echo "========== Running analysis pipeline for 0.8 =========="


/Library/Developer/CommandLineTools/usr/bin/python3 uniprot3.py \
    --indir 80/ted \
    --output 80/biochemical/ \
    --iupred_path /Users/shubham/Downloads/iupred2a/iupred2a.py

echo "----- uniprot finished -----"
echo "-----Biochemical properties extraction of proteins with CHX for 80% count is completed-----"

echo "========== Running analysis pipeline for 0.7 =========="


/Library/Developer/CommandLineTools/usr/bin/python3 uniprot3.py \
    --indir 70/ted \
    --output 70/biochemical/ \
    --iupred_path /Users/shubham/Downloads/iupred2a/iupred2a.py

echo "----- uniprot finished -----"
echo "-----Biochemical properties extraction of proteins with CHX for 70% count is completed-----"

echo "========== Running analysis pipeline for 0.6 =========="


/Library/Developer/CommandLineTools/usr/bin/python3 uniprot3.py \
    --indir 60/ted \
    --output 60/biochemical/ \
    --iupred_path /Users/shubham/Downloads/iupred2a/iupred2a.py

echo "----- uniprot finished -----"
echo "-----Biochemical properties extraction of proteins with CHX for 60% count is completed-----"

echo "========== Running analysis pipeline for 0.5 =========="


/Library/Developer/CommandLineTools/usr/bin/python3 uniprot3.py \
    --indir 50/ted \
    --output 50/biochemical/ \
    --iupred_path /Users/shubham/Downloads/iupred2a/iupred2a.py

echo "----- uniprot finished -----"
echo "-----Biochemical properties extraction of proteins with CHX for 50% count is completed-----"