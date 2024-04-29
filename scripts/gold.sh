mkdir gold_time 2>/dev/null
touch gold_time/start.time
python goldfinder/goldfinder/goldfinder.py -i data/gene_presence_absence.csv -f roary -t data/core-gps_fasttree.newick -pcor bonferroni -o gold_time/ -O
touch gold_time/end.time
