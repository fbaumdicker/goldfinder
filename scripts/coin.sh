mkdir coin_time 2>/dev/null
touch coin_time/start.time
coinfinder -i data/gene_presence_absence.csv -I -p data/core-gps_fasttree.newick --all -o coin_time/coin
touch coin_time/end.time
