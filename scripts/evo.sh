mkdir evo_time 2>/dev/null
base_dir=$(pwd)
cd evo_time
touch start.time
# EvoScope does not work with ".."
$base_dir/EvoScope/evoscope -t $base_dir/data/EvoScope_input.tsv -T $base_dir/data/core-gps_fasttree.newick -f -o evo
touch end.time