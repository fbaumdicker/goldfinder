# Phage host interaction
This folder includes a preprocessing script to combine a phage host interaction matrix with bacterial and phage pangenome.

Example usage for adding phage presence to a bacterial pangenome:

```commandline
python pan_interactions.py -i test_small/matrix.txt -m test_small/bact_map.txt -b test_small/roary_mini_example.csv -o test_small/roary_output.csv -r test_small/roary_output_rev.csv 

python pan_interactions.py -i test_small/matrix2.txt -j -m test_small/bact_map.txt -b test_small/roary_mini_example.csv -o test_small/roary_output2.csv 
```

Example usage for adding phage genes to a bacterial pangenome:

```commandline
python pan_interactions.py -i test_small/matrix.txt -m test_small/bact_map.txt -n test_small/phage_map.txt -b test_small/roary_mini_example.csv -p test_small/roary_phage.csv -o test_small/roary_output_genes.csv 
```

Script usage:
```commandline
usage: pan_interactions.py [-h] -i INTER [-s SEP] [-j] -b BACT [-p PHAG] [-m BACT_MAP] [-n PHAGE_MAP] -o OUTF [-r OUTF_REV]

Combine pangenomes and interaction matrix

options:
  -h, --help            show this help message and exit
  -i INTER, --interactions INTER
                        Interaction matrix (empty cell or 0 means negative interactions, all other entries indicate positive interaction)
  -s SEP, --separator SEP
                        Separator in interactions matrix, default ','
  -j, --phage_column    Interaction matrix has phages as columns and bacteria as lines instead of phages as lines and bacteria as columns (default)
  -b BACT, --bact_pan BACT
                        Bacterial pangenome in roary format
  -p PHAG, --phage_pan PHAG
                        Phage pangenome in roary format (optional)
  -m BACT_MAP, --bact_map BACT_MAP
                        Map (tab-separated file) with bacteria names in interaction matrix (first column) and pangenome (second column)
  -n PHAGE_MAP, --phage_map PHAGE_MAP
                        Map (tab-separated file) with phage names in interaction matrix (first column) and pangenome (second column)
  -o OUTF, --out OUTF   output file
  -r OUTF_REV, --out_rev OUTF_REV
                        output file for reverse encoding of phages
```
