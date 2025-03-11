# Phage host interaction
This folder includes a preprocessing script to combine a phage host interaction matrix with bacterial and phage pangenome.

## Input
Please see the files in `test_small/` for example formats. 
The pangenomes need to be in roary format. 
The interaction matrix is given with phages as lines and bacteria as columns, the inverse encoding is also allowed with the `-j` option.

## Output
The output is written in roary format and includes all lines included in the given bacterial pangenome.

In the **phage presence** mode, one line is added for each phage, where infecting phages are encoded as present and non-infecting phages as absent.
The reverse encoding will encode infecting phages as absent and non-infecting phages as present.
To distinguish phages from bacterial genes, `_Phag` is added to the name of each phage.

In the **phage gene** mode, one line is added for each phage gene, where all genes that are present in phages that infect a bacterium are marked as present in that bacterium.
To distinguish genes names between phages and bacteria, `_Phag` is added to the name of each phage gene.
The reverse encoding is also available in this mode for completeness, but is less relevant.

## Usage
Example usage for adding **phage presence** to a bacterial pangenome:

```commandline
python pan_interactions.py -i test_small/matrix.txt -m test_small/bact_map.txt -b test_small/roary_mini_example.csv -o test_small/roary_output.csv -r test_small/roary_output_rev.csv 

python pan_interactions.py -i test_small/matrix2.txt -j -m test_small/bact_map.txt -b test_small/roary_mini_example.csv -o test_small/roary_output2.csv 
```

Example usage for adding **phage genes** to a bacterial pangenome:

```commandline
python pan_interactions.py -i test_small/matrix.txt -m test_small/bact_map.txt -n test_small/phage_map.txt -b test_small/roary_mini_example.csv -p test_small/roary_phage.csv -o test_small/roary_output_genes.csv 
```

By default, any gene that is present in any phage that can infect a particular bacterium is encoded as present in that bacterium. With `-a` it can be enforced that only genes present in ALL phages that infect a bacterium are encoded as present

```commandline
python pan_interactions.py -a -i test_small/matrix.txt -m test_small/bact_map.txt -n test_small/phage_map.txt -b test_small/roary_mini_example.csv -p test_small/roary_phage.csv -o test_small/roary_output_genesall.csv 
```

Script usage:
```commandline
usage: pan_interactions.py [-h] -i INTER [-s SEP] [-j] -b BACT [-p PHAG] [-m BACT_MAP] [-n PHAGE_MAP] [-a] -o OUTF [-r OUTF_REV]

Combine pangenomes and interaction matrix

optional arguments:
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
  -a, --phage_gene_AND  Only mark phage genes as present if they are present in ALL phages that can infect a bacterium (default:off)
  -o OUTF, --out OUTF   output file
  -r OUTF_REV, --out_rev OUTF_REV
                        output file for reverse encoding of phages
```
