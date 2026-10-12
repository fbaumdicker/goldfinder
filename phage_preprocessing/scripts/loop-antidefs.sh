#!/bin/bash

# Überprüfen, ob die Datei mit der Liste als Argument übergeben wurde
if [ $# -ne 1 ]; then
    echo "Usage: $0 <file_with_accs>"
    exit 1
fi

file_list="$1"

# Überprüfen, ob die Datei existiert
if [ ! -f "$file_list" ]; then
    echo "Error: File '$file_list' not found!"
    exit 1
fi

# Datei Zeile für Zeile einlesen und das Programm ausführen
while IFS= read -r acc; do
    echo "Processing: $acc"
    gfilename=$(ls gbff-files-phages/ncbi_dataset/data/ | grep "$acc\.")
    pfilename=$(find anti-defenseFinderResults_phages/$acc.*_genomic/$acc.*.prt)
    dfilename=$(find anti-defenseFinderResults_phages/$acc.*_genomic/$acc.*genes.tsv)
    python scripts/identify_ANTIDEF_genes.py -g "gbff-files-phages/ncbi_dataset/data/$gfilename/genomic.gbff" -p "$pfilename" -d "$dfilename" -o "antidefense_genes_csv/$acc-antidefense_genes.csv"

done < "$file_list"
