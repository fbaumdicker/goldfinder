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
    gfilename=$(ls gbff-files-bacteria/ncbi_dataset/data/ | grep "$acc\.")
    pfilename=$(find DefenseFinder_Output/$acc.*_genomic/$acc.*.prt)
    dfilename=$(find DefenseFinder_Output/$acc.*_genomic/$acc.*genes.tsv)
    python scripts/identify_DEF_genes.py -g "gbff-files-bacteria/ncbi_dataset/data/$gfilename/genomic.gbff" -p "$pfilename" -d "$dfilename" -o "defense_genes_csv/$acc-defense_genes.csv"

done < "$file_list"
