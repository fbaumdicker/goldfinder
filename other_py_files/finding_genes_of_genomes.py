import glob
from Bio import SeqIO





files = glob.glob("/home/chris/Documents/amy_test/**/*.faa", recursive=True)

for file in files:
    fasta_sequences = SeqIO.parse(open(file), 'fasta')
    
    count = 0
    for fasta in fasta_sequences:
        count += 1
    print(count)
