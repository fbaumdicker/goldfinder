from Bio import SeqIO
import csv
import argparse

# Paths to input files
#gbff_file = "/home/franz/tempDEFENSETAG/ncbi_dataset/data/GCA_002874555.1/genomic.gbff"
#defensefinder_results = "/home/franz/tempDEFENSETAG/transportfolder/GCA_002874555.1_ASM287455v1_genomic/GCA_002874555.1_ASM287455v1_genomic_defense_finder_genes.tsv"
#prt_filepath = "/home/franz/tempDEFENSETAG/transportfolder/GCA_002874555.1_ASM287455v1_genomic/GCA_002874555.1_ASM287455v1_genomic.prt"



def parse_prt(file_path, identifiers):
    extracted_data = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                parts = line.split(" # ")
                header_id = parts[0][1:]  # Remove '>'
                
                if header_id in identifiers:
                    numbers = parts[1:4]  # Extract the first three numbers after '#'
                    extracted_data[header_id] = numbers
    
    return extracted_data


def save_to_csv(data, filename):
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)




def main():
    parser = argparse.ArgumentParser(description="Identify CDS with DefenseFinder results")
    parser.add_argument("-g", "--gbff_file", required=True, help="Path to the input gbff file.")
    parser.add_argument("-p", "--prt_filepath", required=True, help="Path to the input prt file as given by DefenseFinder.")
    parser.add_argument("-d", "--defensefinder_results", required=True, help="Path to the file with DefenseFinder results")
    parser.add_argument("-o", "--output", required=True, help="Path to save the list of defense genes")

    args = parser.parse_args()

    # Load defense gene names from DefenseFinder results
    defense_hits = []
    with open(args.defensefinder_results, 'r') as df_file:
        next(df_file)  # Skip the header line
        for line in df_file:
            columns = line.strip().split("\t")
            if columns[-1] == "Defense":  # Ensure the type is 'Defense'
                replicon = columns[0].split('.')[0]  # Contig name (replicon column) removing the suffix
                hit_id = columns[1]
                #hit_pos = int(columns[3])  # Position in prt file
                gene_name = columns[2]  # Gene name for reference
                def_type = columns[-3]
                defense_hits.append((replicon, hit_id, gene_name, def_type))

    identifiers = {entry[1] for entry in defense_hits} # Extract second element from each tuple
    prt_dict = parse_prt(args.prt_filepath, identifiers)

    #Modify the .gbff file
    #output_gbff_file = "/home/guerrero/data/gbff_modified_files/modified_genome_DEF.gbff"
    print(prt_dict)


    defense_genes = []
    #all_genes = []

    with open(args.gbff_file, "r") as gbff:
        for record in SeqIO.parse(gbff, "genbank"):
            print("---------------------------------------------------------------------")
            contig_name = record.name  # The contig name (from the LOCUS line in the GBFF)
            for feature in record.features:
                if feature.type == "CDS":
                    start = int(feature.location.start)
                    end = int(feature.location.end)

                    # Check for overlaps with defense genes
                    for replicon, hit_id, gene_name, def_type in defense_hits:
                        if replicon == contig_name:
                            positions = prt_dict[hit_id]
                            hit_start, hit_end =  int(positions[0]), int(positions[1])
                            if (start  <= hit_end and end >= hit_start) or (start <= hit_start and hit_start <= end) or (start <= hit_end and hit_end <= end):
                                #print("found one!")
                                protein_id = "NA"
                                locus_tag = "NA"

                                if "protein_id" in feature.qualifiers:
                                    protein_id = feature.qualifiers["protein_id"][0]
                                if "locus_tag" in feature.qualifiers:
                                    locus_tag = feature.qualifiers["locus_tag"][0]

                                defense_genes.append((locus_tag, protein_id, "defense", def_type, gene_name, hit_id))
                                #all_genes.append((locus_tag, protein_id, "defense", def_type, gene_name, hit_id))



    print(defense_genes)

    save_to_csv(defense_genes, args.output)


if __name__ == "__main__":
    main()

