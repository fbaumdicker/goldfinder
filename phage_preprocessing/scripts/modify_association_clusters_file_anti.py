import re

def load_defense_systems(defense_file):
    """Lädt die Zuordnungen von Gruppen zu Verteidigungssystemen aus der zweiten Datei."""
    defense_dict = {}
    with open(defense_file, 'r') as df:
        for line in df:
            parts = line.strip().split(',')
            if len(parts) == 2:
                defense_dict[parts[0]] = parts[1]
    return defense_dict

def modify_first_file(input_file, output_file, defense_dict):
    """Ersetzt die Gruppenbezeichnungen in der ersten Datei basierend auf der Zuordnungstabelle."""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            match = re.match(r'^(group_\d+)$', line.strip())
            if match:
                group_name = match.group(1)
                if group_name in defense_dict:
                    line = f"{group_name}_ANTIDEF_{defense_dict[group_name]}\n"
            outfile.write(line)

if __name__ == "__main__":
    group_defense_file = "phage_antidefensegroups2.csv"
    cluster_file = "test/mgenes/association_clusters.txt"  
    newcluster_file = "test/mgenes/association_clusters_withANTIDEFtags.txt"  # Ausgabe-Datei
    
    defense_dict = load_defense_systems(group_defense_file)
    modify_first_file(cluster_file, newcluster_file, defense_dict)
    
    print(f"Die Datei wurde erfolgreich verarbeitet und als {newcluster_file} gespeichert.")
