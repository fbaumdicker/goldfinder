from Bio import SeqIO
import csv

def search_and_extract(csv_file, search_word):
    with open("/home/franz/syncedprojects/RoyalSocietyManuskriptData/phageroary.csv", newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        for row in reader:
            if any(search_word in cell for cell in row):  # Prüft, ob Suchwort in einer beliebigen Spalte vorkommt
                first_entry = row[0].split(',')[0]  # Nimmt den ersten Eintrag der Zeile vor dem ersten Komma
                return(first_entry)

file_path = "allphageantidefensegenes.csv"

defensegroups = {}

def add_to_defensegroups(key, value):
    if key not in defensegroups:
        defensegroups[key] = set()  # Neue Menge erstellen, falls Schlüssel nicht existiert
    defensegroups[key].add(value)

with open(file_path, 'r') as file:
    for line in file:
        parts = line.split(",")
        protein_id = "cds-" + parts[1]
        defsystem = parts[3]
        group = search_and_extract(file_path, protein_id)
        print([group, defsystem])
        add_to_defensegroups(group, defsystem)
        
# Save to CSV
with open("phage_antidefensegroups.csv", "w", newline="") as f:
    writer = csv.writer(f)
    for key, values in defensegroups.items():
        writer.writerow([key, "+".join(values)])  # Convert set to string
