import json

with open('uniprot_json.json') as uniprot_json:
    uniprot_data = json.load(uniprot_json)

print(uniprot_data.keys())