from Bio import Entrez
Entrez.email = "darsh.shri123@gmail.com"

try:
    handle = Entrez.esearch(
        db="nucleotide", 
        term="Influenza A virus[Organism]", 
        retmax=1
    )
    record = Entrez.read(handle)
    print(f"NCBI connection successful! Found {record['Count']} records")
except Exception as e:
    print(f"NCBI connection failed: {str(e)}")
