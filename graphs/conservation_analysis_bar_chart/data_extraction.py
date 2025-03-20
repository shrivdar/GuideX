from Bio import AlignIO
from collections import defaultdict

def calculate_conservation(alignment_file):
    alignment = AlignIO.read(alignment_file, "fasta")
    conservation = defaultdict(int)
    
    for rec in alignment:
        for pos, nt in enumerate(rec.seq):
            conservation[pos] += (nt != '-')
    
    total_seqs = len(alignment)
    return {pos: (count/total_seqs)*100 for pos, count in conservation.items()}

# Usage
conservation_data = calculate_conservation("alignments/cbsv_alignment.fasta/ALIGNMENT_OUT.fasta")
