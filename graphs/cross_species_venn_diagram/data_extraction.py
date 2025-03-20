# Using sequence alignment
from Bio import pairwise2

cbsv_targets = ["GGGACCCAAAGGGACCCAAAGGGACCCAAA", ...] 
ucbsv_seqs = ["GGGACCCAAAGGGACCCAAAGTAACCAAA", ...]  # Load from FASTA

matches = []
for grna in cbsv_targets:
    for seq in ucbsv_seqs:
        aln = pairwise2.align.globalxx(grna, seq)
        score = aln[0].score / len(grna)
        if score > 0.9:
            matches.append(grna)
