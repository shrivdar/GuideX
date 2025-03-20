from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
import numpy as np

Entrez.email = "your.email@example.com"  # NCBI requirement

CBSV_ISOLATES = [
    "NC_012698.2", "GU563327.1", "FN434436.1", "MK955888.1",
    "OQ335847.1", "KY563367.1", "MG019915.1", "MK103393.1",
    "MW961170.2", "MZ362877.1", "GU563326.1", "KY290025.1",
    "LT577538.1", "FN434437.1"
]

def fetch_conserved_regions(pipeline_output: str) -> dict:
    """Parse pipeline's conserved regions output"""
    regions = {}
    with open(pipeline_output) as f:
        for line in f:
            if line.startswith("CR_"):
                cr_id, start, end, score = line.strip().split('\t')
                regions[cr_id] = {
                    'start': int(start),
                    'end': int(end),
                    'score': float(score)
                }
    return regions

def fetch_genome_sequences(accessions: list) -> dict:
    """Retrieve full genome sequences from NCBI"""
    handle = Entrez.efetch(
        db="nucleotide",
        id=accessions,
        rettype="fasta",
        retmode="text"
    )
    return {rec.id: rec.seq for rec in SeqIO.parse(handle, "fasta")}

def analyze_conservation(regions: dict, genomes: dict):
    """Comparative conservation analysis using BLAST"""
    results = []
    
    for cr_id, cr_data in regions.items():
        # Extract conserved region sequence from reference
        ref_seq = genomes["NC_012698.2"][cr_data['start']:cr_data['end']]
        
        # Run BLAST against CBSV database
        result = NCBIWWW.qblast(
            program="blastn",
            database="nr",
            sequence=ref_seq,
            entrez_query="Cassava brown streak virus[organism]",
            hitlist_size=50,
            expect=1e-10
        )
        
        # Parse BLAST results
        blast_rec = NCBIXML.read(result)
        for alignment in blast_rec.alignments:
            for hsp in alignment.hsps:
                results.append({
                    'cr_id': cr_id,
                    'accession': alignment.accession,
                    'identity': hsp.identities / hsp.align_length * 100,
                    'evalue': hsp.expect,
                    'coverage': (hsp.query_end - hsp.query_start) / len(ref_seq) * 100
                })
    
    return pd.DataFrame(results)

def generate_conservation_matrix(blast_results: pd.DataFrame) -> pd.DataFrame:
    """Create presence/absence matrix across isolates"""
    matrix = pd.pivot_table(
        blast_results,
        index='cr_id',
        columns='accession',
        values='identity',
        aggfunc=np.max,
        fill_value=0
    )
    return matrix.reindex(columns=CBSV_ISOLATES, fill_value=0)
