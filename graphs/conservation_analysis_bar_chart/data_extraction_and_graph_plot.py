# guidex/viz/conservation.py
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
from typing import List, Tuple

def calculate_conservation(alignment_path: str) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
    """
    Calculate conservation scores from MUSCLE alignment
    Returns:
        - conservation_scores: Array of conservation percentages (0-100)
        - conserved_regions: List of (start, end) tuples for conserved regions
    """
    alignment = AlignIO.read(alignment_path, "fasta")
    num_seqs = len(alignment)
    length = alignment.get_alignment_length()
    
    # Calculate position conservation
    conservation = np.zeros(length)
    for i in range(length):
        column = alignment[:, i]
        conservation[i] = (sum(1 for base in column if base != '-') / num_seqs) * 100
    
    # Find conserved regions (threshold=80% identity)
    conserved_regions = []
    in_region = False
    start = 0
    for i, score in enumerate(conservation):
        if score >= 80 and not in_region:
            start = i
            in_region = True
        elif score < 80 and in_region:
            conserved_regions.append((start, i-1))
            in_region = False
    if in_region:
        conserved_regions.append((start, len(conservation)-1))
    
    return conservation, conserved_regions

def plot_conservation(conservation: np.ndarray, 
                     regions: List[Tuple[int, int]],
                     invalid_regions: List[Tuple[int, int]] = None,
                     window_size: int = 50):
    """
    Generate conservation bar chart with sliding window averaging
    """
    plt.figure(figsize=(15, 6))
    
    # Apply sliding window average
    window = np.ones(window_size)/window_size
    smoothed = np.convolve(conservation, window, mode='same')
    
    # Plot conservation profile
    plt.bar(range(len(smoothed)), smoothed, 
           color='steelblue', width=1, alpha=0.7,
           label='Conservation Score')
    
    # Highlight conserved regions
    for start, end in regions:
        plt.axvspan(start, end, color='forestgreen', alpha=0.3, 
                   label='Conserved Region' if start == regions[0][0] else "")
    
    # Mark invalid regions
    if invalid_regions:
        for start, end in invalid_regions:
            plt.axvspan(start, end, color='red', alpha=0.2, 
                       label='Invalid Region' if start == invalid_regions[0][0] else "")
    
    plt.xlabel("Genomic Position", fontsize=12)
    plt.ylabel("Conservation (%)", fontsize=12)
    plt.title("CBSV Conservation Analysis", fontsize=16)
    plt.legend()
    plt.grid(alpha=0.2)
    plt.tight_layout()
    return plt.gcf()
