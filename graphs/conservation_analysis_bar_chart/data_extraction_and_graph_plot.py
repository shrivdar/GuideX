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
    with plt.rc_context(POSTER_STYLE):
        fig, ax = plt.subplots(figsize=(18, 6))
        
        # Smoothed conservation
        window = np.hanning(window_size)
        smoothed = np.convolve(conservation, window, mode='same') / window.sum()
        
        # Main conservation plot
        ax.fill_between(range(len(smoothed)), smoothed, 
                       color='#2c7bb6', alpha=0.3, label='Conservation')
        ax.plot(smoothed, color='#2c7bb6', lw=1.5)
        
        # Region highlighting
        for start, end in regions:
            ax.axvspan(start, end, color='#2ca25f', alpha=0.2, lw=0)
            
        if invalid_regions:
            for start, end in invalid_regions:
                ax.axvspan(start, end, color='#d7191c', alpha=0.3, lw=0,
                          hatch='//', label='Invalid' if start == invalid_regions[0][0] else "")
        
        # Formatting
        ax.set_xlabel("Genomic Position (bp)", labelpad=10)
        ax.set_ylabel("Conservation (%)", labelpad=10)
        ax.set_title("CBSV Genome Conservation Analysis", pad=20)
        
        # Unified legend
        handles, labels = ax.get_legend_handles_labels()
        unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
        ax.legend(*zip(*unique), loc='upper right', frameon=True)
        
        # Annotation panel
        ax.text(0.01, 0.95, 'A', transform=ax.transAxes,
               fontsize=24, weight='bold', va='top')
        ax.text(0.05, 0.95, f"Conserved Regions: {len(regions)}\nWindow Size: {window_size}bp",
               transform=ax.transAxes, va='top', ha='left',
               bbox=dict(facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        return fig
