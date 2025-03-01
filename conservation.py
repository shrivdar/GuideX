import os
import subprocess
import numpy as np
from pathlib import Path
from typing import List
from Bio import SeqIO  # Added import
from Bio.Seq import Seq
from skbio import DNA, TabularMSA
from scipy.spatial.distance import jensenshannon
import plotly.express as px
from .utils.logger import setup_logger

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Conservation analysis with MAFFT integration"""
    
    def __init__(self, window_size: int = 30):
        self.window_size = window_size

    def align_genomes(self, genomes, output_dir="alignments"):
        """Fixed MAFFT implementation"""
        # Create directory if needed
        os.makedirs(output_dir, exist_ok=True)

        # Write combined input file
        input_path = os.path.join(output_dir, "combined_input.fasta")
        with open(input_path, "w") as f:
            for genome in genomes:
                if not isinstance(genome.seq, Seq):  # Prevent deprecation
                    genome.seq = Seq(str(genome.seq))
                SeqIO.write(genome, f, "fasta")

        # Run MAFFT safely
        output_path = os.path.join(output_dir, "aligned.fasta")
        cmd = [
            "mafft",
            "--auto",
            "--thread", "1",
            input_path
        ]
        
        with open(output_path, "w") as outfile:
            result = subprocess.run(
                cmd, 
                stdout=outfile, 
                stderr=subprocess.PIPE,
                text=True
            )
            
        if result.returncode != 0:
            raise RuntimeError(f"MAFFT failed: {result.stderr}")
            
        return output_path
    
    def calculate_jsd(self, aligned_file: Path) -> List[float]:
        """Calculate conservation scores using Jensen-Shannon divergence."""
        msa = TabularMSA.read(aligned_file, constructor=DNA)
        return self._sliding_window_analysis(msa)

    def _sliding_window_analysis(self, msa):
        """Perform windowed conservation analysis."""
        scores = []
        for i in range(0, len(msa[0]), self.window_size):
            window = msa[:, i:i+self.window_size]
            scores.extend([self._window_jsd(window)] * self.window_size)
        return scores

    def _window_jsd(self, window):
        """Calculate JSD for a single window."""
        freq_matrix = []
        for seq in window:
            freq = [seq.frequencies().get(nt, 0) for nt in 'ACGT']
            freq_matrix.append(freq)
        
        if len(freq_matrix) < 2:
            return 0.0  # Handle single-sequence edge case
            
        return np.mean([
            jensenshannon(freq_matrix[0], freq) ** 2
            for freq in freq_matrix[1:]
        ])

    def plot_conservation(self, scores: List[float], output_file: Path) -> None:
        """Generate interactive conservation plot."""
        fig = px.line(
            x=list(range(len(scores))),
            y=scores,
            labels={"x": "Genomic Position", "y": "Conservation Score"},
            title="GuideX Conservation Profile"
        )
        fig.write_html(str(output_file))
