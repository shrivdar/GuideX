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
        """Fixed MAFFT implementation with strict input validation"""
        # Convert Path objects to strings
        output_dir = str(output_dir)
        os.makedirs(output_dir, exist_ok=True)

        # Ensure we have SeqRecords with Seq objects
        valid_genomes = []
        for g in genomes:
            if not isinstance(g.seq, Seq):
                g.seq = Seq(str(g.seq))
            valid_genomes.append(g)
    
        # Write combined input (REPLACE existing files)
        input_path = os.path.join(output_dir, "MAFFT_IN.fasta")
        with open(input_path, "w") as f:
            SeqIO.write(valid_genomes, f, "fasta")

        # Verify input file contains sequences
        if os.path.getsize(input_path) == 0:
            raise ValueError("Empty MAFFT input file - check genome data")

        # Run MAFFT with strict error checking
        output_path = os.path.join(output_dir, "MAFFT_OUT.fasta")
        cmd = [
            "mafft",
            "--auto",
            "--thread", "1",
            "--quiet",  # Suppress help text
            input_path
        ]
        
        try:
            with open(output_path, "w") as outfile:
                result = subprocess.run(
                    cmd,
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True  # Raise error on non-zero exit
                )
        except subprocess.CalledProcessError as e:
            error_msg = f"""
            MAFFT failed with code {e.returncode}
            Command: {' '.join(cmd)}
            Error output: {e.stderr}
            """
            raise RuntimeError(error_msg) from None
    
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
