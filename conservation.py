import os
import subprocess
import numpy as np
from pathlib import Path
from typing import List
from Bio import SeqIO
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
        """MAFFT alignment with rigorous path handling"""
        # Convert output_dir to absolute path
        output_dir = Path(output_dir).absolute()
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Validate input genomes
        if not isinstance(genomes, list) or len(genomes) < 2:
            raise ValueError("Need ≥2 SeqRecords for alignment")
            
        # Sanitize sequences
        valid_genomes = []
        for idx, rec in enumerate(genomes):
            if not isinstance(rec.seq, Seq):
                rec.seq = Seq(str(rec.seq))
            seq_str = str(rec.seq).upper().replace("-", "")
            if len(seq_str) < 100:
                continue
            valid_genomes.append(SeqRecord(Seq(seq_str), id=f"Genome_{idx+1}"))

        # Write combined input (atomic write)
        input_path = output_dir/"MAFFT_IN.fasta"
        with open(input_path, "w") as f:
            count = SeqIO.write(valid_genomes, f, "fasta")
            logger.info(f"Wrote {count} sequences to {input_path}")

        # Validate input file
        if not input_path.exists():
            raise FileNotFoundError(f"MAFFT input file not created: {input_path}")
        if count < 2:
            raise ValueError(f"Only {count} valid sequences for alignment")

        # Build MAFFT command
        output_path = output_dir/"MAFFT_OUT.fasta"
        cmd = [
            "mafft",
            "--auto",
            "--thread", "1",
            "--quiet",
            str(input_path)
        ]

        # Execute with error capture
        try:
            with open(output_path, "w") as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
            logger.info(f"MAFFT alignment successful: {output_path}")
        except subprocess.CalledProcessError as e:
            error_msg = (
                f"MAFFT failed (code {e.returncode})\n"
                f"Command: {' '.join(cmd)}\n"
                f"Error: {e.stderr.strip()}"
            )
            logger.error(error_msg)
            output_path.unlink(missing_ok=True)  # Clean failed output
            raise RuntimeError(error_msg)

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
