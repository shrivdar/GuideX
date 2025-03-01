import subprocess
import os
import numpy as np
from pathlib import Path
from typing import List
from skbio import DNA, TabularMSA
from scipy.spatial.distance import jensenshannon
import plotly.express as px
from Bio.SeqRecord import SeqRecord
from .utils.logger import setup_logger

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Identify conserved regions using Jensen-Shannon divergence."""
    
    def __init__(self, algorithm: str = "mafft", window_size: int = 30):
        self.algorithm = algorithm
        self.window_size = window_size

    def align_genomes(self, genomes: List[SeqRecord], output_dir: Path = Path("alignments")) -> Path:
        """Robust MAFFT alignment with comprehensive error handling."""
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # 1. Generate clean input files with sanitized names
        input_files = []
        for idx, genome in enumerate(genomes):
            # Sanitize filename and sequence
            clean_id = genome.id.replace(" ", "_").replace(".", "_").replace("/", "-")[:30]
            clean_seq = str(genome.seq).replace(" ", "").upper()
            
            if not clean_seq:
                raise ValueError(f"Empty sequence for genome {clean_id}")
            
            # Create input file
            input_path = output_dir / f"input_{idx}_{clean_id}.fasta"
            with open(input_path, "w") as f:
                f.write(f">{clean_id}\n{clean_seq}")
            input_files.append(input_path.resolve())

        # 2. Build MAFFT command with absolute paths
        output_path = (output_dir / "aligned.fasta").resolve()
        cmd = [
            "mafft",
            "--auto",
            "--quiet",
            "--thread", "1",  # Single thread for stability
            "--out", str(output_path)
        ] + [str(fp) for fp in input_files]

        # 3. Execute with detailed error handling
        try:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            return output_path
            
        except subprocess.CalledProcessError as e:
            error_msg = f"""
            MAFFT alignment failed!
            Command: {' '.join(cmd)}
            Error Output:
            {e.stderr}
            Input Files: {[str(fp) for fp in input_files]}
            """
            raise RuntimeError(error_msg) from None

    def calculate_jsd(self, aligned_file: Path) -> List[float]:
        """Calculate Jensen-Shannon divergence scores."""
        msa = TabularMSA.read(aligned_file, constructor=DNA)
        scores = []
        for i in range(0, len(msa[0]), self.window_size):
            window = msa[:, i:i+self.window_size]
            window_score = self._jsd_score(window)
            scores.extend([window_score] * self.window_size)
        return scores

    def _jsd_score(self, window):
        """Calculate Jensen-Shannon divergence for a window."""
        freq_matrix = []
        for seq in window:
            freq = [seq.frequencies().get(nt, 0) for nt in 'ACGT']
            freq_matrix.append(freq)
        
        return np.mean([
            jensenshannon(freq_matrix[0], freq) ** 2
            for freq in freq_matrix[1:]
        ])

    def plot_conservation(self, scores: List[float], output_file: Path) -> None:
        """Generate interactive conservation plot."""
        fig = px.line(
            x=list(range(len(scores))), 
            y=scores, 
            labels={"x": "Genome Position", "y": "Conservation Score"},
            title="GuideX Conservation Analysis"
        )
        fig.write_html(output_file)
