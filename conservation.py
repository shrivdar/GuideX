import subprocess
import numpy as np
from skbio import DNA, TabularMSA
from skbio.diversity.alpha import jsd
import plotly.express as px
from pathlib import Path
from typing import List
from .utils.logger import setup_logger

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Identify conserved regions using JSD and entropy."""
    
    def __init__(self, algorithm: str = "mafft", window_size: int = 30):
        self.algorithm = algorithm
        self.window_size = window_size

    def align_genomes(self, genomes: List, output_dir: Path = Path("alignments")) -> Path:
        """Run MAFFT/MUSCLE with auto-parameters."""
        output_dir.mkdir(exist_ok=True)
        fasta_paths = [self._save_temp(genome, output_dir) for genome in genomes]
        
        # Auto-select MAFFT parameters based on genome size
        avg_length = np.mean([len(g.seq) for g in genomes])
        threads = 4 if avg_length > 1e5 else 1  # Use multithreading for large genomes
        
        cmd = (
            f"mafft --thread {threads} --auto {' '.join(map(str, fasta_paths))} "
            f"> {output_dir}/aligned.fasta"
        )
        subprocess.run(cmd, shell=True, check=True)
        return output_dir / "aligned.fasta"

    def calculate_jsd(self, aligned_file: Path) -> List[float]:
        """Calculate Jensen-Shannon divergence scores."""
        msa = TabularMSA.read(aligned_file, constructor=DNA)
        scores = []
        for i in range(0, len(msa[0]), self.window_size):
            window = msa[:, i:i+self.window_size]
            freq_matrix = [
                [seq.frequencies().get(nt, 0) for nt in 'ACGT']
                for seq in window
            ]
            window_jsd = np.mean([jsd(freq_matrix[0], freq) for freq in freq_matrix[1:]])
            scores.extend([window_jsd] * self.window_size)  # Expand to nucleotide level
        return scores

    def plot_conservation(self, scores: List[float], output_file: Path) -> None:
        """Generate interactive conservation plot."""
        fig = px.line(
            x=list(range(len(scores))), 
            y=scores, 
            labels={"x": "Genome Position", "y": "Conservation Score"},
            title="GuideX Conservation Analysis"
        )
        fig.write_html(output_file)
