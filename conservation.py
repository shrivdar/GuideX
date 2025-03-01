import subprocess
import numpy as np
from skbio import DNA, TabularMSA
from scipy.spatial.distance import jensenshannon  # correct import?
import plotly.express as px
from pathlib import Path
from typing import List
from .utils.logger import setup_logger

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Identify conserved regions using Jensen-Shannon divergence."""
    
    def __init__(self, algorithm: str = "mafft", window_size: int = 30):
        self.algorithm = algorithm
        self.window_size = window_size

    def align_genomes(self, genomes: List, output_dir: Path = Path("alignments")) -> Path:
        output_dir.mkdir(exist_ok=True)
        fasta_paths = [self._save_temp(genome, output_dir) for genome in genomes]
        
        avg_length = np.mean([len(g.seq) for g in genomes])
        threads = 4 if avg_length > 1e5 else 1
        
        cmd = (
            f"mafft --thread {threads} --auto {' '.join(map(str, fasta_paths))} "
            f"> {output_dir}/aligned.fasta"
        )
        subprocess.run(cmd, shell=True, check=True)
        return output_dir / "aligned.fasta"

    def calculate_jsd(self, aligned_file: Path) -> List[float]:
        msa = TabularMSA.read(aligned_file, constructor=DNA)
        scores = []
        for i in range(0, len(msa[0]), self.window_size):
            window = msa[:, i:i+self.window_size]
            window_score = self._jsd_score(window)
            scores.extend([window_score] * self.window_size)
        return scores

    def _jsd_score(self, window):
        freq_matrix = []
        for seq in window:
            freq = [seq.frequencies().get(nt, 0) for nt in 'ACGT']
            freq_matrix.append(freq)
        
        return np.mean([
            jensenshannon(freq_matrix[0], freq) ** 2
            for freq in freq_matrix[1:]
        ])

    def plot_conservation(self, scores: List[float], output_file: Path) -> None:
        fig = px.line(
            x=list(range(len(scores))), 
            y=scores, 
            labels={"x": "Genome Position", "y": "Conservation Score"},
            title="GuideX Conservation Analysis"
        )
        fig.write_html(output_file)

    def _save_temp(self, genome, output_dir: Path) -> Path:
        path = output_dir / f"{genome.id}.fasta"
        with open(path, "w") as f:
            f.write(f">{genome.id}\n{genome.seq}")
        return path
