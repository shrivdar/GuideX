import numpy as np
from pathlib import Path
from typing import List
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.SeqRecord import SeqRecord
from skbio import DNA, TabularMSA
from scipy.spatial.distance import jensenshannon
import plotly.express as px
from .utils.logger import setup_logger

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Conservation analysis with pure Python alignment."""
    
    def __init__(self, window_size: int = 30):
        self.window_size = window_size
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'  # Global alignment strategy
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -0.5
        self.aligner.extend_gap_score = -0.1

    def align_genomes(self, genomes: List[SeqRecord], output_dir: Path = Path("alignments")) -> Path:
        """Align sequences using Biopython's pairwise aligner."""
        output_dir.mkdir(exist_ok=True, parents=True)
        
        if len(genomes) < 1:
            raise ValueError("No genomes provided for alignment")

        # Convert sequences to strings to avoid deprecation warnings
        sequences = [str(g.seq).upper() for g in genomes]
        
        # Create reference-based alignment
        alignment = MultipleSeqAlignment([])
        alignment.append(sequences[0])
        
        # Align subsequent sequences to reference
        for seq in sequences[1:]:
            aligned = self.aligner.align(alignment[0], seq)[0]
            alignment.append(aligned[1])

        # Save aligned sequences
        output_file = output_dir / "aligned.fasta"
        with open(output_file, "w") as f:
            for idx, seq in enumerate(alignment):
                f.write(f">genome_{idx}\n{seq}\n")
                
        return output_file

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
