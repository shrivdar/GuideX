import numpy as np
from pathlib import Path
from typing import List
from scikit_bio import DNA, TabularMSA
from scipy.spatial.distance import jensenshannon
import plotly.express as px
from guidex.utils.logger import setup_logger

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Conservation analysis pipeline for aligned genomes"""
    
    def __init__(self, window_size: int = 30):
        self.window_size = window_size
        self.min_conservation = 0.8

    def calculate_jsd(self, aligned_file: Path) -> List[float]:
        """Calculate Jensen-Shannon Divergence scores"""
        msa = self._load_alignment(aligned_file)
        return self._windowed_analysis(msa)

    def plot_conservation(self, scores: List[float], output_file: Path) -> None:
        """Generate interactive conservation plot"""
        fig = px.line(
            x=list(range(len(scores))),
            y=scores,
            labels={"x": "Position", "y": "Conservation"},
            title=f"Conservation Profile (Window={self.window_size})"
        )
        fig.write_html(str(output_file))

    def _load_alignment(self, path: Path) -> TabularMSA:
        """Load and validate alignment file"""
        if not path.exists():
            raise FileNotFoundError(f"Alignment file missing: {path}")
        return TabularMSA.read(path, constructor=DNA)

    def _windowed_analysis(self, msa: TabularMSA) -> List[float]:
        """Sliding window conservation scoring"""
        scores = []
        for i in range(0, len(msa[0]), self.window_size):
            window = msa[:, i:i+self.window_size]
            window_score = self._calculate_window_jsd(window)
            scores.extend([window_score] * self.window_size)
        return scores

    def _calculate_window_jsd(self, window) -> float:
        """Window-specific JSD calculation"""
        if window.shape[1] < self.window_size:
            return 0.0
            
        freq_matrix = []
        for seq in window:
            counts = {nt: seq.frequencies().get(nt, 0) for nt in 'ACGT'}
            freq_matrix.append([counts['A'], counts['C'], counts['G'], counts['T']])
            
        return np.mean([
            jensenshannon(freq_matrix[0], freq) ** 2
            for freq in freq_matrix[1:]
        ]) if len(freq_matrix) > 1 else 0.0
