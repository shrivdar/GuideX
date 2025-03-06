import numpy as np
from pathlib import Path
from typing import List
from skbio import DNA, TabularMSA
from scipy.spatial.distance import jensenshannon
import plotly.express as px
from guidex.utils.logger import setup_logger

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Conservation analysis pipeline for aligned genomes"""
    
    def __init__(self, window_size: int = 30):
        self.window_size = window_size
        self.min_conservation = 0.7

    def calculate_jsd(self, aligned_file: Path) -> List[float]:
        """Calculate Jensen-Shannon Divergence scores for aligned sequences"""
        # Validate alignment file
        if not aligned_file.exists():
            raise FileNotFoundError(f"Alignment file missing: {aligned_file}")
        
        try:
            # Load and validate alignment
            msa = self._load_alignment(aligned_file)
            if len(msa) < 2:
                raise ValueError("At least 2 sequences required for conservation analysis")
            
            # Calculate conservation scores
            logger.info(f"Analyzing conservation for {len(msa)} sequences")
            return self._windowed_analysis(msa)
            
        except Exception as e:
            logger.error(f"Conservation analysis failed: {str(e)}")
            raise

    def plot_conservation(self, scores: List[float], output_file: Path) -> None:
        """Generate interactive conservation plot"""
        try:
            fig = px.line(
                x=list(range(len(scores))[:len(scores)//self.window_size*self.window_size],
                y=scores,
                labels={"x": "Position", "y": "Conservation"},
                title=f"Conservation Profile (Window={self.window_size})"
            )
            output_file.parent.mkdir(parents=True, exist_ok=True)
            fig.write_html(str(output_file))
            logger.info(f"Saved conservation plot to {output_file}")
        except Exception as e:
            logger.error(f"Failed to generate conservation plot: {str(e)}")
            raise

    def _load_alignment(self, path: Path) -> TabularMSA:
        """Load and validate alignment file"""
        logger.debug(f"Loading alignment from {path}")
        msa = TabularMSA.read(str(path), constructor=DNA)
        
        # Validate alignment content
        if len(msa) == 0:
            raise ValueError("Empty alignment file")
        if any(len(seq) != len(msa[0]) for seq in msa):
            raise ValueError("Inconsistent sequence lengths in alignment")
            
        return msa

    def _windowed_analysis(self, msa: TabularMSA) -> List[float]:
        """Sliding window conservation scoring"""
        scores = []
        seq_length = len(msa[0])
        
        for i in range(0, seq_length - self.window_size + 1):
            window_scores = []
            for j in range(i, min(i+self.window_size, seq_length)):
                position_scores = []
                for seq1 in msa:
                    for seq2 in msa:
                        if seq1 != seq2:
                            p = self._position_frequencies(seq1[j])
                            q = self._position_frequencies(seq2[j])
                            position_scores.append(jensenshannon(p, q) ** 2)
                window_scores.append(np.mean(position_scores) if position_scores else 0.0)
            scores.extend(window_scores)
            
        return scores

    def _position_frequencies(self, nucleotide: str) -> List[float]:
        """Calculate normalized nucleotide frequencies at a position"""
        valid_nt = {'A', 'C', 'G', 'T'}
        counts = {
            'A': 0.25,  # Pseudocounts to avoid zero probabilities
            'C': 0.25,
            'G': 0.25,
            'T': 0.25
        }
        
        if nucleotide.upper() in valid_nt:
            counts[nucleotide.upper()] += 1.0
            
        total = sum(counts.values())
        return [v/total for v in counts.values()]
