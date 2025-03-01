import subprocess
import numpy as np
from pathlib import Path
from typing import List
from Bio.Align.Applications import MafftCommandline
from skbio import DNA, TabularMSA
from scipy.spatial.distance import jensenshannon
import plotly.express as px
from Bio.SeqRecord import SeqRecord
from .utils.logger import setup_logger

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Robust conservation analysis using MAFFT alignment."""
    
    def __init__(self, window_size: int = 30):
        self.window_size = window_size

    def align_genomes(self, genomes: List[SeqRecord], output_dir: Path = Path("alignments")) -> Path:
        """Alignment using Biopython's MAFFT wrapper."""
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # Create input file with all sequences
        input_file = output_dir / "input.fasta"
        with open(input_file, "w") as f:
            for genome in genomes:
                clean_id = genome.id.replace(" ", "_").replace(".", "_")[:30]
                f.write(f">{clean_id}\n{str(genome.seq).upper()}\n")

        # Use Biopython's MAFFT interface
        output_file = output_dir / "aligned.fasta"
        mafft_cline = MafftCommandline(
            input=str(input_file.resolve()),
            auto=True,
            thread=1  # Force single thread for stability
        )

        # Execute with error handling
        try:
            stdout, stderr = mafft_cline()
            with open(output_file, "w") as f:
                f.write(stdout)
            return output_file
            
        except Exception as e:
            logger.error(f"MAFFT alignment failed: {str(e)}")
            raise RuntimeError(f"Alignment failed. MAFFT error: {stderr}") from None

    def calculate_jsd(self, aligned_file: Path) -> List[float]:
        """Jensen-Shannon divergence calculation."""
        msa = TabularMSA.read(aligned_file, constructor=DNA)
        return self._sliding_jsd(msa)

    def _sliding_jsd(self, msa):
        """Calculate JSD scores across the alignment."""
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
        
        return np.mean([
            jensenshannon(freq_matrix[0], freq) ** 2 
            for freq in freq_matrix[1:]
        ])

    def plot_conservation(self, scores: List[float], output_file: Path) -> None:
        """Generate interactive conservation plot."""
        fig = px.line(
            x=list(range(len(scores))), 
            y=scores,
            labels={"x": "Position", "y": "Conservation Score"},
            title="Conservation Analysis"
        )
        fig.write_html(str(output_file))
