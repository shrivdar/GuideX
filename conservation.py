import numpy as np
from pathlib import Path
from typing import List
from Bio.Align.Applications import MafftCommandline
from Bio.SeqRecord import SeqRecord
from skbio import DNA, TabularMSA
from scipy.spatial.distance import jensenshannon
import plotly.express as px
from .utils.logger import setup_logger

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Accurate conservation analysis with MAFFT integration."""
    
    def __init__(self, window_size: int = 30):
        self.window_size = window_size

    def align_genomes(self, genomes: List[SeqRecord], output_dir: Path = Path("alignments")) -> Path:
        """MAFFT alignment using Biopython's interface (single input file)."""
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # 1. Create single input file for all sequences
        input_file = output_dir / "input.fasta"
        with open(input_file, "w") as f:
            for genome in genomes:
                clean_id = genome.id.replace("|", "_").replace(" ", "_")[:30]
                clean_seq = str(genome.seq).upper()
                f.write(f">{clean_id}\n{clean_seq}\n")

        # 2. Configure MAFFT command
        mafft_cline = MafftCommandline(
            input=str(input_file.resolve()),
            auto=True,
            thread=1,
            quiet=True
        )
        
        # 3. Execute and save output
        output_file = output_dir / "aligned.fasta"
        try:
            stdout, stderr = mafft_cline()
            with open(output_file, "w") as f:
                f.write(stdout)
            return output_file
        except Exception as e:
            logger.error(f"MAFFT failed: {stderr}")
            raise RuntimeError(f"Alignment error: {e}") from None

    def calculate_jsd(self, aligned_file: Path) -> List[float]:
        """Calculate Jensen-Shannon divergence scores."""
        msa = TabularMSA.read(aligned_file, constructor=DNA)
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
        return np.mean([jensenshannon(freq_matrix[0], freq) ** 2 for freq in freq_matrix[1:]])

    def plot_conservation(self, scores: List[float], output_file: Path) -> None:
        """Generate conservation plot."""
        fig = px.line(
            x=list(range(len(scores))), 
            y=scores,
            labels={"x": "Position", "y": "Conservation Score"},
            title="Conservation Analysis"
        )
        fig.write_html(str(output_file))
