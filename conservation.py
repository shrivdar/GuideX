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
        """Safe MAFFT alignment using Biopython's interface."""
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # 1. Create single input file with all sequences
        input_file = output_dir / "input.fasta"
        with open(input_file, "w") as f:
            for genome in genomes:
                clean_id = genome.id.replace("|", "_").replace(" ", "_")[:30]
                f.write(f">{clean_id}\n{str(genome.seq).upper()}\n")

        # 2. Configure MAFFT command using Biopython
        mafft_cline = MafftCommandline(
            input=str(input_file.resolve()),
            auto=True,
            thread=1,
            quiet=True
        )
        
        # 3. Execute and handle output
        output_file = output_dir / "aligned.fasta"
        try:
            stdout, stderr = mafft_cline()
            with open(output_file, "w") as f:
                f.write(stdout)
            return output_file
        except Exception as e:
            logger.error(f"MAFFT failed: {stderr}")
            raise RuntimeError(f"Alignment error: {e}") from None

    # Keep existing calculate_jsd, _window_jsd, plot_conservation methods
    # ... (unchanged from previous implementation)
