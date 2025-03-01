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
        """Robust MAFFT alignment with comprehensive error handling."""
        output_dir.mkdir(exist_ok=True, parents=True)
    
        # 1. Generate safe filenames and validate sequences
        fasta_paths = []
        for idx, genome in enumerate(genomes):
            # Sanitize filename and sequence
            clean_id = genome.id.replace(" ", "_").replace(".", "_").replace("/", "-")
            clean_seq = str(genome.seq).replace(" ", "").upper()
        
            if not clean_seq:
                raise ValueError(f"Empty sequence for genome {clean_id}")
                
            path = output_dir / f"input_{idx}_{clean_id[:30]}.fasta"
            with open(path, "w") as f:
                f.write(f">{clean_id}\n{clean_seq}")
            fasta_paths.append(path)

        # 2. Verify input files
        missing_files = [str(fp) for fp in fasta_paths if not fp.exists()]
        if missing_files:
            raise FileNotFoundError(f"Missing input files: {', '.join(missing_files)}")

        # 3. Build MAFFT command
        output_path = output_dir / "aligned.fasta"
        cmd = [
            "mafft",
            "--auto",
            "--quiet",
            "--thread", str(min(os.cpu_count(), 4)),  # Optimal thread count
            "--out", str(output_path.resolve())  # Absolute path for safety
        ] + [str(fp.resolve()) for fp in fasta_paths]  # Absolute paths for inputs

        # 4. Execute with robust error handling
        try:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
            )
            if "error" in result.stderr.lower():
                raise RuntimeError(f"MAFFT warning: {result.stderr}")
            
        except subprocess.CalledProcessError as e:
            error_msg = (
                f"MAFFT failed with code {e.returncode}\n"
                f"Command: {' '.join(cmd)}\n"
                f"Error: {e.stderr}"
            )
            raise RuntimeError(error_msg) from None

        return output_path
    
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
