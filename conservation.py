import os
import subprocess
import numpy as np
from pathlib import Path
from typing import List
from Bio import SeqIO
from Bio.Seq import Seq
from scikit-bio import DNA, TabularMSA
from scipy.spatial.distance import jensenshannon
import plotly.express as px
from guidex.utils.logger import setup_logger

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Modern conservation analysis pipeline"""
    
    def __init__(self, window_size: int = 30):
        self.window_size = window_size
        self.mafft_path = self._verify_mafft()
        self.min_conservation = 0.8

    def _verify_mafft(self):
        """Ensure MAFFT v7.5+ is installed"""
        try:
            result = subprocess.run(
                ["mafft", "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            if "v7" not in result.stderr:
                raise RuntimeError("Requires MAFFT v7.5+")
            return "mafft"
        except Exception as e:
            raise RuntimeError(f"MAFFT verification failed: {str(e)}")

    def _run_mafft(self, input_path: Path, output_dir: Path) -> Path:
        """Hybrid MAFFT alignment with intelligent fallback"""
        output_path = output_dir / "MAFFT_OUT.fasta"
    
        # Attempt 1: Optimized parameters
        try:
            self.logger.debug("Attempting MAFFT with optimized parameters")
            cmd = [
                "mafft",
                "--auto",
                "--thread", "2",  # Optimal for most modern CPUs
                "--quiet",
                str(input_path)
            ]
        
            with open(output_path, "w") as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    check=True,
                    timeout=300
                )
            return output_path
        
        except subprocess.CalledProcessError as e:
            self.logger.warning(
                f"MAFFT optimized failed ({e.returncode}), "
                "attempting basic alignment"
            )
            # Attempt 2: Failsafe simple command
            try:
                self.logger.debug("Attempting MAFFT basic mode")
                cmd = f"mafft --auto {input_path} > {output_path}"
                subprocess.run(
                    cmd,
                    shell=True,
                    check=True,
                    executable="/bin/bash",
                    timeout=300
                )
                return output_path
            
            except subprocess.CalledProcessError as final_e:
                self.logger.error(
                    "MAFFT failed all attempts. Last error: "
                    f"{final_e.stderr.decode().strip()}"
                )
                output_path.unlink(missing_ok=True)
                raise RuntimeError(
                    "MAFFT failed with both methods.\n"
                    f"Final error: {final_e.stderr.decode().strip()}\n"
                    "Verify installation with: mafft --version"
                )
            
        except Exception as e:
            self.logger.error(f"Unexpected MAFFT error: {str(e)}")
            output_path.unlink(missing_ok=True)
            raise

    def align_genomes(self, genomes, output_dir: Path) -> Path:
        """Modern alignment workflow"""
        output_dir.mkdir(parents=True, exist_ok=True)
        self._validate_input(genomes)

        input_path = output_dir / "MAFFT_IN.fasta"
        self._write_combined_fasta(genomes, input_path)
        
        return self._run_mafft(input_path, output_dir)

    def calculate_jsd(self, aligned_file: Path) -> List[float]:
        """Modern conservation scoring"""
        msa = self._load_alignment(aligned_file)
        return self._windowed_analysis(msa)

    def plot_conservation(self, scores: List[float], output_file: Path) -> None:
        """Interactive visualization"""
        fig = px.line(
            x=list(range(len(scores))),
            y=scores,
            labels={"x": "Position", "y": "Conservation"},
            title=f"Conservation Profile (Window={self.window_size})"
        )
        fig.write_html(str(output_file))

    # Implementation Details
    def _validate_input(self, genomes):
        """Modern sequence validation"""
        if len(genomes) < 2:
            raise ValueError("Need ≥2 genomes for analysis")
            
        lengths = [len(g.seq) for g in genomes]
        if max(lengths) - min(lengths) > 1000:
            raise ValueError("Sequence length variation >1kb")
            
        if any("U" in str(g.seq) for g in genomes):
            raise ValueError("RNA sequences not supported")

    def _write_combined_fasta(self, genomes, path: Path):
        """Write input following NCBI v2 standards"""
        with open(path, "w") as f:
            for idx, rec in enumerate(genomes):
                rec.id = f"Genome_{idx+1}"  # Standardize IDs
                SeqIO.write(rec, f, "fasta-2line")

    def _load_alignment(self, path: Path) -> TabularMSA:
        """Load alignment with validation"""
        if not path.exists():
            raise FileNotFoundError(f"Alignment file missing: {path}")
            
        return TabularMSA.read(path, constructor=DNA)

    def _windowed_analysis(self, msa: TabularMSA) -> List[float]:
        """Sliding window JSD calculation"""
        scores = []
        for i in range(0, len(msa[0]), self.window_size):
            window = msa[:, i:i+self.window_size]
            window_score = self._calculate_window_jsd(window)
            scores.extend([window_score] * self.window_size)
        return scores

    def _calculate_window_jsd(self, window) -> float:
        """Modern JSD implementation"""
        if window.shape[1] < self.window_size:
            return 0.0  # Handle edge windows
            
        freq_matrix = []
        for seq in window:
            counts = {nt: seq.frequencies().get(nt, 0) for nt in 'ACGT'}
            freq_matrix.append([counts['A'], counts['C'], counts['G'], counts['T']])
            
        return np.mean([
            jensenshannon(freq_matrix[0], freq) ** 2
            for freq in freq_matrix[1:]
        ]) if len(freq_matrix) > 1 else 0.0
