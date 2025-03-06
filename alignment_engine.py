import os
import subprocess
from pathlib import Path
from typing import List
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from guidex.utils.logger import setup_logger

logger = setup_logger(__name__)

class AlignmentEngine:
    """High-performance genome alignment engine"""
    
    def __init__(self, window_size: int = 30, max_threads: int = 8):
        self.window_size = window_size
        self.max_threads = max_threads
        self.aligner_path = self._verify_aligner()
        self.min_conservation = 0.8

    def _verify_aligner(self):
        """Verify MUSCLE v5+ installation"""
        try:
            result = subprocess.run(
                ["muscle", "-version"],
                capture_output=True,
                text=True,
                check=True
            )
            if "MUSCLE v5" not in result.stdout:
                raise RuntimeError("Requires MUSCLE v5+")
            return "muscle"
        except Exception as e:
            raise RuntimeError(f"Aligner verification failed: {str(e)}")

    def _run_alignment(self, input_path: Path, output_dir: Path) -> Path:
        """Optimized alignment for large datasets"""
        output_path = output_dir / "ALIGNMENT_OUT.fasta"
    
        cmd = MuscleCommandline(
            input=str(input_path.resolve()),
            out=str(output_path),
            maxiters=2,
            diags=True,
            sv=True,
            threads=self.max_threads
            # Remove invalid distance1 parameter
        )
    
        try:
            stdout, stderr = cmd()
            if stderr:  # Add error checking
                logger.error(f"MUSCLE warnings: {stderr}")
            return output_path
        except Exception as e:
            logger.error(f"Alignment failed: {stderr}")  # Use actual error output
            raise RuntimeError(f"Alignment error: {stderr}")

    def align_genomes(self, genomes: List[SeqRecord], output_dir: Path) -> Path:
        """Parallelized alignment workflow"""
        output_dir.mkdir(parents=True, exist_ok=True)
        self._validate_input(genomes)

        input_path = output_dir / "INPUT.fasta"
        self._write_chunked_fasta(genomes, input_path, chunk_size=10)
        
        return self._run_alignment(input_path, output_dir)

    def _execute_alignment(self, input_path: Path, output_dir: Path) -> Path:
        """MAFFT alignment with fallback strategies"""
        output_path = output_dir / "MAFFT_OUT.fasta"
        output_path.unlink(missing_ok=True)
    
        try:
            logger.debug("MAFFT Attempt 1: Optimized parameters")
            cmd = [
                self.mafft_path,
                "--auto",
                "--thread", str(self.max_threads),
                "--quiet",
                "--inputorder",
                str(input_path.resolve())
            ]
            
            with open(output_path, "w") as f:
                subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    check=True,
                    timeout=300,
                    env={"PATH": os.getenv("PATH")}
                )
            return output_path
        
        except subprocess.CalledProcessError as e:
            logger.warning(f"MAFFT Attempt 1 failed: {e.stderr.decode().strip()}")
            return self._fallback_alignment(input_path, output_path)

    def _fallback_alignment(self, input_path: Path, output_path: Path) -> Path:
        """Simplified alignment fallback"""
        try:
            logger.debug("MAFFT Attempt 2: Basic mode")
            subprocess.run(
                f"{self.mafft_path} --auto {input_path} > {output_path}",
                shell=True,
                check=True,
                executable="/bin/bash",
                timeout=300,
                env={"PATH": os.getenv("PATH")}
            )
            return output_path
        except subprocess.CalledProcessError as e:
            output_path.unlink(missing_ok=True)
            error_msg = f"MAFFT failed: {e.stderr.decode().strip()}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)

    def _validate_input(self, genomes: List[SeqRecord]):
        """Sequence validation"""
        if len(genomes) < 2:
            raise ValueError("Need ≥2 genomes for analysis")
            
        lengths = [len(g.seq) for g in genomes]
        if max(lengths) - min(lengths) > 1000:
            raise ValueError("Sequence length variation >1kb")
            
        if any("U" in str(g.seq) for g in genomes):
            raise ValueError("RNA sequences not supported") 

    # Optimized Helper Methods
    def _window_jsd(self, start: int) -> float:
        """Vectorized JSD calculation"""
        end = start + self.window_size
        window = self.msa[:, start:end]
        
        # Vectorized frequency calculation
        bases = np.array([list(rec.seq) for rec in window])
        counts = np.stack([
            (bases == 'A').sum(axis=1),
            (bases == 'C').sum(axis=1),
            (bases == 'G').sum(axis=1),
            (bases == 'T').sum(axis=1)
        ], axis=1)
        
        freqs = counts / counts.sum(axis=1, keepdims=True)
        return np.mean([jensenshannon(freqs[0], f)**2 for f in freqs[1:]])

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
      
        
        return output_path
