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
    """High-performance genome alignment engine using MUSCLE"""
    
    def __init__(self, max_threads: int = 8):
        self.max_threads = max_threads
        self.aligner_path = self._verify_aligner()

    def _verify_aligner(self) -> str:
        """Verify MUSCLE v5+ installation"""
        try:
            result = subprocess.run(
                ["muscle", "-version"],
                capture_output=True,
                text=True,
                check=True
            )
            # Check if output contains "muscle 5.x.x"
            if "muscle 5" in result.stdout:
                return "muscle"
            else:
                raise RuntimeError("Requires MUSCLE v5+")
        except Exception as e:
            logger.error(f"Aligner verification failed: {str(e)}")
            raise

    def align(self, genomes: List[SeqRecord], output_dir: Path) -> Path:
        """Optimized alignment workflow for large datasets"""
        output_dir.mkdir(parents=True, exist_ok=True)
        self._validate_input(genomes)

        input_path = output_dir / "INPUT.fasta"
        self._write_chunked_fasta(genomes, input_path, chunk_size=10)
        
        return self._run_alignment(input_path, output_dir)

    def _run_alignment(self, input_path: Path, output_dir: Path) -> Path:
        """Execute MUSCLE alignment with optimized parameters"""
        output_path = output_dir / "ALIGNMENT_OUT.fasta"
        
        cmd = [
            "muscle",
            "-align", str(input_path),
            "-output", str(output_path),
            "-maxiters", "2",
            "-diags",
            "-sv",
            "-threads", str(self.max_threads)
        ]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            if result.stderr:
                logger.warning(f"MUSCLE warnings: {result.stderr}")
            return output_path
        except subprocess.CalledProcessError as e:
            logger.error(f"Alignment failed: {str(e)}")
            raise RuntimeError(f"MUSCLE error: {str(e)}")

    def _validate_input(self, genomes: List[SeqRecord]):
        """Validate input genome sequences"""
        if len(genomes) < 2:
            raise ValueError("Need ≥2 genomes for analysis")
            
        lengths = [len(g.seq) for g in genomes]
        if max(lengths) - min(lengths) > 1000:
            raise ValueError("Sequence length variation >1kb")
            
        if any("U" in str(g.seq) for g in genomes):
            raise ValueError("RNA sequences not supported")

    def _write_chunked_fasta(self, genomes: List[SeqRecord], path: Path, chunk_size: int = 10):
        """Write input in chunks for memory efficiency"""
        with open(path, "w") as f:
            for i in range(0, len(genomes), chunk_size):
                chunk = genomes[i:i+chunk_size]
                for idx, rec in enumerate(chunk):
                    rec.id = f"Genome_{i+idx+1}"
                    SeqIO.write(rec, f, "fasta-2line")
