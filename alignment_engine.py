import os
import subprocess
import logging
from pathlib import Path
from typing import List
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from guidex.utils.logger import setup_logger

logger = setup_logger(__name__)

class AlignmentEngine:
    """Memory-optimized alignment engine for large plant genomes"""
    
    def __init__(self, max_threads: int = 8):
        self.max_threads = max_threads
        self.max_input_size = 50_000_000  # 50MB
        self.aligner_path = self._verify_aligner()
        self._setup_parameters()

    def _setup_parameters(self):
        """Configure MUSCLE parameters for cassava genomes"""
        self.muscle_params = [
            "-maxiters", "1",       # Reduced iterations for speed
            "-diags",               # Faster for similar sequences
            "-sv",                  # Use sparse vertical optimization
            "-quiet",
            "-threads", str(self.max_threads)
        ]

    def _verify_aligner(self) -> str:
        """Verify MUSCLE v5+ installation with large sequence support"""
        try:
            result = subprocess.run(
                ["muscle", "-version"],
                capture_output=True,
                text=True,
                check=True
            )
            if "muscle 5" not in result.stdout:
                raise RuntimeError("MUSCLE v5+ required for plant genome alignment")
            return "muscle"
        except Exception as e:
            logger.error(f"Aligner verification failed: {str(e)}")
            raise

    def align(self, genomes: List[SeqRecord], output_dir: Path) -> Path:
        """Safe alignment with memory constraints"""
        output_dir.mkdir(parents=True, exist_ok=True)
        self._validate_input(genomes)
    
        input_path = output_dir / "INPUT.fasta"
        self._write_chunked_fasta(genomes, input_path)
        
        # Check memory requirements before alignment
        total_size = sum(len(g.seq) for g in genomes)
        if total_size > self.max_input_size:
            raise MemoryError(
                f"Input too large ({total_size/1e6:.1f}MB). "
                "Consider subsetting or using progressive alignment."
            )

        return self._run_alignment(input_path, output_dir)

    def _run_alignment(self, input_file: Path, output_dir: Path) -> Path:
        """Core alignment logic with optimized parameters"""
        output_file = output_dir / "ALIGNMENT_OUT.fasta"
        
        try:
            command = ["muscle", "-align", str(input_file), 
                      "-output", str(output_file)] + self.muscle_params
            logger.info(f"Executing: {' '.join(command)}")
            
            subprocess.run(
                command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            return output_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Alignment failed: {e.stderr.decode()}")
            raise RuntimeError("MUSCLE alignment crashed") from e

    def _validate_input(self, genomes: List[SeqRecord]):
        """Enhanced validation for plant genome characteristics"""
        if len(genomes) < 2:
            raise ValueError("Need ≥2 genomes for analysis")
            
        lengths = [len(g.seq) for g in genomes]
        if max(lengths) - min(lengths) > 1500:
            raise ValueError("Sequence length variation >1.5kb - check genome versions")
            
        n_count = sum(str(g.seq).upper().count('N') for g in genomes)
        if n_count > 1000:
            raise ValueError(f"Excessive ambiguous bases ({n_count} Ns)")

    def _write_chunked_fasta(self, genomes: List[SeqRecord], path: Path, chunk_size: int = 5):
        """Memory-safe FASTA writing for large sequences"""
        try:
            with open(path, "w") as f:
                for i in range(0, len(genomes), chunk_size):
                    chunk = genomes[i:i+chunk_size]
                    for idx, rec in enumerate(chunk):
                        sanitized_seq = str(rec.seq).replace('?', 'N')  # Handle ambiguous chars
                        clean_rec = SeqRecord(
                            Seq(sanitized_seq),
                            id=f"Cassava_{i+idx+1}",
                            description=""
                        )
                        SeqIO.write(clean_rec, f, "fasta-2line")
            logger.info(f"Wrote {len(genomes)} sequences to {path}")
        except IOError as e:
            logger.error(f"Failed to write input file: {str(e)}")
            raise
