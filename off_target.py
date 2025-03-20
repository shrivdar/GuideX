import logging
import json
from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from pydantic import BaseModel
from tqdm import tqdm

logger = logging.getLogger(__name__)

@dataclass(frozen=True)
class OffTarget:
    sequence: str
    chromosome: str
    position: int
    strand: str
    mismatches: int
    gene_context: str = ""

class OffTargetAnalyzer:
    """Direct sequence-based off-target analysis without external dependencies"""

    def __init__(
        self,
        reference_genome: Path,
        max_mismatches: int = 3,
        output_dir: Path = Path("results/off_targets"),
        window_size: int = 28,
        batch_size: int = 1000000
    ):
        self.reference_genome = Path(reference_genome)
        self.output_dir = Path(output_dir)
        self.max_mismatches = max_mismatches
        self.window_size = window_size
        self.batch_size = batch_size

        self._validate_dependencies()
        self.genome = self._load_genome()
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def analyze(self, spacer: str) -> List[dict]:
        """Analyze off-targets with robust path handling"""
        targets = self._find_offtargets(spacer)

        # Create unique directory using MD5 hash (deterministic)
        import hashlib
        spacer_hash = hashlib.md5(spacer.encode()).hexdigest()[:8]  # Deterministic
        spacer_dir = self.output_dir / f"{spacer[:8]}_{spacer_hash}"

        # Create directory before any file operations
        spacer_dir.mkdir(parents=True, exist_ok=True)

        self._save_results(spacer, targets, spacer_dir)
        self._save_metadata(spacer, targets, spacer_dir)
        return targets

    def _normalize_spacer(self, spacer: str) -> str:
        """Sanitize and trim spacer sequence"""
        return spacer.upper()[:self.window_size].strip()

    def _validate_spacer(self, spacer: str) -> bool:
        """Validate spacer sequence format"""
        valid_chars = set('ACGT')
        return (20 <= len(spacer) <= self.window_size and
                all(c in valid_chars for c in spacer))

    def _create_output_dir(self, spacer: str) -> Path:
        """Create unique output directory per spacer"""
        spacer_hash = abs(hash(spacer)) % 1000
        return self.output_dir / f"{spacer[:8]}_{spacer_hash:03d}"

    def _validate_dependencies(self):
        """Check for required Python packages"""
        try:
            import numpy as np  # noqa: F401
            from Bio import SeqIO  # noqa: F401
        except ImportError as e:
            raise RuntimeError(f"Missing required dependency: {e}")

        if not self.reference_genome.exists():
            raise FileNotFoundError(f"Reference genome {self.reference_genome} not found")

    def _load_genome(self) -> Dict[str, str]:
        """Load genome into memory with progress tracking"""
        genome = {}
        with open(self.reference_genome) as f:
            for record in tqdm(SeqIO.parse(f, "fasta"),
                             desc="Loading genome"):
                genome[record.id] = str(record.seq).upper()
        return genome

    def _find_matches(self, search_seq: str, chromosome: str,
                    sequence: str, strand: str) -> List[OffTarget]:
        """Optimized sliding window search with vectorized operations"""
        targets = []
        seq_len = len(search_seq)
        if seq_len == 0:
            return targets

        # Convert to numpy arrays for vectorized operations
        search_arr = np.frombuffer(search_seq.encode(), dtype='S1')
        seq_arr = np.frombuffer(sequence.encode(), dtype='S1')

        # Calculate mismatches in sliding windows
        for i in tqdm(range(len(seq_arr) - seq_len + 1),
                    desc=f"Scanning {chromosome[:15]}...",
                    leave=False):
            window = seq_arr[i:i+seq_len]
            mismatches = np.sum(window != search_arr)

            if 0 < mismatches <= self.max_mismatches:
                targets.append(OffTarget(
                    sequence=window.tobytes().decode(),
                    chromosome=chromosome,
                    position=i,
                    strand=strand,
                    mismatches=int(mismatches)
                ))

        return targets

    def _reverse_complement(self, seq: str) -> str:
        """Generate reverse complement using BioPython"""
        return str(Seq(seq).reverse_complement())

    def _save_results(self, spacer: str, targets: List[dict], output_dir: Path) -> None:
        """Save results without recreating directories"""
        df = pd.DataFrame(targets)
        df.to_csv(output_dir / "offtargets.csv", index=False)
        (output_dir / "spacer.txt").write_text(spacer)

    def _generate_plots(self, spacer: str, targets: List[OffTarget],
                      output_dir: Path):
        """Generate visualization of mismatch distribution"""
        try:
            if not targets:
                return

            plt.figure(figsize=(10, 6))
            df = pd.DataFrame([self._target_to_dict(t) for t in targets])

            ax = sns.histplot(
                data=df,
                x='mismatches',
                bins=range(self.max_mismatches + 2),
                kde=True,
                edgecolor='black'
            )
            ax.set_title(f"Off-target Distribution: {spacer[:8]}...")
            plt.savefig(output_dir / "mismatch_distribution.png",
                      bbox_inches='tight')
            plt.close()

        except Exception as e:
            logger.error(f"Plot generation failed: {str(e)}")
        finally:
            plt.close('all')

    def _save_metadata(self, spacer: str, targets: List[dict], output_dir: Path) -> None:
        """Save metadata with existing directory validation"""
        metadata = {
            "spacer": spacer,
            "analysis_date": datetime.now().isoformat(),
            "offtarget_count": len(targets)
        }
        (output_dir / "metadata.json").write_text(json.dumps(metadata))

    def _get_mismatch_counts(self, targets: List[OffTarget]) -> Dict[int, int]:
        return pd.Series([t.mismatches for t in targets]).value_counts().to_dict()

    def _target_to_dict(self, target: OffTarget) -> Dict:
        return {
            "chromosome": target.chromosome,
            "position": target.position,
            "strand": target.strand,
            "mismatches": target.mismatches,
            "sequence": target.sequence
        }

class OffTargetResult(BaseModel):
    spacer: str
    total_offtargets: int
    analysis_directory: Path
    summary_path: Path
    plot_path: Path

    @classmethod
    def summarize(cls, spacer: str, targets: List[OffTarget],
                output_dir: Path):
        return cls(
            spacer=spacer,
            total_offtargets=len(targets),
            analysis_directory=output_dir,
            summary_path=output_dir / "summary.json",
            plot_path=output_dir / "mismatch_distribution.png"
        )







    def analyze(self, spacer: str) -> List[OffTarget]:
        """Main analysis workflow"""
        spacer = self._normalize_spacer(spacer)
        if not self._validate_spacer(spacer):
            return []
            
        spacer_dir = self._create_output_dir(spacer)
        targets = []
        
        # Search both forward and reverse strands
        for strand, seq in [('+', spacer), 
                          ('-', self._reverse_complement(spacer))]:
            for chrom, chrom_seq in self.genome.items():
                targets += self._find_matches(
                    search_seq=seq,
                    chromosome=chrom,
                    sequence=chrom_seq,
                    strand=strand
                )
        
        if targets:
            self._save_results(spacer, targets, spacer_dir)
            self._generate_plots(spacer, targets, spacer_dir)
            
        self._save_metadata(spacer, targets, spacer_dir)
        return targets

    def _save_results(self, spacer: str, targets: List[OffTarget], 
                    output_dir: Path):
        """Save results to CSV and JSON"""
        # CSV output
        df = pd.DataFrame([self._target_to_dict(t) for t in targets])
        df.to_csv(output_dir / "offtargets.csv", index=False)
        
        # JSON summary
        summary = {
            "spacer": spacer,
            "total_offtargets": len(targets),
            "mismatch_distribution": self._get_mismatch_counts(targets),
            "top_hits": [self._target_to_dict(t) for t in targets[:5]]
        }
        (output_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    def _save_metadata(self, spacer: str, targets: List[OffTarget], 
                     output_dir: Path):
        """Save analysis metadata"""
        metadata = {
            "spacer": spacer,
            "valid": self._validate_spacer(spacer),
            "offtarget_count": len(targets),
            "plot_exists": (output_dir / "mismatch_distribution.png").exists(),
            "analysis_dir": str(output_dir)
        }
        (output_dir / "metadata.json").write_text(json.dumps(metadata))
