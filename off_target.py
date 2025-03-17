import shutil
import subprocess
import logging
import json
import os
from pathlib import Path
from typing import List, Dict
from dataclasses import dataclass
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pydantic import BaseModel

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
    """Bowtie-based off-target analysis with enhanced error handling"""
    
    def __init__(
        self,
        genome_index: Path,
        reference_genome: Path,
        max_mismatches: int = 3,
        output_dir: Path = Path("results/off_targets"),
        threads: int = 4,
        timeout: int = 300
    ):
        self.genome_index = genome_index
        self.reference_genome = reference_genome
        self.max_mismatches = max_mismatches
        self.output_dir = output_dir
        self.threads = threads
        self.timeout = timeout
        
        self._validate_dependencies()
        self._validate_index()
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def analyze(self, spacer: str) -> List[OffTarget]:
        """Run Bowtie analysis and generate results"""
        spacer = self._normalize_spacer(spacer)
        if not self._validate_spacer(spacer):
            return []
            
        spacer_dir = self._create_output_dir(spacer)
        targets = self._run_bowtie(spacer, spacer_dir)
        
        if targets:
            self._save_results(spacer, targets, spacer_dir)
            self._generate_plots(spacer, targets, spacer_dir)
            
        self._save_metadata(spacer, targets, spacer_dir)
        return targets

    def _normalize_spacer(self, spacer: str) -> str:
        """Clean and validate spacer sequence"""
        return spacer.upper()[:28].strip()

    def _validate_spacer(self, spacer: str) -> bool:
        """Validate spacer sequence requirements"""
        valid_chars = set('ACGT')
        return (20 <= len(spacer) <= 28 and 
                all(c in valid_chars for c in spacer))

    def _create_output_dir(self, spacer: str) -> Path:
        """Create unique output directory for spacer"""
        spacer_hash = abs(hash(spacer)) % 1000
        return self.output_dir / f"{spacer[:8]}_{spacer_hash:03d}"

    def _validate_dependencies(self):
        """Verify required executables are available"""
        if not shutil.which("bowtie2"):
            raise RuntimeError("bowtie2 not found in PATH")
        if not shutil.which("samtools"):
            raise RuntimeError("samtools required for sequence extraction")
        if not self.reference_genome.exists():
            raise FileNotFoundError(f"Reference genome {self.reference_genome} not found")

    def _validate_index(self):
        """Validate Bowtie2 index files exist"""
        suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
        missing = [s for s in suffixes if not self.genome_index.with_suffix(s).exists()]
        if missing:
            raise FileNotFoundError(f"Missing Bowtie2 index files: {missing[:3]}...")

    def _run_bowtie(self, spacer: str, output_dir: Path) -> List[OffTarget]:
        """Execute Bowtie2 alignment with optimized parameters"""
        try:
            output_dir.mkdir(exist_ok=True)
            sam_path = output_dir / "alignment.sam"
            
            cmd = [
                "bowtie2",
                "-x", str(self.genome_index),
                "-U", "-",
                "--local",
                "-N", "0",       # No mismatches in seed
                "-L", "20",      # Seed length
                "-i", "C,1",     # Interval between seed indexes
                "--score-min", f"G,1,{self.max_mismatches*4}",
                "-a",            # Report all alignments
                "--sam",
                "--quiet",
                "-p", str(self.threads)
            ]

            with open(sam_path, 'w') as sam_file:
                process = subprocess.run(
                    cmd,
                    input=f">spacer\n{spacer}\n",
                    encoding='utf-8',
                    stdout=sam_file,
                    stderr=subprocess.PIPE,
                    timeout=self.timeout
                )

            if process.returncode != 0:
                raise RuntimeError(f"Bowtie failed: {process.stderr}")
                
            return self._parse_sam_output(sam_path, len(spacer))
            
        except subprocess.TimeoutExpired:
            logger.error("Bowtie analysis timed out")
            return []
        except Exception as e:
            logger.error(f"Alignment failed: {str(e)}")
            return []

    def _parse_sam_output(self, sam_path: Path, spacer_len: int) -> List[OffTarget]:
        """Parse SAM file and extract off-target information"""
        targets = []
        with open(sam_path) as f:
            for line in f:
                if line.startswith('@') or not line.strip():
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue

                flag = int(parts[1])
                chrom = parts[2]
                pos = int(parts[3])
                cigar = parts[5]
                md_tag = next((p.split(':')[-1] for p in parts if p.startswith('MD:Z:')), "")
                
                # Calculate mismatches from MD tag
                mismatches = sum(1 for c in md_tag if c.isupper())
                
                if mismatches == 0 or mismatches > self.max_mismatches:
                    continue  # Skip perfect matches and over-threshold
                    
                # Extract sequence from reference
                seq = self._extract_sequence(chrom, pos, spacer_len, '-' if (flag & 16) else '+')
                
                targets.append(OffTarget(
                    sequence=seq,
                    chromosome=chrom,
                    position=pos,
                    strand='-' if (flag & 16) else '+',
                    mismatches=mismatches
                ))
                
        return targets

    def _extract_sequence(self, chrom: str, start: int, length: int, strand: str) -> str:
        """Extract sequence using samtools faidx"""
        try:
            end = start + length - 1
            cmd = [
                "samtools", "faidx",
                str(self.reference_genome),
                f"{chrom}:{start}-{end}"
            ]
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            seq = ''.join(result.stdout.split('\n')[1:]).upper()
            return self._reverse_complement(seq) if strand == '-' else seq
        except Exception as e:
            logger.error(f"Sequence extraction failed: {str(e)}")
            return "N/A"

    def _reverse_complement(self, seq: str) -> str:
        """Generate reverse complement of DNA sequence"""
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(comp.get(base, base) for base in reversed(seq))

    def _save_results(self, spacer: str, targets: List[OffTarget], output_dir: Path):
        """Save analysis results to files"""
        # Save CSV
        df = pd.DataFrame([self._target_to_dict(t) for t in targets])
        df.to_csv(output_dir / "offtargets.csv", index=False)
        
        # Save JSON summary
        summary = {
            "spacer": spacer,
            "total_offtargets": len(targets),
            "mismatch_distribution": self._get_mismatch_counts(targets),
            "top_hits": [self._target_to_dict(t) for t in targets[:5]]
        }
        (output_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    def _generate_plots(self, spacer: str, targets: List[OffTarget], output_dir: Path):
        """Generate mismatch distribution plot"""
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
            plt.savefig(output_dir / "mismatch_distribution.png", bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            logger.error(f"Plot generation failed: {str(e)}")
        finally:
            plt.close('all')

    def _save_metadata(self, spacer: str, targets: List[OffTarget], output_dir: Path):
        """Save analysis metadata"""
        metadata = {
            "spacer": spacer,
            "valid": self._validate_spacer(spacer),
            "offtarget_count": len(targets),
            "plot_exists": (output_dir / "mismatch_distribution.png").exists(),
            "analysis_dir": str(output_dir)
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
    def summarize(cls, spacer: str, targets: List[OffTarget], output_dir: Path):
        return cls(
            spacer=spacer,
            total_offtargets=len(targets),
            analysis_directory=output_dir,
            summary_path=output_dir / "summary.json",
            plot_path=output_dir / "mismatch_distribution.png"
        )
