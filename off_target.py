import shutil
import subprocess
import logging
import json
import os
from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pydantic import BaseModel
from ..utils.visualization import PlotTracker

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
    """CRISPR off-target analysis with CRISPRitz/Bowtie fallback and enhanced error handling"""
    
    def __init__(
        self,
        genome_index: Path,
        reference_genome: Optional[Path] = None,  # Made optional
        max_mismatches: int = 3,
        output_dir: Path = Path("results/off_targets"),
        bowtie_path: Path = Path("bowtie"),
        threads: int = 4,
        timeout: int = 300
    ):
        self.genome_index = genome_index
        self.reference_genome = reference_genome
        self.max_mismatches = max_mismatches
        self.output_dir = output_dir
        self.bowtie_path = bowtie_path
        self.threads = threads
        self.timeout = timeout
        self.crispritz_path = Path(__file__).parent.parent.parent / "CRISPRitz/crispritz.py"
        
        self._validate_dependencies()
        self._validate_index()
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def analyze(self, spacer: str) -> list:
        """Run analysis with CRISPRitz/Bowtie and track visualization status"""
        spacer_hash = abs(hash(spacer)) % 1000
        spacer_dir = self.output_dir / f"{spacer[:8]}_{spacer_hash:03d}"
        spacer_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize plot tracking variables
        plot_generated = False
        plot_path = spacer_dir / "mismatch_distribution.png"
        targets = []
        
        try:
            logger.debug(f"🔧 Starting off-target analysis for {spacer[:12]}...")
            logger.debug(f"📂 Using output directory: {spacer_dir}")
            
            if self.crispritz_available:
                logger.debug("⚡ Using CRISPRitz engine")
                targets = self._run_crispritz(spacer, spacer_dir)
                if targets is None:
                    raise RuntimeError("CRISPRitz returned empty results")
            else:
                logger.debug("⚡ CRISPRitz unavailable, using Bowtie fallback")
                targets = self._run_bowtie(spacer, spacer_dir)
                
        except Exception as e:
            logger.error(f"🔴 Primary analysis failed: {str(e)}")
            if self.bowtie_available:
                logger.info("🔄 Attempting Bowtie fallback...")
                try:
                    targets = self._run_bowtie(spacer, spacer_dir)
                except Exception as bowtie_error:
                    logger.critical(f"🔥 Complete analysis failure: {bowtie_error}")
                    targets = []
            else:
                logger.critical("❌ No available analysis engines")
                targets = []
    
        try:
            if targets:
                logger.debug(f"📦 Saving {len(targets)} results...")
                self._save_results(spacer, targets, spacer_dir)
                
                logger.debug("🎨 Generating visualization...")
                self._generate_plots(spacer, targets, spacer_dir)
                plot_generated = plot_path.exists()
                
                # Add debug output for plot paths
                if plot_generated:
                    logger.debug(f"🖼️ Plot saved to: {plot_path}")
                else:
                    logger.warning("⚠️ Plot file missing despite generation attempt")
        except Exception as save_error:
            logger.error(f"📦 Results handling failed: {save_error}")
            plot_generated = False
    
        # Final status report
        status_msg = [
            f"\n🔍 OFF-TARGET ANALYSIS REPORT: {spacer[:8]}...",
            f"▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔",
            f"• Total off-targets: {len(targets)}",
            f"• Output directory: {spacer_dir}",
            f"• Visualization: {'SUCCESS' if plot_generated else 'FAILED'}",
        ]
        
        if plot_generated:
            status_msg.append(f"• Plot path: {plot_path}")
        
        logger.info("\n".join(status_msg))
        
        return targets

    def _validate_dependencies(self):
        """Check availability of required tools"""
        self.crispritz_available = self.crispritz_path.exists()
        self.bowtie_available = False

        if self.crispritz_available:
            try:
                subprocess.run(
                    [str(self.crispritz_path), "--help"],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )
                logger.info("CRISPRitz available and functional")
                return
            except Exception as e:
                logger.warning(f"CRISPRitz check failed: {str(e)}")

        # Bowtie fallback checks
        self.crispritz_available = False
        if not shutil.which(str(self.bowtie_path)):
            raise RuntimeError(f"Bowtie not found at {self.bowtie_path}")
        if not shutil.which("samtools"):
            raise RuntimeError("samtools required for Bowtie fallback not found in PATH")
        if not self.reference_genome:
            raise ValueError("reference_genome must be provided for Bowtie fallback")
        if not self.reference_genome.exists():
            raise FileNotFoundError(f"Reference genome {self.reference_genome} not found")
        
        self.bowtie_available = True
        logger.info("Using Bowtie fallback configuration")

    def _validate_index(self):
        """Validate genome index exists"""
        suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
        missing = [s for s in suffixes if not self.genome_index.with_suffix(s).exists()]
        if missing:
            raise FileNotFoundError(f"Missing index files: {missing[:3]}...")

    def _run_crispritz(self, spacer: str, output_dir: Path) -> Optional[List[OffTarget]]:
        try:
            # Add spacer validation
            if not self._is_valid_spacer(spacer):
                logger.error(f"Invalid spacer: {spacer} (length={len(spacer)})")
                return []
                
            # Add debug logging for raw command
            logger.debug(f"CRISPRitz command: {' '.join(cmd)}")
            
            # Run with output capture for debugging
            result = subprocess.run(
                cmd,
                check=True,
                timeout=self.timeout,
                capture_output=True,
                text=True
            )
            logger.debug(f"CRISPRitz raw output:\n{result.stdout[:200]}...")  # First 200 chars
            
            # Parse output
            return self._parse_crispritz_output(output_file)
            
        except subprocess.CalledProcessError as e:
            logger.error(f"CRISPRitz failed with code {e.returncode}")
            logger.debug(f"STDERR: {e.stderr}")  # Critical debug info
            return None
    
    def _is_valid_spacer(self, spacer: str) -> bool:
        """Validate spacer sequence before analysis"""
        return len(spacer) >= 20 and all(c in 'ACGT' for c in spacer.upper())

    def _run_bowtie(self, spacer: str, output_dir: Path) -> List[OffTarget]:
        """Execute Bowtie alignment with sequence extraction"""
        try:
            sam_path = output_dir / "bowtie_results.sam"
            cmd = [
                str(self.bowtie_path),
                "-v", str(self.max_mismatches),
                "-a",
                "--best",
                "--strata",
                "--sam",
                "-p", str(self.threads),
                str(self.genome_index),
                "-c", spacer,
                str(sam_path)
            ]
            
            subprocess.run(
                cmd,
                check=True,
                timeout=self.timeout,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
            
            return self._parse_bowtie_output(sam_path, len(spacer))
        
        except Exception as e:
            logger.error(f"Bowtie failed: {str(e)}")
            return []

    def _parse_crispritz_output(self, output_file: Path) -> List[OffTarget]:
        """Parse CRISPRitz output file"""
        with open(output_file) as f:
            return [
                OffTarget(
                    sequence=parts[4],
                    chromosome=parts[0],
                    position=int(parts[1]),
                    strand=parts[2],
                    mismatches=int(parts[3])
                )
                for line in f if not line.startswith('#')
                if len(parts := line.strip().split('\t')) >= 5
            ]

    def _parse_bowtie_output(self, sam_path: Path, spacer_len: int) -> List[OffTarget]:
        """Parse Bowtie SAM output and extract sequences"""
        targets = []
        with open(sam_path) as f:
            for line in f:
                if line.startswith('@'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue

                flag = int(parts[1])
                chrom = parts[2]
                pos = int(parts[3])
                strand = '-' if (flag & 16) else '+'
                
                # Calculate mismatches from MD tag
                md_tag = next((p for p in parts if p.startswith('MD:Z:')), None)
                mismatches = sum(1 for c in md_tag.split(':')[-1] if c.isupper()) if md_tag else 0
                
                # Extract sequence from reference genome
                seq = self._extract_sequence(chrom, pos, spacer_len, strand)
                
                targets.append(OffTarget(
                    sequence=seq,
                    chromosome=chrom,
                    position=pos,
                    strand=strand,
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
        """Save analysis artifacts"""
        # Save JSON summary
        summary = {
            "spacer": spacer,
            "total_offtargets": len(targets),
            "top_hits": [self._target_to_dict(t) for t in targets[:5]],
            "mismatch_distribution": self._get_mismatch_counts(targets)
        }
        (output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
        
        # Save raw CSV
        pd.DataFrame([self._target_to_dict(t) for t in targets]).to_csv(
            output_dir / "offtargets.csv", index=False
        )

    def _generate_plots(self, spacer: str, targets: List[OffTarget], output_dir: Path):
        """Generate visualization plots with path reporting"""
        try:
            if not targets:
                return
                
            sns.set_theme(style="whitegrid")
            df = pd.DataFrame([self._target_to_dict(t) for t in targets])
            
            plt.figure(figsize=(10, 6))
            ax = sns.histplot(
                data=df,
                x='mismatches',
                bins=range(self.max_mismatches + 2),
                kde=True,
                edgecolor='black'
            )
            ax.set_title(f"Off-target Distribution: {spacer[:8]}...")
            plot_path = output_dir / "mismatch_distribution.png"
            plt.savefig(plot_path)
            plt.close()
            logger.info(f"Successfully generated plot at: {plot_path}")
            
        except Exception as e:
            logger.error(f"Plot generation failed: {str(e)}")
            raise

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

    def _clean_temp_files(self, output_dir: Path):
        """Remove intermediate files"""
        for f in output_dir.glob("*.tmp"):
            try:
                f.unlink()
            except Exception as e:
                logger.debug(f"Failed to delete {f}: {str(e)}")


class OffTargetResult(BaseModel):
    spacer: str
    total_offtargets: int
    analysis_directory: Path
    summary_path: Path
    plot_path: Path
    raw_results_path: Path

    @property
    def plot_generated(self) -> bool:
        return self.plot_path.exists()

    @classmethod
    def summarize(cls, spacer: str, targets: List[OffTarget], output_dir: Path):
        result = cls(
            spacer=spacer,
            total_offtargets=len(targets),
            analysis_directory=output_dir,
            summary_path=output_dir / "summary.json",
            plot_path=output_dir / "mismatch_distribution.png",
            raw_results_path=output_dir / "raw_results.txt"
        )
        logger.debug(f"Plot generation status: {'Success' if result.plot_generated else 'Failed'}")
        return result
