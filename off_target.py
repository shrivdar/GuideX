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
    """CRISPRitz-based off-target analysis with validation and enhanced output handling"""
    
    def __init__(self, genome_index: Path, max_mismatches: int = 3, output_dir: Path = Path("results/off_targets")):
        self.genome_index = genome_index
        self.max_mismatches = max_mismatches
        self.output_dir = output_dir
        self._validate_index()
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def analyze(self, spacer: str) -> list:
        """Run CRISPRitz analysis and save full results with visualization"""
        try:
            # Create unique output directory with spacer hash
            spacer_hash = abs(hash(spacer)) % 1000
            spacer_dir = self.output_dir / f"{spacer[:8]}_{spacer_hash:03d}"
            spacer_dir.mkdir(exist_ok=True)
            
            # Run analysis with strict output handling
            targets = self._run_crispritz(spacer, spacer_dir)
            
            if targets:
                self._save_results(spacer, targets, spacer_dir)
                self._generate_plots(spacer, targets, spacer_dir)
            
            return targets
            
        except subprocess.TimeoutExpired:
            logger.error(f"CRISPRitz timed out for spacer: {spacer[:12]}...")
            return []
        finally:
            self._clean_temp_files(spacer_dir)

    def _run_crispritz(self, spacer: str, output_dir: Path) -> List[OffTarget]:
        """Execute CRISPRitz command with improved error handling"""
        try:
            crispritz_path = Path(__file__).parent.parent.parent / "CRISPRitz/crispritz.py"
            if not crispritz_path.exists():
                raise FileNotFoundError(f"CRISPRitz not found at {crispritz_path}")

            output_file = output_dir / "crispritz_results.txt"
            
            cmd = [
                "python3",
                str(crispritz_path),
                "search",
                str(self.genome_index),
                spacer,
                str(self.max_mismatches),
                "-o", str(output_file.with_suffix('')),  # CRISPRitz auto-adds extension
                "-th", "4",
                "-quiet"  # Suppress interactive prompts
            ]
            
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=300,
                text=True,
                input='\n'  # Force feed newline to prevent hangs
            )
            
            # Debug logging
            logger.debug(f"CRISPRitz stdout: {result.stdout[:200]}")
            if result.stderr:
                logger.warning(f"CRISPRitz stderr: {result.stderr[:200]}")

            # Handle output file naming variations
            for ext in ['.targets.txt', '.results.txt']:
                candidate_file = output_file.with_suffix(ext)
                if candidate_file.exists():
                    with open(candidate_file) as f:
                        return self._parse_output(f.read())
            
            logger.warning(f"No output file found for spacer: {spacer}")
            return []

        except subprocess.CalledProcessError as e:
            logger.error(f"CRISPRitz failed with code {e.returncode}: {e.stderr[:200]}")
            raise RuntimeError("Off-target analysis failed") from e

    def _save_results(self, spacer: str, targets: List[OffTarget], output_dir: Path):
        """Save analysis results with validation"""
        try:
            raw_path = output_dir / "raw_results.txt"
            with open(raw_path, 'w') as f:
                f.write("\n".join(f"{t.chromosome}\t{t.position}\t{t.strand}\t{t.mismatches}\t{t.sequence}" 
                          for t in targets))
            
            summary = {
                "spacer": spacer,
                "total_offtargets": len(targets),
                "mismatch_distribution": self._get_mismatch_counts(targets),
                "top_hits": [self._target_to_dict(t) for t in targets[:5]]
            }
            (output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
            
        except Exception as e:
            logger.error(f"Failed to save results: {str(e)}")

    def _generate_plots(self, spacer: str, targets: List[OffTarget], output_dir: Path):
        """Generate visualization plots with error handling"""
        try:
            if not targets:
                return
                
            sns.set_theme(style="whitegrid", palette="muted")
            df = pd.DataFrame([self._target_to_dict(t) for t in targets])
            
            plt.figure(figsize=(10, 6))
            ax = sns.histplot(
                data=df, 
                x='mismatches', 
                bins=range(self.max_mismatches + 2),
                kde=True,
                edgecolor='black'
            )
            ax.set_title(f"Off-target Mismatch Distribution: {spacer[:8]}...")
            ax.set_xlabel("Number of Mismatches")
            ax.set_ylabel("Count")
            plt.tight_layout()
            plt.savefig(output_dir / "mismatch_distribution.png")
            plt.close()
            
        except Exception as e:
            logger.error(f"Plot generation failed: {str(e)}")

    def _clean_temp_files(self, output_dir: Path):
        """Clean up intermediate files"""
        try:
            for f in output_dir.glob("*.tmp"):
                f.unlink()
        except Exception as e:
            logger.warning(f"Temp file cleanup failed: {str(e)}")

    def _get_mismatch_counts(self, targets: List[OffTarget]) -> Dict[int, int]:
        counts = {}
        for t in targets:
            counts[t.mismatches] = counts.get(t.mismatches, 0) + 1
        return counts

    def _target_to_dict(self, target: OffTarget) -> Dict:
        return {
            "chromosome": target.chromosome,
            "position": target.position,
            "strand": target.strand,
            "mismatches": target.mismatches,
            "sequence": target.sequence
        }

    def _parse_output(self, output: str) -> List[OffTarget]:
        """Parse and validate CRISPRitz output"""
        targets = []
        for line in output.split('\n'):
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            try:
                targets.append(OffTarget(
                    sequence=parts[4],
                    chromosome=parts[0],
                    position=int(parts[1]),
                    strand=parts[2],
                    mismatches=int(parts[3])
                ))
            except (IndexError, ValueError) as e:
                logger.warning(f"Invalid line format: {line}")
                
        return targets

    def _validate_index(self):
        """Verify genome index exists"""
        required_files = [
            self.genome_index.with_suffix(suf) 
            for suf in ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
        ]
        
        missing = [f for f in required_files if not f.exists()]
        if missing:
            raise FileNotFoundError(
                f"Missing genome index files: {missing[:3]}..."
            )

class OffTargetResult(BaseModel):
    # Add output paths to results
    spacer: str
    total_offtargets: int
    analysis_directory: Path
    summary_path: Path
    plot_path: Path
    raw_results_path: Path
    
    @classmethod
    def summarize(cls, spacer: str, targets: List[OffTarget], output_dir: Path):
        return cls(
            spacer=spacer,
            total_offtargets=len(targets),
            analysis_directory=output_dir,
            summary_path=output_dir / "summary.json",
            plot_path=output_dir / "mismatch_distribution.png",
            raw_results_path=output_dir / "raw_results.txt"
        )
