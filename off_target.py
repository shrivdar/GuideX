import shutil
import subprocess
import logging
import json
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
            # Create spacer-specific output directory
            spacer_dir = self.output_dir / f"{spacer[:8]}_{hash(spacer) % 1000:03d}"
            spacer_dir.mkdir(exist_ok=True)
            
            # Run CRISPRitz analysis
            targets = self._run_crispritz(spacer, spacer_dir)
            
            # Save and visualize results
            self._save_results(spacer, targets, spacer_dir)
            self._generate_plots(spacer, targets, spacer_dir)
            
            return targets
            
        except subprocess.TimeoutExpired:
            logger.error(f"CRISPRitz timed out for spacer: {spacer[:12]}...")
            return []

    def _run_crispritz(self, spacer: str, output_dir: Path) -> List[OffTarget]:
        """Execute CRISPRitz command with proper output handling"""
        try:
            crispritz_path = Path(__file__).parent.parent.parent / "CRISPRitz/crispritz.py"
            
            cmd = [
                "python3",
                str(crispritz_path),
                "search",
                str(self.genome_index),
                spacer,
                str(self.max_mismatches),
                "-o", str(output_dir),
                "-th", "4"
            ]
            
            subprocess.run(cmd, check=True, timeout=300)
            
            # Parse and return targets
            output_file = output_dir / f"{spacer}.targets.txt"
            if output_file.exists():
                with open(output_file) as f:
                    return self._parse_output(f.read())
            return []
            
        except subprocess.CalledProcessError as e:
            logger.error(f"CRISPRitz failed: {e.stderr}")
            raise RuntimeError("Off-target analysis failed") from e

    def _save_results(self, spacer: str, targets: List[OffTarget], output_dir: Path):
        """Save analysis results in multiple formats"""
        # Save raw CRISPRitz output
        raw_path = output_dir / "raw_results.txt"
        with open(raw_path, 'w') as f:
            f.write("\n".join(str(t) for t in targets))
            
        # Save JSON summary
        summary = {
            "spacer": spacer,
            "total_offtargets": len(targets),
            "mismatch_distribution": self._get_mismatch_counts(targets),
            "top_hits": [self._target_to_dict(t) for t in targets[:5]]
        }
        json_path = output_dir / "summary.json"
        with open(json_path, 'w') as f:
            json.dump(summary, f, indent=2)

    def _generate_plots(self, spacer: str, targets: List[OffTarget], output_dir: Path):
        """Generate visualization plots for off-target analysis"""
        plt.style.use('seaborn')
        
        # Mismatch distribution plot
        df = pd.DataFrame([self._target_to_dict(t) for t in targets])
        if not df.empty:
            plt.figure(figsize=(10, 6))
            sns.histplot(data=df, x='mismatches', bins=range(self.max_mismatches+2))
            plt.title(f"Off-target Mismatch Distribution: {spacer[:8]}...")
            plt.xlabel("Number of Mismatches")
            plt.ylabel("Count")
            plt.savefig(output_dir / "mismatch_distribution.png")
            plt.close()

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
