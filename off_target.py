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
    """High-performance off-target analysis using Cas-OFFinder"""
    
    def __init__(self, 
                 genome_index: Path, 
                 max_mismatches: int = 3,
                 output_dir: Path = Path("results/off_targets"),
                 cas_offinder_path: Path = Path("cas-offinder")):
        """
        genome_index: Path to reference genome FASTA
        max_mismatches: Maximum allowed mismatches
        cas_offinder_path: Path to cas-offinder executable
        """
        self.genome_index = genome_index
        self.max_mismatches = max_mismatches
        self.output_dir = output_dir
        self.cas_offinder = cas_offinder_path
        self._validate_dependencies()
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def analyze(self, spacer: str) -> List[OffTarget]:
        """Run full off-target analysis pipeline"""
        spacer = spacer.upper()
        spacer_dir = self._create_output_dir(spacer)
        
        try:
            # Generate query file
            query_file = self._create_query_file(spacer, spacer_dir)
            
            # Run Cas-OFFinder
            raw_results = spacer_dir / "raw_results.txt"
            self._run_cas_offinder(query_file, raw_results)
            
            # Parse and process results
            targets = self._parse_results(raw_results)
            
            if targets:
                self._save_results(spacer, targets, spacer_dir)
                self._generate_plots(spacer, targets, spacer_dir)
            
            return targets
            
        except subprocess.TimeoutExpired:
            logger.error(f"Analysis timed out for spacer: {spacer[:12]}...")
            return []
        finally:
            self._clean_temp_files(spacer_dir)

    def _create_output_dir(self, spacer: str) -> Path:
        """Create unique output directory for each spacer"""
        spacer_hash = abs(hash(spacer)) % 1000
        return self.output_dir / f"{spacer[:8]}_{spacer_hash:03d}"

    def _create_query_file(self, spacer: str, output_dir: Path) -> Path:
        """Generate Cas-OFFinder query file"""
        output_dir.mkdir(exist_ok=True)
        query_file = output_dir / "query.txt"
        
        # Cas-OFFinder query format: [sequence] [mismatches] [genome_path]
        query_content = f"{spacer} {self.max_mismatches} {self.genome_index}"
        query_file.write_text(query_content)
        
        return query_file

    def _run_cas_offinder(self, query_file: Path, output_file: Path):
        """Execute Cas-OFFinder with error handling"""
        cmd = [
            str(self.cas_offinder),
            str(query_file),
            "C",  # Search both strands
            str(output_file)
        ]
        
        try:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=300,
                text=True
            )
            logger.debug(f"Cas-OFFinder stdout: {result.stdout[:200]}")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Cas-OFFinder failed: {e.stderr}")
            raise RuntimeError("Off-target analysis failed") from e

    def _parse_results(self, results_file: Path) -> List[OffTarget]:
        """Parse Cas-OFFinder output file"""
        if not results_file.exists():
            return []
            
        targets = []
        with open(results_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                    
                try:
                    targets.append(OffTarget(
                        chromosome=parts[0],
                        position=int(parts[1]),
                        sequence=parts[2],
                        strand=parts[3],
                        mismatches=int(parts[4])
                    ))
                except (IndexError, ValueError) as e:
                    logger.warning(f"Invalid result line: {line}")
                    
        return targets

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
        """Generate visualization plots"""
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
        plt.savefig(output_dir / "mismatch_distribution.png")
        plt.close()

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

    def _validate_dependencies(self):
        """Verify system dependencies"""
        if not shutil.which(str(self.cas_offinder)):
            raise FileNotFoundError(
                f"Cas-OFFinder not found at {self.cas_offinder}. "
                "Install from: https://github.com/snugel/cas-offinder"
            )
            
        if not self.genome_index.exists():
            raise FileNotFoundError(f"Genome FASTA not found: {self.genome_index}")

class OffTargetResult(BaseModel):
    spacer: str
    total_offtargets: int
    analysis_directory: Path
    summary_path: Path
    plot_path: Path
    csv_path: Path

    @classmethod
    def summarize(cls, spacer: str, output_dir: Path):
        return cls(
            spacer=spacer,
            total_offtargets=len(pd.read_csv(output_dir/"offtargets.csv")),
            analysis_directory=output_dir,
            summary_path=output_dir/"summary.json",
            plot_path=output_dir/"mismatch_distribution.png",
            csv_path=output_dir/"offtargets.csv"
        )
