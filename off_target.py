import subprocess
import logging
from pathlib import Path
from typing import List, Dict
from dataclasses import dataclass
from pydantic import BaseModel, ValidationError

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
    """CRISPRitz-based off-target analysis with validation"""
    
    def __init__(self, genome_index: Path, max_mismatches: int = 3):
        self.genome_index = genome_index
        self.max_mismatches = max_mismatches
        self._validate_index()

    def analyze(self, spacer: str) -> list:
        """Run CRISPRitz analysis with proper path handling"""
        try:
            # Path to CRISPRitz script in your project structure
            crispritz_path = Path(__file__).parent.parent.parent / "CRISPRitz/crispritz.py"
            
            # Create temporary working directory
            temp_dir = Path("crispritz_temp")
            temp_dir.mkdir(exist_ok=True)
            
            cmd = [
                "python3",
                str(crispritz_path),
                "search",
                str(self.genome_index),  # Path to bowtie2 index
                spacer,
                str(self.max_mismatches),
                "-o", str(temp_dir),
                "-th", "4"  # Use 4 threads
            ]
            
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            # Parse output from generated file
            output_file = temp_dir / f"{spacer}.targets.txt"
            if output_file.exists():
                with open(output_file) as f:
                    return self._parse_output(f.read())
            return []
            
        except subprocess.TimeoutExpired:
            logger.error(f"CRISPRitz timed out for spacer: {spacer[:12]}...")
            return []
        finally:
            # Cleanup temp files
            shutil.rmtree(temp_dir, ignore_errors=True)

    def _run_crispritz(self, spacer: str) -> str:
        """Execute CRISPRitz command"""
        try:
            cmd = [
                "crispritz",
                "search",
                str(self.genome_index),
                spacer,
                str(self.max_mismatches)
            ]
            
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                encoding='utf-8'
            )
            return result.stdout
            
        except subprocess.CalledProcessError as e:
            logger.error(f"CRISPRitz failed: {e.stderr}")
            raise RuntimeError("Off-target analysis failed") from e

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
    spacer: str
    total_offtargets: int
    top_hits: List[dict]
    
    @classmethod
    def summarize(cls, targets: List[OffTarget], top_n: int = 5):
        return cls(
            spacer=targets[0].sequence if targets else "",
            total_offtargets=len(targets),
            top_hits=[
                {
                    "position": f"{t.chromosome}:{t.position}",
                    "mismatches": t.mismatches,
                    "strand": t.strand
                } for t in targets[:top_n]
            ]
        )
