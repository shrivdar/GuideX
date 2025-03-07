import re
import logging
import subprocess
from pathlib import Path
from typing import List, Dict, Tuple, Generator
from concurrent.futures import ThreadPoolExecutor
import yaml
from Bio.Seq import Seq
from pydantic import BaseModel, ValidationError, confloat, conint
from dataclasses import dataclass
from datetime import datetime
from pydantic import BaseModel
import json
from typing import Optional

logger = logging.getLogger(__name__)


# Configuration Models
class Cas13Config(BaseModel):
    """Pydantic model for configuration validation"""
    spacer_length: conint(ge=20, le=30) = 28
    gc_min: confloat(ge=0.0, le=1.0) = 0.4
    gc_max: confloat(ge=0.0, le=1.0) = 0.6
    mfe_threshold: float = -2.5
    homopolymer_max: conint(ge=3, le=6) = 3
    batch_size: conint(ge=1, le=100) = 50
    max_workers: conint(ge=1, le=16) = 4


# Data Structures
@dataclass(frozen=True)
class gRNACandidate:
    """Immutable data structure for gRNA results"""
    sequence: str
    start: int
    end: int
    gc_content: float
    mfe: float
    passes_checks: bool

class gRNAResult(BaseModel):
    sequence: str
    start: int
    end: int
    gc_content: float
    mfe: float
    offtargets: Optional[List[dict]] = None
    timestamp: datetime = datetime.now()

    @classmethod
    def from_candidate(cls, candidate: gRNACandidate):
        return cls(
            sequence=candidate.sequence,
            start=candidate.start,
            end=candidate.end,
            gc_content=candidate.gc_content,
            mfe=candidate.mfe
        )

class ResultsHandler:
    """Handles serialization/deserialization of results"""
    
    @staticmethod
    def save_to_json(results: List[gRNACandidate], path: Path):
        """Save results with full metadata"""
        data = [gRNAResult.from_candidate(c).dict() for c in results]
        with open(path, 'w') as f:
            json.dump(data, f, indent=2, default=str)

    @staticmethod
    def load_from_json(path: Path) -> List[gRNACandidate]:
        """Load results with validation"""
        with open(path) as f:
            raw = json.load(f)
            
        return [
            gRNACandidate(
                sequence=item['sequence'],
                start=item['start'],
                end=item['end'],
                gc_content=item['gc_content'],
                mfe=item['mfe'],
                passes_checks=True
            ) for item in raw
        ]


# Core Designer Class
class Cas13gRNADesigner:
    """High-performance Cas13 gRNA designer with parallel MFE prediction"""
    
    def __init__(self, config_path: Path = Path("config/cas13_rules.yaml")):
        self.config = self._load_and_validate_config(config_path)
        self._verify_rnafold()

    def design(self, sequence: str, regions: List[Tuple[int, int]]) -> List[gRNACandidate]:
        """Main design pipeline with parallel processing"""
        sequence = sequence.upper()
        candidates = self._generate_candidates(sequence, regions)
        
        with ThreadPoolExecutor(self.config.max_workers) as executor:
            # Process in batches to balance memory/performance
            batches = [
                candidates[i:i+self.config.batch_size] 
                for i in range(0, len(candidates), self.config.batch_size)
            ]
            processed = []
            for batch in batches:
                processed.extend(
                    executor.map(self._process_single_grna, batch)
                )
        
        return [
            candidate for candidate in processed
            if candidate.passes_checks
        ]


    # Internal Methods
    def _load_and_validate_config(self, config_path: Path) -> Cas13Config:
        """Load and validate YAML configuration"""
        try:
            with open(config_path) as f:
                raw_config = yaml.safe_load(f)
            return Cas13Config(**raw_config)
        except (FileNotFoundError, ValidationError, yaml.YAMLError) as e:
            logger.error(f"Configuration error: {str(e)}")
            raise ValueError("Invalid configuration") from e

    def _verify_rnafold(self):
        """Ensure RNAfold is installed and accessible"""
        try:
            result = subprocess.run(
                ["RNAfold", "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            if "ViennaRNA" not in result.stdout:
                raise RuntimeError("RNAfold not properly installed")
        except Exception as e:
            logger.critical("RNAfold dependency missing")
            raise

    def _generate_candidates(self, sequence: str, regions: List[Tuple[int, int]]]) -> Generator[gRNACandidate, None, None]:
        """Generate candidate gRNAs from conserved regions"""
        for start, end in regions:
            region_seq = sequence[start:end]
            for i in range(len(region_seq) - self.config.spacer_length + 1):
                spacer = region_seq[i:i+self.config.spacer_length]
                gc = self._calculate_gc(spacer)
                
                yield gRNACandidate(
                    sequence=spacer,
                    start=start + i,
                    end=start + i + self.config.spacer_length,
                    gc_content=gc,
                    mfe=0.0,
                    passes_checks=self._validate_basic(spacer, gc)
                )

    def _process_single_grna(self, candidate: gRNACandidate) -> gRNACandidate:
        """Process individual gRNA (MFE calculation + final validation)"""
        if not candidate.passes_checks:
            return candidate
        
        try:
            mfe = self._calculate_mfe(candidate.sequence)
            passes = self._validate_final(mfe)
        except subprocess.CalledProcessError as e:
            logger.warning(f"MFE calculation failed for {candidate.sequence}: {e.stderr}")
            passes = False
            mfe = 0.0

        return gRNACandidate(
            **{**vars(candidate), "mfe": mfe, "passes_checks": passes}
        )

    def _calculate_gc(self, spacer: str) -> float:
        """Calculate GC content with vectorized operations"""
        gc = spacer.count("G") + spacer.count("C")
        return gc / len(spacer)

    def _validate_basic(self, spacer: str, gc: float) -> bool:
        """Fast pre-filter before MFE calculation"""
        return all([
            self.config.gc_min <= gc <= self.config.gc_max,
            not re.search(r"(.)\1{%d}" % self.config.homopolymer_max, spacer)
        ])

    def _validate_final(self, mfe: float) -> bool:
        """Post-MFE validation"""
        return mfe < self.config.mfe_threshold

    def _calculate_mfe(self, spacer: str) -> float:
        """Calculate MFE using RNAfold with direct stdin/stdout"""
        process = subprocess.run(
            ["RNAfold", "--noPS"],
            input=f">{spacer}\n{spacer}\n",
            capture_output=True,
            text=True,
            check=True,
            encoding="utf-8"
        )
        
        # Extract MFE from RNAfold output
        for line in process.stdout.split("\n"):
            if spacer in line:
                return float(line.split()[-1].strip("()"))
        
        raise ValueError(f"MFE parsing failed for {spacer}")


# Exception Classes
class GrnaDesignError(Exception):
    """Base exception for gRNA design failures"""

class InvalidSequenceError(GrnaDesignError):
    """Raised when input sequence is invalid"""

class ConfigurationError(GrnaDesignError):
    """Raised for configuration issues"""
