import os
import re
import logging
import subprocess
import shutil
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
@dataclass  # Remove frozen=True
class gRNACandidate:
    """Mutable data structure for gRNA results"""
    sequence: str
    start: int
    end: int
    gc_content: float
    mfe: float
    passes_checks: bool
    offtargets: Optional[List[dict]] = None  # New field
    offtarget_score: Optional[int] = None   # New field

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
                passes_checks=True,
                offtargets=item.get('offtargets'),        # Add these
                offtarget_score=item.get('offtarget_score')
            ) for item in raw
        ]


# Core Designer Class
class Cas13gRNADesigner:
    """High-performance Cas13 gRNA designer with parallel MFE prediction"""
    
    logger.debug("Current environment PATH: %s", os.environ.get('PATH', ''))
    logger.debug("RNAfold path: %s", '/opt/homebrew/Caskroom/miniforge/base/bin/RNAfold')
    
    def __init__(self, config_path: Path = Path("config/cas13_rules.yaml")):
        self.config = self._load_and_validate_config(config_path)
        self._verify_rnafold()

    def design(self, sequence: str, regions: List[Tuple[int, int]]) -> List[gRNACandidate]:
        """Main design pipeline with parallel processing"""
        sequence = sequence.upper()
        candidates = list(self._generate_candidates(sequence, regions))
        
        if not candidates:
            logger.warning("No candidate gRNAs generated")
            return []
        
        with ThreadPoolExecutor(self.config.max_workers) as executor:
            batches = [
                candidates[i:i+self.config.batch_size] 
                for i in range(0, len(candidates), self.config.batch_size)
            ]
            processed = []
            for batch in batches:
                processed.extend(executor.map(self._process_single_grna, batch))
        
        return [candidate for candidate in processed if candidate.passes_checks]
        unique_candidates = {c.sequence: c for c in processed if c.passes_checks}
        return list(unique_candidates.values())

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
            rnafold_path = shutil.which("RNAfold")
            
            if not rnafold_path:
                raise RuntimeError("RNAfold not found in PATH")
                
            result = subprocess.run(
                [rnafold_path, "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Handle different version formats
            version_pattern = r"(?:ViennaRNA|RNAfold)[ \t]+(\d+\.\d+\.\d+)"
            match = re.search(version_pattern, result.stdout)
            
            if not match:
                raise RuntimeError(f"Unexpected version format: {result.stdout[:50]}...")
                
            version = tuple(map(int, match.group(1).split('.')))
            logger.info(f"Detected RNAfold version: {'.'.join(map(str, version))}")
            
            if version < (2, 4, 0):
                raise RuntimeError(f"Requires ViennaRNA ≥2.4.0, found {'.'.join(map(str, version))}")
                
            self.rnafold_path = rnafold_path
            logger.info(f"RNAfold verified at: {rnafold_path}")
    
        except Exception as e:
            logger.critical(f"RNAfold verification failed: {str(e)}")
            logger.info("💡 Verify installation with: RNAfold --version")
            raise

    def _generate_candidates(
        self, 
        sequence: str, 
        regions: List[Tuple[int, int]]
    ) -> Generator[gRNACandidate, None, None]:
        """Generate candidate gRNAs from conserved regions with validation"""
        valid_chars = {'A', 'T', 'C', 'G'}
        seq_len = len(sequence)
        
        for start, end in regions:
            # Validate region boundaries
            if start < 0 or end > seq_len or start >= end:
                logger.warning(f"Skipping invalid region: ({start}, {end})")
                continue
                
            region_seq = sequence[start:end].upper()
            
            # Calculate maximum valid starting position
            max_start = len(region_seq) - self.config.spacer_length
            if max_start < 0:
                continue
                
            for i in range(max_start + 1):
                # Extract and sanitize spacer
                spacer = region_seq[i:i+self.config.spacer_length]
                spacer = ''.join([c for c in spacer if c in valid_chars])
                
                # Skip invalid spacers
                if len(spacer) != self.config.spacer_length:
                    logger.debug(f"Skipping invalid spacer: {spacer}")
                    continue
                    
                # Calculate metrics
                gc = self._calculate_gc(spacer)
                absolute_start = start + i
                absolute_end = absolute_start + self.config.spacer_length
                
                yield gRNACandidate(
                    sequence=spacer,
                    start=absolute_start,
                    end=absolute_end,
                    gc_content=gc,
                    mfe=0.0,  # Placeholder until MFE calculation
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
        """Calculate MFE using RNAfold with hybrid parsing and enhanced diagnostics"""
        try:
            # Validate input and convert to RNA format
            if not all(c in 'ATCG' for c in spacer):
                raise ValueError(f"Invalid characters in spacer: {spacer}")
            
            rna_spacer = spacer.replace('T', 'U')
            input_data = f">{spacer}\n{spacer}\n"
    
            # Execute RNAfold with full diagnostics
            process = subprocess.run(
                [self.rnafold_path, "--noPS"],
                input=input_data,
                capture_output=True,
                text=True,
                check=True,
                encoding="utf-8",
                env=os.environ
            )
            
            # Attempt 1: Precise line-based parsing
            mfe_value = None
            for line in process.stdout.split('\n'):
                if line.startswith(rna_spacer):
                    parts = line.rsplit(' ', 1)
                    if len(parts) == 2:
                        try:
                            mfe_value = float(parts[-1].strip("()"))
                            break
                        except (ValueError, IndexError):
                            logger.debug(f"Line parsing failed for: {line}")
    
            # Attempt 2: Regex fallback parsing
            if mfe_value is None:
                mfe_match = re.search(r"\(\s*([-]?\d+\.\d+)\s*\)", process.stdout)
                if mfe_match:
                    return float(mfe_match.group(1))
                    logger.debug(f"Regex fallback parsed MFE: {mfe_value}")
                else:
                    raise ValueError("No MFE found in RNAfold output")
    
            return mfe_value
    
        except subprocess.CalledProcessError as e:
            logger.error("RNAfold execution failed for spacer: %s", spacer)
            logger.error("Command: %s", " ".join(e.cmd))
            logger.error("Return code: %d", e.returncode)
            logger.error("Output:\n%s", e.stdout)
            logger.error("Error:\n%s", e.stderr)
            raise GrnaDesignError(f"RNAfold failed with code {e.returncode}")
    
        except Exception as e:
            logger.error("MFE calculation failed for spacer: %s", spacer)
            logger.error("DNA sequence: %s", spacer)
            logger.error("RNA sequence: %s", rna_spacer)
            logger.error("RNAfold input:\n%s", input_data)
            logger.error("RNAfold output:\n%s", process.stdout if 'process' in locals() else 'No output')
            raise GrnaDesignError(f"MFE calculation failed: {str(e)}")

    def _dna_to_rna(self, sequence: str) -> str:
        """Convert DNA spacer to RNA format for output matching"""
        return sequence.upper().replace('T', 'U')
    
    def _validate_rnafold_output(self, spacer: str, output: str) -> bool:
        """Verify RNAfold output contains expected RNA sequence"""
        rna_seq = self._dna_to_rna(spacer)
        return any(line.startswith(rna_seq) for line in output.split('\n'))


# Exception Classes
class GrnaDesignError(Exception):
    pass
    """Base exception for gRNA design failures"""

class InvalidSequenceError(GrnaDesignError):
    pass
    """Raised when input sequence is invalid"""

class ConfigurationError(GrnaDesignError):
    """Raised for configuration issues"""
