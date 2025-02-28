import re
import yaml
from pathlib import Path
from typing import List, Dict
from Bio.Seq import Seq
import subprocess
from .utils.logger import setup_logger
from .utils.exceptions import GrnaDesignError

logger = setup_logger(__name__)

class GuideXGrnaDesigner:
    """Design Cas13 gRNAs with subtype-specific rules."""
    
    def __init__(self, subtype: str = "LwaCas13a"):
        self.subtype = subtype
        self.config = self._load_config(subtype)
        self.spacer_length = self.config["spacer_length"]

    def _load_config(self, subtype: str) -> Dict:
        """Load parameters from YAML config."""
        config_path = Path(__file__).parent.parent / "config/cas13_subtypes.yaml"
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
        return config.get(subtype, config["LwaCas13a"])  # Fallback to default

    def design(self, sequence: str, conserved_regions: List[tuple]) -> List[Dict]:
        """Design gRNAs from conserved regions."""
        grnas = []
        for start, end in conserved_regions:
            region_seq = str(sequence[start:end])
            for i in range(0, len(region_seq) - self.spacer_length + 1):
                spacer = region_seq[i:i+self.spacer_length]
                if self._is_valid(spacer):
                    grna_data = {
                        "start": start + i,
                        "end": start + i + self.spacer_length,
                        "spacer": spacer,
                        "gc": self._gc_content(spacer),
                        "mfe": self._rnafold_mfe(spacer),
                    }
                    grnas.append(grna_data)
        return grnas

    def _is_valid(self, spacer: str) -> bool:
        """Validate spacer against Cas13 rules."""
        if not (self.config["gc_min"] <= self._gc_content(spacer) <= self.config["gc_max"]):
            return False
        if re.search(r"(.)\1{3}", spacer):  # Homopolymers (4+ repeats)
            return False
        if self._rnafold_mfe(spacer) > self.config["mfe_threshold"]:
            return False
        return True

    def _gc_content(self, spacer: str) -> float:
        return (spacer.count("G") + spacer.count("C")) / len(spacer)

    def _rnafold_mfe(self, spacer: str) -> float:
        """Predict minimum free energy (MFE) using RNAfold."""
        try:
            with open("temp.fa", "w") as f:
                f.write(f">spacer\n{spacer}")
            subprocess.run(
                "RNAfold --infile=temp.fa --outfile=temp.out --noPS", 
                shell=True, 
                check=True
            )
            with open("temp.out", "r") as f:
                return float(f.readlines()[1].split()[-1].strip("()"))
        except Exception as e:
            logger.error(f"RNAfold failed: {e}")
            raise GrnaDesignError(f"RNAfold error: {e}")
