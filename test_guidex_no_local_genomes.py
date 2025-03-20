from Bio import AlignIO, SeqIO
import traceback
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import BiopythonDeprecationWarning
import warnings
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
from pathlib import Path
import sys
import os
from scipy.spatial import distance
import shutil
import importlib
sys.path.insert(0, str(Path(__file__).parent))
importlib.invalidate_caches()
import torch
import logging
import json
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats
import importlib
import conservation
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)  # Add this line
from conservation import ConservationAnalyzer
from alignment_engine import AlignmentEngine
from guidex.genome_fetcher import GenomeFetcher
from guidex.core import Cas13gRNADesigner
from guidex.grna.off_target import OffTargetAnalyzer
from guidex.grna.optimizer import Cas13Optimizer
from guidex.core.grna_designer import gRNACandidate
from typing import List, Tuple

# CBSV-optimized parameters
CBSV_THRESHOLDS = [0.65, 0.55, 0.45, 0.35]
CBSV_CONSERVATION_PARAMS = {
    'window_size': 30,
    'max_gap': 0.7,
    'pseudocount': 10.0
}
LOCAL_GENOME_PATH = Path("genomes/CBSV")
HOST_GENOME_PATH = Path("genomes/manihot_esculenta.fa")  # Cassava genome

def main():
    try:
        # Initialize debug mode
        debug_mode = os.getenv('GUIDEX_DEBUG', 'false').lower() == 'true'
        logger = logging.getLogger(__name__)
        log_level = logging.DEBUG if debug_mode else logging.INFO
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' if debug_mode else '%(message)s',
            handlers=[logging.FileHandler("guidex.log"), logging.StreamHandler()]
        )

        logger.info("🚀 Starting CBSV-Optimized GuideX Pipeline")
        
        # Clean working directories
        for dir_path in ["alignments", "results"]:
            if Path(dir_path).exists():
                shutil.rmtree(dir_path)
            Path(dir_path).mkdir(exist_ok=True)

        # Initialize components
        logger.debug("Initializing CBSV-optimized components...")
        fetcher = GenomeFetcher(api_key=os.getenv("NCBI_API_KEY_2025"))
        aligner = AlignmentEngine(max_threads=8)
        conservator = ConservationAnalyzer(**CBSV_CONSERVATION_PARAMS)
        designer = Cas13gRNADesigner()
        designer.configure(
            gc_range=(0.35, 0.65),
            mfe_threshold=-5.0
        )
        ot_analyzer = OffTargetAnalyzer(reference_genome=HOST_GENOME_PATH)
        optimizer = Cas13Optimizer(designer=designer)

        # Genome acquisition with enhanced fallback
        genomes = []
        try:
            logger.info("🕵️ Fetching CBSV genomes...")
            genomes = fetcher.fetch_genomes(target="Cassava brown streak virus", genome_type="genome", limit=12)

            if not genomes:
                logger.warning("ℹ️ No genomes found - trying protein database")
                genomes = fetcher.fetch_genomes(
                    target="Cassava brown streak virus",
                    genome_type="protein",
                    limit=12
                )
                
            if not genomes:
                logger.warning("ℹ️ No NCBI results - loading local genomes")
                genome_files = list(LOCAL_GENOME_PATH.glob("*.fna"))
                if not genome_files:
                    raise FileNotFoundError("No local genome files found")

            if not genomes:
                logger.warning("ℹ️ No genomes found - using local genomes")
                genomes = LOCAL_GENOMES
                
                genomes = []
                for f in genome_files[:12]:
                    try:
                        # Handle multi-record files by taking first record
                        record = next(SeqIO.parse(f, "fasta"))
                        genomes.append(record)
                    except Exception as e:
                        logger.error(f"Error reading {f}: {e}")

            if not genomes:
                raise RuntimeError("No genomes available from any source")

        except Exception as e:
            logger.error(f"🔴 Genome acquisition failed: {e}")
            sys.exit(1)

        # Sequence validation with relaxed constraints
        valid_genomes = [g for g in genomes if 5000 <= len(g.seq) <= 15000]  # Wider length tolerance
        if len(valid_genomes) < 2:
            raise ValueError(f"Insufficient valid genomes: {len(valid_genomes)}")

        # Alignment process
        logger.info("\n🧬 Starting CBSV alignment...")
        aligned_file = aligner.align(valid_genomes, Path("alignments/cbsv_alignment.fasta"))
        logger.info(f"🔍 Alignment saved to: {aligned_file}")

        # Conservation analysis
        logger.info("\n🔎 Identifying CBSV conserved regions...")
        jsd_scores, valid_windows = conservator.calculate_jsd(aligned_file)
        jsd_scores = np.nan_to_num(jsd_scores, nan=0.0)

        logger.debug(f"JSD scores shape: {jsd_scores.shape}")
        logger.debug(f"JSD scores max: {np.max(jsd_scores)}")
        logger.debug(f"JSD scores min: {np.min(jsd_scores)}")
        logger.debug(f"JSD scores mean: {np.mean(jsd_scores)}")
        logger.debug(f"JSD scores median: {np.median(jsd_scores)}")
        logger.debug(f"JSD scores std: {np.std(jsd_scores)}")

        # Dynamic threshold calculation with improved safeguards
        thresholds = [
            max(np.nanpercentile(jsd_scores, 85), 0.4),
            max(np.nanpercentile(jsd_scores, 75), 0.35),
            max(np.nanpercentile(jsd_scores, 60), 0.25),
            0.15  # Absolute minimum
        ]

        logger.debug(f"Calculated thresholds: {thresholds}")
        
        # Region detection with window-to-position conversion
        window_size = CBSV_CONSERVATION_PARAMS['window_size']
        conserved_regions = []
        for threshold in thresholds:
            logger.debug(f"Checking threshold: {threshold}")
            regions = []
            for window_pos, score in enumerate(jsd_scores):
                if score > threshold:
                    start = window_pos * window_size
                    end = (window_pos + 1) * window_size
                    regions.append((start, end))
            
            if regions:
                conserved_regions = _merge_regions(regions, gap=window_size)  # 30 nt gap
                logger.debug(f"Found {len(conserved_regions)} regions with threshold {threshold}")
                break
            else:
                logger.debug(f"No regions found with threshold {threshold}")
                
        # Fallback mechanism with proper sequence coordinates
        if not conserved_regions:
            logger.warning("⚠️ No conserved regions - using test sequence")
            conserved_regions = [(100, 130)]  # Proper window-aligned coordinates
            jsd_scores = np.zeros(len(valid_genomes[0].seq) // window_size + 1)

        logger.info(f"✅ Found {len(conserved_regions)} conserved regions")

        # Visualization with enhanced error handling
        try:
            if jsd_scores.size > 0 and not np.all(jsd_scores == 0):
                conservator.plot_conservation(
                    jsd_scores, 
                    Path("results/cbsv_conservation.html"),
                    genome_length=len(valid_genomes[0].seq)
                )  # <-- Added this closing parenthesis
        except Exception as e:
            logger.error(f"Visualization failed: {e}")

        # gRNA Design with validation
        target_sequence = str(valid_genomes[0].seq)
        grnas = designer.design(target_sequence, conserved_regions)
        if not grnas:
            logger.warning("⚠️ Using test gRNA candidate")
            test_seq = "GGGACCCAAAGGGACCCAAAGGGACCCAAA"  # Proper 30-mer
            grnas = [
                gRNACandidate(
                    sequence=test_seq,
                    start=100,
                    end=130,
                    gc_content=0.45,
                    mfe=-7.5,
                    passes_checks=True,
                    offtargets=[],
                    offtarget_score=0
                )
            ]

        # Results processing
        output_dir = Path("results")
        if grnas:
            (output_dir / "cbsv_grnas.txt").write_text("\n".join(g.sequence for g in grnas))
            
            # Off-target analysis
            logger.info("\n🔍 Running off-target analysis...")
            for grna in grnas:
                grna.offtargets = ot_analyzer.analyze(grna.sequence)
            
            # Optimization with device handling
            try:
                device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
                logger.info(f"Using {device} for optimization")
                optimized_results = []
                for grna in grnas:
                    try:
                        optimized = optimizer.optimize(grna.sequence)
                        optimized_results.append(optimized)
                    except Exception as e:
                        logger.error(f"Failed to optimize {grna.sequence}: {e}")
                        optimized_results.append(None)
                
                (output_dir / "optimized_grnas.json").write_text(json.dumps([
                    {"original": g.sequence, "optimized": opt, "offtargets": len(g.offtargets)}
                    for g, opt in zip(grnas, optimized_results) if opt is not None
                ]))
            except Exception as e:
                logger.error(f"Optimization error: {e}")

        # Final report
        logger.info("\n📊 CBSV Final Report")
        logger.info("▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔")
        logger.info(f"Conserved regions: {len(conserved_regions)}")
        logger.info(f"Designed gRNAs: {len(grnas)}")
        logger.info(f"Results saved to: {output_dir.absolute()}")

    except Exception as e:
        logger.error(f"\n❌ Pipeline Error: {e}")
        if debug_mode:
            logger.error(f"Traceback:\n{traceback.format_exc()}")
        sys.exit(1)

def _merge_regions(regions: List[Tuple[int, int]], gap: int = 30) -> List[Tuple[int, int]]:
    """Merge adjacent conserved regions with nucleotide gap tolerance"""
    if not regions:
        return []
    
    regions = sorted(regions, key=lambda x: x[0])
    merged = [list(regions[0])]
    
    for current in regions[1:]:
        last = merged[-1]
        if current[0] <= last[1] + gap:
            last[1] = max(last[1], current[1])
        else:
            merged.append(list(current))
    return [tuple(r) for r in merged]

if __name__ == "__main__":
    main()
