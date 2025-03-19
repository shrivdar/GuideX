from Bio import AlignIO
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
from guidex.genome_fetcher import GenomeFetcher
from conservation import ConservationAnalyzer
from alignment_engine import AlignmentEngine
from guidex.core import Cas13gRNADesigner
from guidex.grna.off_target import OffTargetAnalyzer
from guidex.grna.optimizer import Cas13Optimizer
import torch
import logging
import json
from datetime import datetime
import pandas as pd
import numpy as np
import importlib
import conservation
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)  # Add this line

logger.debug(f"Current working directory: {os.getcwd()}")
logger.debug(f"Python path: {sys.path}")

# CBSV-optimized parameters
CBSV_THRESHOLDS = [0.65, 0.55, 0.45, 0.35]
CBSV_CONSERVATION_PARAMS = {
    'window_size': 30,
    'max_gap': 0.7,
    'pseudocount': 10.0
}

def main():
    try:
        # Debug and logging initialization
        debug_mode = os.getenv('GUIDEX_DEBUG', 'false').lower() == 'true'
        log_level = logging.DEBUG if debug_mode else logging.INFO
        
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' if debug_mode else '%(message)s',
            handlers=[
                logging.FileHandler("guidex.log"),
                logging.StreamHandler()
            ]
        )
        logger = logging.getLogger(__name__)
        logger.info("🚀 Starting CBSV-Optimized GuideX Pipeline")
        
        if debug_mode:
            logger.debug("🐛 DEBUG MODE ACTIVATED")
            logger.debug(f"System version: {sys.version}")
            logger.debug(f"NumPy: {np.__version__}, PyTorch: {torch.__version__}")

        # Directory cleanup
        for dir_path in ["alignments", "results"]:
            if Path(dir_path).exists():
                shutil.rmtree(dir_path)
                logger.debug(f"Cleaned directory: {dir_path}")
        Path("alignments").mkdir(exist_ok=True)
        Path("results").mkdir(exist_ok=True)

        # Component initialization
        logger.debug("Initializing CBSV-optimized components...")
        fetcher = GenomeFetcher(api_key=os.getenv("NCBI_API_KEY_2025"))
        aligner = AlignmentEngine(max_threads=8)
        logger.debug(f"Conservation module path: {conservation.__file__}")
        conservator = conservation.ConservationAnalyzer(**CBSV_CONSERVATION_PARAMS)
        designer = Cas13gRNADesigner()
        designer.configure(
            gc_range=(0.35, 0.65),
            mfe_threshold=-5.0
        )
        ot_analyzer = OffTargetAnalyzer(reference_genome=Path("genomes/hg38.fa"))
        optimizer = Cas13Optimizer(designer=designer)

        try:
            importlib.reload(conservation)
            logger.debug("Successfully reloaded conservation module")
        except Exception as e:
            logger.error(f"Module reload failed: {e}")
        # Genome acquisition
        genomes = []
        try:
            logger.info("🕵️ Fetching CBSV genomes...")
            genomes = fetcher.fetch_genomes(
                target="Cassava brown streak virus",
                genome_type="genome",
                limit=12
            )
            
            if debug_mode and genomes:
                logger.debug("NCBI fetch results:")
                logger.debug(f"First genome ID: {genomes[0].id}")
                logger.debug(f"First 50 bases: {str(genomes[0].seq[:50])}")


            if not genomes:
                logger.warning("ℹ️ No genomes found - trying protein database")
                genomes = fetcher.fetch_genomes(
                    target="Cassava brown streak virus",
                    genome_type="protein",
                    limit=12
                )

        except Exception as e:
            logger.error(f"🔴 Genome fetch error: {e}")
            sys.exit(1)

        # Sequence validation
        valid_genomes = [g for g in genomes if 8000 <= len(g.seq) <= 10000]
        logger.debug(f"Valid genomes after filtering: {len(valid_genomes)}")
        
        if len(valid_genomes) < 2:
            raise ValueError(f"Only {len(valid_genomes)} valid CBSV genomes (need ≥2)")

        # Alignment process
        logger.info("\n🧬 Starting CBSV alignment...")
        aligned_file = None
        try:
            aligned_file = aligner.align(valid_genomes, Path("alignments/cbsv_alignment.fasta"))
            if not aligned_file.exists():
                raise FileNotFoundError(f"Alignment failed at {aligned_file}")
                
            logger.info(f"🔍 Alignment saved to: {aligned_file}")
            
            if debug_mode:
                try:
                    alignment = AlignIO.read(aligned_file, "fasta")
                    logger.debug(f"Alignment contains {len(alignment)} sequences")
                    logger.debug(f"Alignment length: {alignment.get_alignment_length()} bp")
                    
                    valid_seqs = [seq for seq in alignment if len(seq.seq) == alignment.get_alignment_length()]
                    logger.debug(f"Valid sequences: {len(valid_seqs)}/{len(alignment)}")
                    
                    logger.debug("First 3 headers:")
                    for rec in alignment[:3]:
                        logger.debug(f"- {rec.id[:50]}")

                except Exception as e:
                    logger.error(f"Alignment validation failed: {e}")

        except Exception as e:
            logger.error(f"Alignment failed: {e}")
            if debug_mode:
                logger.error(f"Traceback:\n{traceback.format_exc()}")
            sys.exit(1)

        # Conservation analysis
        logger.info("\n🔎 Identifying CBSV conserved regions...")
        jsd_scores = []
        conserved_regions = []
        
        try:
            conservation_result = conservator.calculate_jsd(aligned_file)
            
            if not (isinstance(conservation_result, tuple) and len(conservation_result) == 2):
                raise ValueError(f"Invalid conservation format: {type(conservation_result)}")
                
            jsd_scores, valid_windows = conservation_result
            jsd_scores = np.nan_to_num(jsd_scores, nan=0.0)
            
            # Dynamic threshold calculation
            thresholds = CBSV_THRESHOLDS.copy()
            try:
                if len(jsd_scores) > 10:
                    thresholds[0] = np.percentile(jsd_scores, 85)
                    thresholds[1] = np.percentile(jsd_scores, 75)
                    thresholds[2] = np.percentile(jsd_scores, 60)
            except Exception as e:
                logger.warning(f"Using default CBSV thresholds: {e}")

            # Region detection
            conserved_regions = []
            for threshold in thresholds:
                regions = [
                    (pos, min(pos+30, len(jsd_scores)-1))
                    for pos, score in enumerate(jsd_scores)
                    if score > threshold
                ]
                if regions:
                    conserved_regions = self._merge_regions(regions)
                    break
                    
            logger.info(f"✅ Found {len(conserved_regions)} CBSV conserved regions")

            if debug_mode:
                logger.debug("Conservation analysis raw results:")
                logger.debug(f"- JSD scores: {len(jsd_scores)}")
                logger.debug(f"- Valid windows: {valid_windows}")
                logger.debug(f"- NaN count: {np.isnan(jsd_scores).sum()}")

        except Exception as e:
            logger.error(f"Conservation failed: {e}")
            sys.exit(1)

        # Visualization and debugging
        try:
            if jsd_scores and len(jsd_scores) > 10:
                conservator.plot_conservation(jsd_scores, Path("results/cbsv_conservation.html"))
                logger.info("📊 CBSV conservation plot generated")
        except Exception as e:
            logger.error(f"Visualization failed: {e}")

        if debug_mode:
            try:
                debug_data = {
                    'metadata': {
                        'timestamp': datetime.now().isoformat(),
                        'version': '1.2'
                    },
                    'jsd_scores': jsd_scores,
                    'thresholds': thresholds,
                    'conserved_regions': conserved_regions
                }
                np.save("results/cbsv_debug.npy", debug_data)
                logger.debug("Debug data persisted")
            except Exception as e:
                logger.error(f"Debug data error: {e}")

        # gRNA Design
        grnas = []
        if conserved_regions:
            try:
                target_sequence = str(valid_genomes[0].seq)
                grnas = designer.design(target_sequence, conserved_regions)
                
                if not grnas:
                    logger.warning("⚠️ Relaxing gRNA constraints")
                    designer.config.mfe_threshold = -1.0
                    grnas = designer.design(target_sequence, conserved_regions)
                    
                logger.info(f"🧬 Designed {len(grnas)} gRNAs")

            except Exception as e:
                logger.error(f"gRNA design failed: {e}")

        # Results handling
        output_dir = Path("results")
        if grnas:
            (output_dir / "cbsv_grnas.txt").write_text("\n".join(g.sequence for g in grnas))
            
            # Off-target analysis
            logger.info("\n🔍 Running CBSV off-target analysis...")
            for grna in grnas:
                grna.offtargets = ot_analyzer.analyze(grna.sequence)
                if debug_mode and grna.offtargets:
                    logger.debug(f"gRNA {grna.sequence} has {len(grna.offtargets)} off-targets")

            # Optimization
            try:
                device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
                logger.info(f"⚙️ Optimizing on {device}")
                optimized_results = []
                for grna in grnas:
                    try:
                        optimized = optimizer.optimize(grna.sequence)
                        optimized_results.append({
                            "original": grna.sequence,
                            "optimized": optimized,
                            "offtargets": len(grna.offtargets)
                        })
                    except Exception as e:
                        logger.warning(f"Optimization failed for {grna.sequence[:8]}: {e}")
                (output_dir / "optimized_grnas.json").write_text(json.dumps(optimized_results))
            except Exception as e:
                logger.error(f"Optimization error: {e}")

        # Final report
        logger.info("\n📊 CBSV Final Report")
        logger.info("▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔")
        logger.info(f"Conserved regions: {len(conserved_regions)}")
        logger.info(f"Valid gRNAs: {len(grnas)}")
        if grnas:
            logger.info(f"gRNAs with off-targets: {sum(1 for g in grnas if g.offtargets)}")
            logger.info(f"Optimized: {len(optimized_results)}")
        
        if debug_mode:
            import psutil
            logger.debug("\nSystem Resources:")
            logger.debug(f"Memory: {psutil.virtual_memory().percent}%")
            logger.debug(f"CPU: {psutil.cpu_percent()}%")

        logger.info(f"\n📁 Results saved to: {output_dir.absolute()}")

    except Exception as e:
        logger.error(f"\n❌ Pipeline Error: {e}")
        if debug_mode:
            logger.error(f"Traceback:\n{traceback.format_exc()}")
        logger.info("\n🔧 Debug Tips:")
        logger.info("1. Check alignment: muscle -in alignments/cbsv_alignment.fasta")
        logger.info("2. Inspect conservation scores: results/cbsv_debug.npy")
        sys.exit(1)

def _merge_regions(regions, gap=10):
    """Merge adjacent conserved regions"""
    if not regions:
        return []
        
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
