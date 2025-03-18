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
import warnings
from guidex.genome_fetcher import GenomeFetcher
from guidex.conservation import ConservationAnalyzer
from alignment_engine import AlignmentEngine
from guidex.core import Cas13gRNADesigner
from guidex.grna.off_target import OffTargetAnalyzer
from guidex.grna.optimizer import Cas13Optimizer
import torch
import logging
import json
import traceback
from datetime import datetime
import pandas as pd
import numpy as np
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)  # Add this line

def main():
    try:
        # Initialize debug mode first
        debug_mode = os.getenv('GUIDEX_DEBUG', 'false').lower() == 'true'
        log_level = logging.DEBUG if debug_mode else logging.INFO
        
        # Configure logging with detailed debug format
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s' if debug_mode else '%(message)s',
            handlers=[
                logging.FileHandler("guidex.log"),
                logging.StreamHandler()
            ]
        )
        logger.info("🚀 Starting GuideX pipeline")
        
        if debug_mode:
            logger.debug("🐛 DEBUG MODE ACTIVATED")
            logger.debug(f"System version: {sys.version}")
            logger.debug(f"NumPy version: {np.__version__}")
            logger.debug(f"PyTorch version: {torch.__version__}")

        # Clean previous runs with verification
        for dir_path in ["alignments", "results"]:
            if Path(dir_path).exists():
                shutil.rmtree(dir_path)
                logger.debug(f"Cleaned directory: {dir_path}")
        Path("alignments").mkdir(exist_ok=True)
        Path("results").mkdir(exist_ok=True)

        # Component initialization with version checks
        logger.debug("Initializing pipeline components...")
        fetcher = GenomeFetcher(api_key=os.getenv("NCBI_API_KEY_2025"))
        aligner = AlignmentEngine(max_threads=8)
        conservator = ConservationAnalyzer(window_size=30)
        designer = Cas13gRNADesigner()
        
        logger.debug(f"Alignment engine: {aligner.__class__.__name__}")
        logger.debug(f"Conservation analyzer: {conservator.__class__.__name__}")

        # Genome acquisition with enhanced debugging
        genomes = []
        try:
            logger.info("🕵️ Attempting NCBI fetch...")
            genomes = fetcher.fetch_genomes(
                target="Influenza A virus",
                gene="HA",
                genome_type="gene",
                limit=10
            )
            
            if debug_mode:
                logger.debug("NCBI fetch results:")
                logger.debug(f"Retrieved {len(genomes)} genomes")
                if genomes:
                    logger.debug(f"First genome ID: {genomes[0].id}")
                    logger.debug(f"First 50 bases: {str(genomes[0].seq[:50])}")

            if not genomes:
                logger.warning("ℹ️ No genomes found - trying protein database")
                genomes = fetcher.fetch_genomes(
                    target="Influenza A virus HA",
                    genome_type="protein",
                    limit=10
                )

        except Exception as e:
            logger.error(f"🔴 NCBI Error: {e}")
            logger.info("🔄 Falling back to local genomes")
            genomes = LOCAL_GENOMES
            if debug_mode:
                logger.debug("Local genomes used:")
                for g in genomes[:3]:
                    logger.debug(f" - {g.id}: {len(g.seq)}bp")

        # Sequence validation with quality checks
        valid_genomes = [g for g in genomes if 5000 <= len(g.seq) <= 15000]
        logger.debug(f"Valid genomes after filtering: {len(valid_genomes)}")
        
        if len(valid_genomes) < 2:
            raise ValueError(f"Only {len(valid_genomes)} valid genomes (need ≥2)")

        # Alignment process with enhanced verification and error containment
        logger.info("\n🧬 Starting alignment...")
        aligned_file = None
        try:
            # Perform alignment with validation
            aligned_file = aligner.align(valid_genomes, Path("alignments"))
            if not aligned_file.exists():
                raise FileNotFoundError(f"Alignment failed - no output file at {aligned_file}")
                
            logger.info(f"🔍 Alignment saved to: {aligned_file}")
            
            if debug_mode:
                logger.debug("Alignment file structure verification:")
                try:
                    from Bio import AlignIO  # Ensure import is available
                    alignment = AlignIO.read(aligned_file, "fasta")
                    logger.debug(f"Alignment contains {len(alignment)} sequences")
                    logger.debug(f"Alignment length: {alignment.get_alignment_length()} bp")
                    
                    # Sequence quality checks
                    valid_seqs = [seq for seq in alignment if len(seq.seq) == alignment.get_alignment_length()]
                    logger.debug(f"Valid sequences: {len(valid_seqs)}/{len(alignment)}")
                    
                    # Sample output with sanitization
                    logger.debug("First 3 sequence headers:")
                    for rec in alignment[:3]:
                        logger.debug(f"- {rec.id[:50]}")  # Truncate long IDs

                except Exception as e:
                    logger.error(f"Alignment validation failed: {str(e)}")
                    if debug_mode:
                        logger.debug("Alignment file content sample:")
                        try:
                            with open(aligned_file) as f:
                                logger.debug(f"First 2 lines: {f.readline()}{f.readline()}")
                        except Exception as fe:
                            logger.debug(f"File read error: {str(fe)}")

        except Exception as e:
            logger.error(f"Alignment failed: {str(e)}")
            if debug_mode:
                logger.error(f"Alignment error traceback:\n{traceback.format_exc()}")
            raise RuntimeError("Critical alignment failure") from e

        # Conservation analysis with robust error handling and numerical stability
        logger.info("\n🔎 Identifying conserved regions...")
        jsd_scores = []
        valid_regions = 0
        conserved_regions = []
        thresholds = [0.6, 0.5, 0.6, 0.5]  # Safe defaults

        try:
            # Validate and unpack conservation results with type checking
            conservation_result = conservator.calculate_jsd(aligned_file)
            
            if not (isinstance(conservation_result, tuple) and len(conservation_result) == 2):
                raise ValueError(
                    f"Invalid conservation result format: {type(conservation_result)}. "
                    f"Expected tuple of (scores, valid_regions)"
                )
            
            jsd_scores, valid_regions = conservation_result
            
            if debug_mode:
                logger.debug("Conservation analysis raw results:")
                logger.debug(f"- Total positions: {len(jsd_scores)}")
                logger.debug(f"- Reported valid regions: {valid_regions}")
                logger.debug(f"- Initial NaN count: {np.isnan(jsd_scores).sum()}")

            # Data quality checks
            if not jsd_scores:
                raise ValueError("Empty conservation scores array")
                
            nan_ratio = np.isnan(jsd_scores).mean()
            if nan_ratio == 1:
                raise ValueError("All conservation scores are NaN")
            elif nan_ratio > 0.5:
                logger.warning(f"High NaN ratio: {nan_ratio:.1%}")

            # Threshold calculation with safe defaults
            clean_scores = [s for s in jsd_scores if not np.isnan(s)]
            try:
                if clean_scores:
                    threshold1 = np.nanpercentile(jsd_scores, 90)
                    threshold2 = np.nanpercentile(jsd_scores, 80)
                    threshold3 = max(np.nanpercentile(jsd_scores, 10) + 0.2, 0.6)
                else:
                    threshold1 = threshold2 = threshold3 = 0.6
                
                thresholds = [threshold1, threshold2, threshold3, 0.5]
                
                if debug_mode:
                    logger.debug("Computed Thresholds:")
                    logger.debug(f"1. {threshold1:.3f} (90th percentile)")
                    logger.debug(f"2. {threshold2:.3f} (80th percentile)")
                    logger.debug(f"3. {threshold3:.3f} (10th pct + 0.2)")
                    logger.debug(f"4. 0.500 (fixed fallback)")

            except Exception as e:
                logger.error(f"Threshold calculation error: {str(e)}")
                thresholds = [0.6, 0.5, 0.6, 0.5]
                if debug_mode:
                    logger.debug("Using safe default thresholds")

            # Region detection with bounds checking and validation
            conserved_regions = []
            for i, threshold in enumerate(thresholds):
                try:
                    current_regions = [
                        (pos, min(pos + 30, len(jsd_scores) - 1))  # Fixed bounds
                        for pos, s in enumerate(jsd_scores)
                        if not np.isnan(s) and s > threshold
                    ]
                    
                    if current_regions:
                        conserved_regions = current_regions
                        if debug_mode:
                            logger.debug(f"Threshold {i+1} ({threshold:.3f}) success:")
                            logger.debug(f"- Regions found: {len(current_regions)}")
                            if current_regions:
                                start, end = current_regions[0]
                                logger.debug(f"- First region: {start}-{end}")
                                logger.debug(f"- Region scores: {jsd_scores[start]:.3f}-{jsd_scores[end]:.3f}")
                        break
                        
                except Exception as e:
                    logger.warning(f"Region detection failed at threshold {i}: {str(e)}")
                    if debug_mode:
                        logger.debug(f"Failure context:")
                        logger.debug(f"- Threshold: {threshold}")
                        logger.debug(f"- JSD scores length: {len(jsd_scores)}")
                        logger.debug(f"- NaN count: {np.isnan(jsd_scores).sum()}")

            logger.info(f"✅ Found {len(conserved_regions)} conserved regions")

        except Exception as e:
            logger.error(f"Conservation analysis failed: {str(e)}")
            if debug_mode:
                logger.error(f"Error traceback:\n{traceback.format_exc()}")
                logger.debug("Alignment file diagnostics:")
                logger.debug(f"- Path: {aligned_file}")
                logger.debug(f"- Size: {aligned_file.stat().st_size} bytes" if aligned_file.exists() else "- File missing")
            conserved_regions = []

        # Visualization with enhanced safeguards and error reporting
        try:
            if jsd_scores and len(jsd_scores) > 10 and not np.all(np.isnan(jsd_scores)):
                conservator.plot_conservation(
                    jsd_scores, 
                    Path("results/conservation.html")
                )
                logger.info("📊 Conservation plot generated")
            else:
                logger.warning("Skipping visualization - insufficient valid data")
                if debug_mode:
                    valid_count = len([s for s in jsd_scores if not np.isnan(s)])
                    logger.debug(f"Valid data: {valid_count}/{len(jsd_scores)}")
                    if valid_count > 0:
                        logger.debug(f"Score range: {np.nanmin(jsd_scores):.3f}-{np.nanmax(jsd_scores):.3f}")

        except Exception as e:
            logger.error(f"Visualization failed: {str(e)}")
            if debug_mode:
                logger.error(f"Failed JSD scores sample:")
                logger.error(f"{jsd_scores[:10]}")

        # Debug data persistence with validation and error containment
        if debug_mode:
            try:
                logger.debug("Persisting debug datasets...")
                
                # Save numerical data with versioning
                debug_data = {
                    'metadata': {
                        'timestamp': datetime.now().isoformat(),
                        'version': '1.2'
                    },
                    'jsd_scores': jsd_scores,
                    'thresholds': thresholds,
                    'conserved_regions': conserved_regions
                }
                
                np.save("results/debug_data.npy", debug_data)
                
                # CSV with data validation
                debug_df = pd.DataFrame({
                    'position': range(len(jsd_scores)),
                    'jsd_score': jsd_scores,
                    'is_nan': np.isnan(jsd_scores)
                })
                debug_df.to_csv("results/conservation_debug.csv", index=False)
                
                # Threshold metadata with validation
                threshold_metadata = {
                    'calculated_thresholds': [float(t) for t in thresholds],
                    'final_threshold_used': float(thresholds[0]) if conserved_regions else None,
                    'conserved_region_count': len(conserved_regions),
                    'nan_ratio': np.isnan(jsd_scores).mean() if len(jsd_scores) > 0 else 1.0
                }
                with open("results/threshold_metadata.json", 'w') as f:
                    json.dump(threshold_metadata, f, indent=2)
                
                logger.debug("Debug data persisted successfully")
                
            except Exception as e:
                logger.error(f"Debug data persistence failed: {str(e)}")
                if debug_mode:
                    logger.error(f"Data dump error details:\n{traceback.format_exc()}")

        # gRNA Design with enhanced feedback
        grnas = []
        if conserved_regions:
            logger.info("\n🔬 Designing Cas13 gRNAs...")
            target_sequence = str(valid_genomes[0].seq)
            
            if debug_mode:
                logger.debug(f"Target sequence length: {len(target_sequence)}")
                logger.debug(f"First conserved region: {conserved_regions[0][0]}-{conserved_regions[0][1]}")
                logger.debug(f"Region sequence: {target_sequence[conserved_regions[0][0]:conserved_regions[0][1]]}")

            grnas = designer.design(target_sequence, conserved_regions)
            
            if not grnas:
                logger.warning("⚠️ Designed 0 gRNAs - relaxing constraints")
                designer.config.mfe_threshold = -1.0
                designer.config.gc_min = 0.3
                grnas = designer.design(target_sequence, conserved_regions)

        # Results handling with validation
        output_dir = Path("results")
        if grnas:
            logger.info(f"\n🎉 Success! Designed {len(grnas)} gRNAs:")
            for i, grna in enumerate(grnas[:5], 1):
                logger.info(f"{i}. {grna.sequence} (pos {grna.start}-{grna.end})")

            (output_dir / "grnas.txt").write_text("\n".join(g.sequence for g in grnas))
            
            # Off-target analysis
            logger.info("\n🔍 Running off-target analysis...")
            ot_results = []
            for grna in grnas:
                grna.offtargets = ot_analyzer.analyze(grna.sequence)
                grna.offtarget_score = len(grna.offtargets)
                ot_results.append({
                    "sequence": grna.sequence,
                    "offtarget_count": grna.offtarget_score
                })
                
                if debug_mode and grna.offtarget_score > 0:
                    logger.debug(f"gRNA {grna.sequence} has {grna.offtarget_score} off-targets")
                    logger.debug(f"First off-target: {grna.offtargets[0]}")

            # Optimization
            logger.info("\n⚙️ Optimizing gRNAs...")
            optimized_results = []
            try:
                device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
                logger.info(f"Using compute device: {device}")
                optimizer = optimizer.to(device)
                
                for grna in grnas:
                    try:
                        optimized = optimizer.optimize(grna.sequence)
                        optimized_results.append({
                            "original": grna.sequence,
                            "optimized": optimized,
                            "offtarget_score": grna.offtarget_score
                        })
                        if debug_mode:
                            logger.debug(f"Optimized {grna.sequence[:8]}... -> {optimized[:8]}...")
                    except Exception as e:
                        logger.warning(f"Optimization failed for {grna.sequence[:8]}: {str(e)}")

                (output_dir / "optimized_grnas.json").write_text(json.dumps(optimized_results))
                
            except Exception as e:
                logger.error(f"Optimization error: {str(e)}")
                if "CUDA" in str(e):
                    logger.info("💡 Try running on CPU: export CUDA_VISIBLE_DEVICES=''")

        # Final report with system stats
        logger.info("\n📊 Final Report")
        logger.info("▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔")
        logger.info(f"Conserved regions identified: {len(conserved_regions)}")
        logger.info(f"gRNAs designed: {len(grnas)}")
        if grnas:
            logger.info(f"gRNAs with off-targets: {sum(1 for g in grnas if g.offtarget_score > 0)}")
            logger.info(f"Successfully optimized: {len(optimized_results)}")
        
        if debug_mode:
            import psutil
            logger.debug("\nSystem Resources:")
            logger.debug(f"Memory usage: {psutil.virtual_memory().percent}%")
            logger.debug(f"CPU usage: {psutil.cpu_percent()}%")
            logger.debug(f"Disk usage: {psutil.disk_usage('/').percent}%")

        logger.info(f"\n📁 Results saved to: {output_dir.absolute()}")

    except Exception as e:
        logger.error(f"\n❌ Pipeline Error: {e}")
        if debug_mode:
            logger.error(f"Full error traceback:\n{traceback.format_exc()}")
        
        logger.info("\n🔧 Debugging Tips:")
        logger.info("1. Check alignment file: muscle -in alignments/INPUT.fasta")
        logger.info("2. Verify conservation scores: results/js_scores.csv")
        logger.info("3. Inspect thresholds: results/thresholds.json")
        
        sys.exit(1)

if __name__ == "__main__":
    main()
