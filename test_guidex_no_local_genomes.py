
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import BiopythonDeprecationWarning
import warnings
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
from pathlib import Path
import sys
import os
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
import pandas as pd
import numpy as np
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)  # Add this line

def main():
    try:
        # Add version debug
        logger.debug(f"Python version: {sys.version}")
        logger.debug(f"PyTorch version: {torch.__version__}")
        logger.debug(f"NumPy version: {np.__version__}")

        # Initialize with debug logging
        debug_mode = os.getenv('GUIDEX_DEBUG', 'false').lower() == 'true'
        if debug_mode:
            logging.basicConfig(
                level=logging.DEBUG,
                format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                handlers=[
                    logging.FileHandler("guidex.debug.log"),
                    logging.StreamHandler()
                ]
            )
            logger.info("🐛 DEBUG MODE ACTIVATED")

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

        # Alignment process with file verification
        logger.info("\n🧬 Starting alignment...")
        aligned_file = aligner.align(valid_genomes, Path("alignments"))
        logger.info(f"🔍 Alignment saved to: {aligned_file}")

        if debug_mode:
            logger.debug("Alignment file inspection:")
            try:
                with open(aligned_file) as f:
                    lines = f.readlines()
                    logger.debug(f"File length: {len(lines)} lines")
                    logger.debug("First 3 sequences:")
                    for line in lines[:6]:
                        if line.startswith(">"):
                            logger.debug(f"Header: {line.strip()}")
            except Exception as e:
                logger.error(f"Alignment file read error: {e}")

        # Conservation analysis with enhanced debugging
        logger.info("\n🔎 Identifying conserved regions...")
        try:
            jsd_scores, valid_regions = conservator.calculate_jsd(aligned_file)
            logger.info(f"Analyzed {valid_regions} valid genomic regions")

            if debug_mode:
                logger.debug("Conservation analysis details:")
                logger.debug(f"JSD scores count: {len(jsd_scores)}")
                logger.debug(f"NaN values: {sum(np.isnan(s) for s in jsd_scores)}")
                logger.debug(f"Score range: {np.nanmin(jsd_scores):.3f}-{np.nanmax(jsd_scores):.3f}")
                
                # Plot conservation scores
                plt.figure(figsize=(10, 4))
                plt.plot(jsd_scores)
                plt.title("Raw Conservation Scores")
                plt.savefig("results/conservation_raw.png")
                plt.close()
                logger.debug("Saved raw conservation plot")

            # Threshold adaptation debugging
            clean_scores = [s for s in jsd_scores if not np.isnan(s)]
            thresholds = [
                (max(clean_scores) * 0.9 if clean_scores else 0,
                (max(clean_scores) * 0.8 if clean_scores else 0,
                max((min(clean_scores) + 0.2 if clean_scores else 0), 0.6),
                0.5
            ]

            conserved_regions = []
            for i, threshold in enumerate(thresholds):
                current_regions = [(pos, pos+30) for pos, s in enumerate(jsd_scores)
                                  if not np.isnan(s) and s > threshold]
                if debug_mode:
                    logger.debug(f"Threshold {i+1} ({threshold:.3f}): {len(current_regions)} regions")
                    if current_regions:
                        logger.debug(f" Example region: {current_regions[0][0]}-{current_regions[0][1]}")
                        logger.debug(f" Region scores: {jsd_scores[current_regions[0][0]:.3f}-{jsd_scores[current_regions[0][1]:.3f}")

                if current_regions:
                    conserved_regions = current_regions
                    break

            logger.info(f"✅ Found {len(conserved_regions)} conserved regions")

        except Exception as e:
            logger.error(f"Conservation analysis failed: {str(e)}")
            if debug_mode:
                logger.error(f"Error details:\n{traceback.format_exc()}")
                logger.debug("Alignment file metadata:")
                logger.debug(f"File size: {aligned_file.stat().st_size} bytes")
                logger.debug(f"Last modified: {aligned_file.stat().st_mtime}")
            conserved_regions = []

        # Visualization with error protection
        try:
            conservator.plot_conservation(jsd_scores, Path("results/conservation.html"))
            logger.info("📊 Conservation plot generated")
        except Exception as e:
            logger.error(f"Visualization failed: {str(e)}")
            if debug_mode:
                logger.error(f"Plot error details:\n{traceback.format_exc()}")

        # Debug data exports
        if debug_mode:
            logger.debug("Saving debug data...")
            np.save("results/js_scores.npy", jsd_scores)
            pd.DataFrame({
                'position': range(len(jsd_scores)),
                'jsd': jsd_scores
            }).to_csv("results/js_scores.csv", index=False)
            
            with open("results/thresholds.json", 'w') as f:
                json.dump({
                    'calculated': [float(t) for t in thresholds],
                    'final_threshold': float(threshold) if conserved_regions else None,
                    'regions_count': len(conserved_regions)
                }, f)

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
