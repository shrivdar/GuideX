
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
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)  # Add this line

# ... (keep existing imports and LOCAL_GENOMES) ...

def main():
    try:
        jsd_scores = []
        valid_regions = 0
        # Initialize with debug logging
        logging.basicConfig(level=logging.INFO)
        logger.info("🚀 Starting GuideX pipeline")
        logging.basicConfig(
            level=logging.INFO,
            format='%(message)s',  # Clean format without module names
            handlers=[
                logging.FileHandler("guidex.log"),  # Save to file
                logging.StreamHandler()  # Console output
            ]
        )
        logger.info("🚀 Starting GuideX pipeline")
        # Clean previous runs
        shutil.rmtree("alignments", ignore_errors=True)
        shutil.rmtree("results", ignore_errors=True)

        # Component initialization
        fetcher = GenomeFetcher(api_key=os.getenv("NCBI_API_KEY_2025"))
        aligner = AlignmentEngine(max_threads=8)
        conservator = ConservationAnalyzer(window_size=30)
        designer = Cas13gRNADesigner()
        ot_analyzer = OffTargetAnalyzer(
            reference_genome=Path("genomes/hg38.fa"),  # Direct genome file
            max_mismatches=3,
            output_dir=Path("results/off_targets"),
            window_size=28
        )
        optimizer = Cas13Optimizer(designer)
        
        # Load model weights if available
        if Path("weights/optimizer.pth").exists():
            optimizer.load_state_dict(torch.load("weights/optimizer.pth"))

        # Genome acquisition
        genomes = []
        try:
            logger.info("🕵️ Attempting NCBI fetch...")
            genomes = fetcher.fetch_genomes(
                target="Influenza A virus",
                gene="HA",
                genome_type="gene",
                limit=10
            )
            if not genomes:
                print("ℹ️ No genomes found - trying protein database")
                genomes = fetcher.fetch_genomes(
                    target="Influenza A virus HA",
                    genome_type="protein",
                    limit=10
                )
            if not genomes:
                raise RuntimeError("NCBI returned 0 genomes")
            print(f"✅ Retrieved {len(genomes)} NCBI genomes")
        except Exception as e:
            print(f"🔴 NCBI Error: {e}")
            print("🔄 Falling back to local genomes")
            genomes = LOCAL_GENOMES

        # Sequence validation
        valid_genomes = [g for g in genomes if 5000 <= len(g.seq) <= 15000]
        if len(valid_genomes) < 2:
            raise ValueError(f"Only {len(valid_genomes)} valid genomes (need ≥2)")
            
        # Alignment
        print("\n🧬 Starting alignment...")
        aligned_file = aligner.align(valid_genomes, Path("alignments"))
        print(f"🔍 Alignment saved to: {aligned_file}")

        # Conservation analysis with adaptive thresholds
        print("\n🔎 Identifying conserved regions...")
        debug_mode = os.getenv('GUIDEX_DEBUG', 'false').lower() == 'true'  # New debug flag
        
        try:
            # Core analysis remains the same
            jsd_scores, valid_regions = conservator.calculate_jsd(aligned_file)
            logger.info(f"Analyzed {valid_regions} valid genomic regions")
            
            if debug_mode:
                logger.debug(f"Raw JSD data count: {len(jsd_scores)}")
                logger.debug(f"Initial NaN count: {sum(np.isnan(s) for s in jsd_scores)}")
        
            if not jsd_scores:
                raise ValueError("All JSD scores were invalid - check alignment quality")
        
            # Enhanced NaN filtering with debug info
            clean_scores = [s for s in jsd_scores if not np.isnan(s)]
            if debug_mode:
                logger.debug(f"Filtered {len(jsd_scores)-len(clean_scores)} NaN values")
                logger.debug(f"Score range: {np.nanmin(jsd_scores):.3f}-{np.nanmax(jsd_scores):.3f}")
        
            # Threshold calculation with debug logging
            thresholds = [
                (max(clean_scores) * 0.9) if clean_scores else 0,
                (max(clean_scores) * 0.8) if clean_scores else 0,
                max((min(clean_scores) + 0.2 if clean_scores else 0), 0.6),
                0.5
            ]
            if debug_mode:
                logger.debug(f"Threshold sequence: {[f'{t:.3f}' for t in thresholds]}")
        
            # Enhanced threshold testing with debug output
            conserved_regions = []
            for i, threshold in enumerate(thresholds):
                current_regions = [(i, i+30) for i, s in enumerate(jsd_scores)
                                  if not np.isnan(s) and s > threshold]
                
                if debug_mode:
                    logger.debug(f"Threshold {i+1} ({threshold:.3f}): {len(current_regions)} candidates")
                
                if current_regions:
                    conserved_regions = current_regions
                    logger.info(f"Selected threshold {threshold:.3f} with {len(conserved_regions)} regions")
                    break
        
            if debug_mode and not conserved_regions:
                logger.debug("No conserved regions found across all thresholds")

        except Exception as e:
            logger.error(f"Conservation analysis failed: {str(e)}")
            if debug_mode:
                # Add alignment file inspection
                from Bio import AlignIO
                logger.debug("Alignment file contents:")
                try:
                    alignment = AlignIO.read(aligned_file, "fasta")
                    logger.debug(f"Sequences: {len(alignment)}, Length: {alignment.get_alignment_length()}")
                    logger.debug(f"First sequence: {str(alignment[0].seq)[:50]}...")
                except Exception as parse_error:
                    logger.debug(f"Failed to read alignment: {str(parse_error)}")
            
            conserved_regions = []
        
        print(f"✅ Found {len(conserved_regions)} conserved regions")
        
        # Enhanced visualization with debug features
        Path("results").mkdir(parents=True, exist_ok=True)
        try:
            conservator.plot_conservation(
                jsd_scores, 
                Path("results/conservation.html"),
            )
            print("📊 Conservation plot generated")
        except Exception as e:
            logger.error(f"Visualization failed: {str(e)}")
            if debug_mode:
                logger.error(f"Plot error details:\n{traceback.format_exc()}")
        
        # New debug output files
        if debug_mode:
            # Save numerical data for inspection
            np.save("results/js_scores.npy", jsd_scores)
            pd.DataFrame({
                'position': range(len(jsd_scores)),
                'jsd': jsd_scores
            }).to_csv("results/js_scores.csv", index=False)
            
            # Save threshold information
            with open("results/thresholds.json", 'w') as f:
                json.dump({
                    'calculated': thresholds,
                    'final': conserved_regions[0][1] if conserved_regions else None
                }, f)
            
            logger.debug("Debug data saved to results/ directory")

        if conserved_regions:
            print("\n🔬 Designing Cas13 gRNAs...")
            target_sequence = str(valid_genomes[0].seq)
            grnas = designer.design(target_sequence, conserved_regions)
            
            if not grnas:
                print("⚠️ Designed 0 gRNAs - relaxing constraints")
                designer.config.mfe_threshold = -1.0  # Relax MFE
                designer.config.gc_min = 0.3
                grnas = designer.design(target_sequence, conserved_regions)
        else:
            print("\n⏭️ Skipping gRNA design - no conserved regions")
            grnas = []
        
        # Output results
        print(f"\n🎉 Success! Designed {len(grnas)} gRNAs:")
        for i, grna in enumerate(grnas[:5], 1):
            print(f"{i}. {grna.sequence} (pos {grna.start}-{grna.end})")

        # Save outputs
        output_dir = Path("results")
        (output_dir / "grnas.txt").write_text("\n".join(g.sequence for g in grnas))  # ✅ Use .sequence attribute
        print(f"\n📁 Results saved to {output_dir}/")

        if grnas:
            logger.info("\n🔍 Running off-target analysis...")
            ot_results = []
            
            for grna in grnas:
                grna.offtargets = ot_analyzer.analyze(grna.sequence)
                grna.offtarget_score = len(grna.offtargets)
                ot_results.append({
                    "sequence": grna.sequence,
                    "offtarget_count": grna.offtarget_score
                })
        
            # Add json import verification
            if 'json' not in globals():
                raise ImportError("JSON module not loaded")
                
            # Now safe to use json.dumps
            analysis_summary = {
                "total_grnas": len(grnas),
                "grnas_with_offtargets": sum(1 for g in grnas if g.offtarget_score > 0),
                "details": ot_results
            }
            (output_dir / "offtarget_summary.json").write_text(json.dumps(analysis_summary))
            
            logger.info("\n⚙️ Optimizing gRNAs...")
            try:
                device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
                logger.info(f"Using compute device: {device}")
                optimizer = optimizer.to(device)
                
                optimized_results = []
                for grna in grnas:
                    try:
                        optimized = optimizer.optimize(grna.sequence)
                        optimized_results.append({
                            "original": grna.sequence,
                            "optimized": optimized,
                            "offtarget_score": grna.offtarget_score
                        })
                    except Exception as e:
                        logger.debug(f"Optimization failed for {grna.sequence[:8]}: {str(e)}")
                        
                # Save optimization results
                (output_dir / "optimized_grnas.json").write_text(json.dumps(optimized_results))
                
            except Exception as e:
                logger.error(f"Optimization error: {str(e)}")
                if "CUDA" in str(e):
                    logger.info("💡 Try running on CPU: export CUDA_VISIBLE_DEVICES=''")

        # Final report
        logger.info("\n📊 Final Report")
        logger.info("▔▔▔▔▔▔▔▔▔▔▔")
        logger.info(f"Total gRNAs designed: {len(grnas)}")
        if grnas:
            logger.info(f"gRNAs with off-targets: {analysis_summary['grnas_with_offtargets']}")
            logger.info(f"Optimized gRNAs: {len(optimized_results)}")
        logger.info(f"Results directory: {output_dir.absolute()}")

    except Exception as e:
        logger.error(f"\n❌ Pipeline Error: {e}")
        logger.info("💡 Debug Steps:")
        logger.info("1. Verify input sequences with: echo 'ATGCGATAGCAT' | RNAfold --noPS")
        logger.info("2. Check alignment file: muscle -in alignments/INPUT.fasta")
        logger.info("3. Test conservation: python3 -m guidex.conservation alignments/INPUT.fasta")
        sys.exit(1)

if __name__ == "__main__":
    main()
