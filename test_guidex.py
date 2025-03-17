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

LOCAL_GENOMES = [
    SeqRecord(
        Seq(("ATGCGATAGCATCGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA"  # Known off-target rich
             "GGTACCGAGCTCGTTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA") * 50),
        id="Test_OffTarget_Rich",
        description="Contains known off-target sequences"
    ),
    SeqRecord(
        Seq(("ATGCGATAGCATCGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA"
             "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC") * 50),
        id="Test_Control",
        description="Low complexity control"
    )
]

# ... (keep existing imports and LOCAL_GENOMES) ...

def main():
    try:
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
        ot_analyzer = OffTargetAnalyzer(Path("genomes/hg38"))
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
        jsd_scores = conservator.calculate_jsd(aligned_file)
        
        # Dynamic threshold adjustment
        max_jsd = max(jsd_scores) if jsd_scores else 0
        conserved_regions = []
        thresholds = [
            0.8, 
            0.7, 
            0.6, 
            max(0.15, max_jsd * 0.9),  # Ensures minimum threshold
            max_jsd - 0.05             # Catch borderline cases
        ]  # More conservative
        
        for threshold in thresholds:
            conserved_regions = [(i, i+30) for i, score in enumerate(jsd_scores) if score > threshold]
            if conserved_regions:
                print(f"✅ Found {len(conserved_regions)} regions at JSD > {threshold:.2f}")
                break
        
        if not conserved_regions:
            print(f"⚠️ No conserved regions found (max JSD: {max_jsd:.2f})")
            print("💡 Try providing more similar sequences or adjust conservation thresholds")
            conserved_regions = []
                    
        print(f"✅ Found {len(conserved_regions)} conserved regions")

        # Visualization with directory check
        Path("results").mkdir(parents=True, exist_ok=True)
        conservator.plot_conservation(jsd_scores, Path("results/conservation.html"))
        print("📊 Conservation plot generated")

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
