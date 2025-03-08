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

LOCAL_GENOMES = [
    SeqRecord(
        Seq("ATGCGATAGCATCGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA" * 100),
        id="Local_HA_1",
        description="Hemagglutinin (Test Strain A)"
    ),
    SeqRecord(
        Seq("ATGCGATAGCATCGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA".replace("C", "T", 10) * 100),
        id="Local_HA_2",
        description="Hemagglutinin (Test Strain B)"
    )
]

# ... (keep existing imports and LOCAL_GENOMES) ...

def main():
    try:
        # Clear previous runs
        shutil.rmtree("alignments", ignore_errors=True)
        shutil.rmtree("results", ignore_errors=True)

        # Initialize components
        fetcher = GenomeFetcher(api_key=os.getenv("NCBI_API_KEY_2025"))
        aligner = AlignmentEngine(max_threads=8)
        conservator = ConservationAnalyzer(window_size=30)
        designer = Cas13gRNADesigner()

        # Genome acquisition
        genomes = []
        try:
            print("🕵️ Attempting NCBI Datasets API v2 fetch...")
                genomes = fetcher.fetch_genomes(
                    target="11320",  # Influenza A virus taxid
                    gene_name="HA",  # Hemagglutinin
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
        thresholds = [0.8, 0.7, 0.6, max(0.3, max_jsd * 0.9)]  # More conservative
        
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
            print(f"{i}. {grna['spacer']} (pos {grna['start']}-{grna['end']})")

        # Save outputs
        output_dir = Path("results")
        (output_dir / "grnas.txt").write_text("\n".join(g['spacer'] for g in grnas))
        print(f"\n📁 Results saved to {output_dir}/")

    except Exception as e:
        print(f"\n❌ Pipeline Error: {e}")
        print("💡 Final Verification Steps:")
        print("1. Check RNAfold output manually:")
        print(f"   echo '>TEST\nATGCGATAGCAT' | RNAfold --noPS")
        print("2. Validate alignment file:")
        print("   muscle -in alignments/INPUT.fasta -out test.fasta")
        print("3. Test conservation analysis:")
        print("   python3 -m guidex.conservation alignments/INPUT.fasta")
        sys.exit(1)

if __name__ == "__main__":
    main()
