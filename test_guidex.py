from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import BiopythonDeprecationWarning
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
from guidex.core import GuideXGrnaDesigner

LOCAL_GENOMES = [
    SeqRecord(
        Seq("ATGCGATAGCATCGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA" * 100),
        id="Local_HA_1",
        description="Hemagglutinin (Local Backup) | 5600bp"
    ),
    SeqRecord(
        Seq("ATGCGATAGCATGGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA" * 100),
        id="Local_NA_2",
        description="Neuraminidase (Local Backup) | 5600bp"
    )
]

def main():
    try:
        # Clear previous runs
        shutil.rmtree("alignments", ignore_errors=True)
        shutil.rmtree("results", ignore_errors=True)

        # Initialize components
        fetcher = GenomeFetcher(api_key=os.getenv("NCBI_API_KEY_2025"))
        aligner = AlignmentEngine(max_threads=8)
        conservator = ConservationAnalyzer(window_size=30)
        designer = GuideXGrnaDesigner(subtype="LwaCas13a")

        # Genome acquisition
        genomes = []
        try:
            print("🕵️ Attempting NCBI Datasets API v2 fetch...")
            genomes = fetcher.fetch_genomes(
                target="Influenza A virus (H1N1) hemagglutinin",  # Specific gene + subtype
                genome_type="protein",  # Switch to protein accessions
                limit=5,
                exclude_partial=True  # Add this if present in your GenomeFetcher class
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
        print(f"🔍 Alignment saved to: {aligned_file}")  # Now works

        # Conservation analysis
        print("\n🔎 Identifying conserved regions...")
        jsd_scores = conservator.calculate_jsd(aligned_file)
        conserved_regions = [(i, i+30) for i, score in enumerate(jsd_scores) if score > 0.8]
        print(f"✅ Found {len(conserved_regions)} conserved regions")

        # Visualization
        conservator.plot_conservation(jsd_scores, Path("results/conservation.html"))

        # gRNA Design
        print("\n🔬 Designing Cas13 gRNAs...")
        grnas = designer.design(str(valid_genomes[0].seq), conserved_regions)
        
        # Output results
        print(f"\n🎉 Success! Designed {len(grnas)} gRNAs:")
        for i, grna in enumerate(grnas[:5], 1):
            print(f"{i}. {grna['spacer']} (pos {grna['start']}-{grna['end']})")

        # Save outputs
        output_dir = Path("results")
        output_dir.mkdir(exist_ok=True)
        (output_dir / "grnas.txt").write_text("\n".join(g['spacer'] for g in grnas))
        print(f"\n📁 Results saved to {output_dir}/")

    except Exception as e:
        print(f"\n❌ Pipeline Error: {e}")
        print("💡 Debug Checklist:")
        print("1. Verify NCBI_API_KEY_2025 environment variable")
        print("2. Check MUSCLE installation: muscle -version")
        print("3. Inspect input files in alignments/ directory")
        sys.exit(1)

if __name__ == "__main__":
    main()
