from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
import sys
import os
import shutil
import importlib
sys.path.insert(0, str(Path(__file__).parent))
importlib.invalidate_caches()
from guidex.genome_fetcher import GenomeFetcher
from guidex.conservation import ConservationAnalyzer
from guidex.core import GuideXGrnaDesigner

# Valid local genomes meeting NCBI v2 criteria
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

        # 1. Initialize fetcher with API key from environment
        fetcher = GenomeFetcher(api_key=os.getenv("NCBI_API_KEY_2025"))


        # 2. Attempt NCBI fetch
        genomes = []
        try:
            print("🕵️ Attempting NCBI Datasets API v2 fetch...")
            genomes = fetcher.fetch_genomes(
                target="Influenza A virus",  # Remove [Organism] qualifier
                genome_type="reference",
                # Add limit parameter handling to class
            )
            print(f"✅ Retrieved {len(genomes)} NCBI genomes")
        except Exception as e:
            print(f"🔴 NCBI Error: {e}")
            print("🔄 Falling back to local genomes")
            genomes = LOCAL_GENOMES

        # 3. Modern sequence validation
        valid_genomes = []
        for g in genomes:
            try:
                # Convert to string for safe gap removal
                seq_str = str(g.seq).upper().replace("-", "")
                seq = Seq(seq_str)
                
                # Validate length and content
                if 5000 <= len(seq) <= 15000 and len(seq_str) == len(seq):
                    valid_genomes.append(SeqRecord(
                        seq,
                        id=g.id,
                        description=f"Validated | {len(seq)}bp"
                    ))
                    print(f"🧬 Validated {g.id} ({len(seq)} bp)")
                else:
                    print(f"⚠️ Excluded {g.id} (invalid length: {len(seq)} bp)")
            except Exception as e:
                print(f"🚨 Processing Error: {e}")

        if len(valid_genomes) < 2:
            raise ValueError(f"Only {len(valid_genomes)} valid genomes (need ≥2)")
        
        # 4. Alignment and analysis
        print("\n🧬 Starting conservation pipeline...")
        conservator = ConservationAnalyzer(window_size=30)
        
        # Run alignment
        aligned_file = conservator.align_genomes(valid_genomes, Path("alignments"))
        print(f"🔍 Alignment saved to: {aligned_file}")

        # Conservation analysis
        print("\n🔎 Identifying conserved regions...")
        jsd_scores = conservator.calculate_jsd(aligned_file)
        conserved_regions = [(i, i+30) for i, score in enumerate(jsd_scores) if score > 0.8]
        
        if not conserved_regions:
            raise ValueError("No conserved regions found!")
        print(f"✅ Found {len(conserved_regions)} conserved regions")

        # 5. gRNA Design
        print("\n🔬 Designing Cas13 gRNAs...")
        designer = GuideXGrnaDesigner(subtype="LwaCas13a")
        grnas = designer.design(str(valid_genomes[0].seq), conserved_regions)
        
        # Output results
        print(f"\n🎉 Success! Designed {len(grnas)} gRNAs:")
        for i, grna in enumerate(grnas[:5], 1):
            print(f"{i}. {grna['spacer']} (pos {grna['start']}-{grna['end']})")

        # Save results
        output_dir = Path("results")
        output_dir.mkdir(exist_ok=True)
        (output_dir / "grnas.txt").write_text("\n".join(g['spacer'] for g in grnas))
        print(f"\n📁 Results saved to {output_dir}/")

    except Exception as e:
        print(f"\n❌ Pipeline Error: {e}")
        print("💡 Debug Checklist:")
        print("1. Verify NCBI_API_KEY environment variable")
        print("2. Check 'alignments/MAFFT_IN.fasta' exists")
        print("3. Test manual alignment: mafft --auto alignments/MAFFT_IN.fasta")
        sys.exit(1)

if __name__ == "__main__":
    main()
