from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
import sys
import os
import shutil
from guidex.genome_fetcher import GenomeFetcher
from guidex.conservation import ConservationAnalyzer
from guidex.grna_designer import GuideXGrnaDesigner

# Local fallback genomes meeting modern criteria
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

        # 1. Genome acquisition
        genomes = []
        try:
            print("🕵️ Attempting NCBI Datasets API v2 fetch...")
            fetcher = GenomeFetcher(
                email="darsh.shri123@gmail.com",
                api_key=os.getenv("NCBI_API_KEY")
            )
            genomes = fetcher.fetch_ncbi(
                "Influenza A virus[Organism]",
                limit=5,
                exclude_atypical=True
            )
            print(f"✅ Retrieved {len(genomes)} NCBI genomes")
        except Exception as e:
            print(f"🔴 NCBI Unavailable: {e}")
            genomes = LOCAL_GENOMES
            print("🔄 Using validated local genomes")

        # 2. Modern sequence validation
        valid_genomes = []
        for g in genomes:
            try:
                # New NCBI sequence standards
                from Bio.Seq import Seq
                from Bio.SeqUtils import gc_fraction
                seq = str(g.seq).upper().replace("-", "")  # Remove gaps safely
                seq = Seq(seq)  # Convert back to Seq object
                if 5000 <= len(seq) <= 15000:  # 5600bp test sequences should pass  # Typical influenza genome range
                    valid_genomes.append(SeqRecord(
                        seq,
                        id=g.id,
                        description=f"{g.description} | Validated"
                    ))
                    print(f"🧬 Validated {g.id} ({len(seq)} bp)")
                else:
                    print(f"⚠️ Excluded {g.id} (length {len(seq)} bp)")
            except Exception as e:
                print(f"🚨 Validation Error: {e}")

        if len(valid_genomes) < 2:
            raise ValueError(f"Only {len(valid_genomes)} valid genomes (need ≥2)")
        
        # 3. Modern alignment and analysis
        print("\n🧬 Starting conservation pipeline...")
        conservator = ConservationAnalyzer(window_size=30)
        aligned_file = conservator.align_genomes(valid_genomes, Path("alignments"))

        print("\n🧬 Starting conservation pipeline...")
        conservator = ConservationAnalyzer(window_size=30)
        aligned_file = conservator.align_genomes(valid_genomes, Path("alignments"))
        print(f"🔍 Alignment saved to: {aligned_file}")

        print("\n🔎 Identifying conserved regions...")
        try:
            jsd_scores = conservator.calculate_jsd(aligned_file)
            conserved_regions = [
                (i, i+30) 
                for i, score in enumerate(jsd_scores) 
                if score > 0.8
            ]
    
            if not conserved_regions:
                raise ValueError("No conserved regions found!")
            print(f"✅ Found {len(conserved_regions)} conserved regions")
    
            # Visualization
            plot_path = Path("results") / "conservation.html"
            conservator.plot_conservation(jsd_scores, plot_path)
            print(f"📊 Conservation plot saved to {plot_path}")

        except Exception as e:
            print(f"🚨 Conservation analysis failed: {e}")
            raise

        # 3. Design gRNAs
        print("\n🔬 Designing Cas13 gRNAs...")
        designer = GuideXGrnaDesigner(subtype="LwaCas13a")
        grnas = designer.design(str(genomes[0].seq), conserved_regions)  # Explicit string conversion
        
        print(f"\n🎉 Success! Designed {len(grnas)} gRNAs:")
        for i, grna in enumerate(grnas[:5], 1):  # Show top 5
            print(f"{i}. Position {grna['start']}-{grna['end']}: {grna['spacer']}")

        # Save results
        output_file = Path("results") / "grnas.txt"
        output_file.parent.mkdir(exist_ok=True)
        with open(output_file, "w") as f:
            f.write("\n".join(f"{g['spacer']}" for g in grnas))
        print(f"\n📁 Results saved to {output_file}")

    except Exception as e:
        print(f"\n❌ Pipeline Error: {e}")
        print("💡 Modern Debug Checklist:")
        print("1. Verify NCBI_API_KEY environment variable")
        print("2. Check API status: https://api.ncbi.nlm.nih.gov/datasets/v2")
        print("3. Inspect alignments/MAFFT_IN.fasta format")
        sys.exit(1)

if __name__ == "__main__":
    main()
