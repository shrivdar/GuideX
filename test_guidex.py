from Bio.SeqRecord import SeqRecord
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
from guidex.genome_fetcher import GenomeFetcher
from guidex.conservation import ConservationAnalyzer
from guidex.grna_designer import GuideXGrnaDesigner
import sys
import shutil

# Local fallback sequences (HA and NA genes from influenza)
LOCAL_GENOMES = [
    SeqRecord(
        Seq(("ATGCGATAGCATCGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA" * 100) +
            ("ATGCGATAGCATCGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA" * 100)),
        id="Local_HA_1",
        description="Hemagglutinin (Local Backup)"
    ),
    SeqRecord(
        Seq(("ATGCGATAGCATGGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA" * 100) +
            ("ATGCGATAGCATCGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA" * 100)),
        id="Local_NA_2",
        description="Neuraminidase (Local Backup)"
    )
]

def main():
    try:
        # Clear previous runs
        shutil.rmtree("alignments", ignore_errors=True)
        shutil.rmtree("results", ignore_errors=True)
        
        # 1. Attempt NCBI fetch with fallback
        genomes = []
        try:
            print("🕵️ Attempting NCBI genome fetch...")
            fetcher = GenomeFetcher(email="darsh.shri123@gmail.com")
            raw_genomes = fetcher.fetch_ncbi("Influenza A virus[Organism]", limit=5)
            
            # Process NCBI genomes
            for g in raw_genomes:
                try:
                    seq = Seq(str(g.seq).upper().ungap())
                    if len(seq) >= 1000:
                        genomes.append(SeqRecord(seq, id=g.id))
                        print(f"✅ NCBI Genome: {g.id} ({len(seq)} bp)")
                except Exception as e:
                    print(f"⚠️ NCBI Processing Error: {str(e)}")
            
            if len(genomes) < 2:
                raise RuntimeError("Insufficient NCBI genomes")
                
        except Exception as e:
            print(f"🔴 NCBI Unavailable: {str(e)}")
            print("🔄 Using local test genomes")
            genomes = LOCAL_GENOMES

        # 2. Validate genomes
        print(f"\n🔬 Final Genome Set ({len(genomes)} sequences)")
        for g in genomes:
            print(f" - {g.id}: {len(g.seq)} bp")
        
        # 3. Alignment and Analysis
        print("\n🧬 Aligning genomes...")
        conservator = ConservationAnalyzer(window_size=30)
        aligned_file = conservator.align_genomes(genomes, "alignments")
        print(f"🔍 Alignment saved to: {aligned_file}")

        print("\n🔎 Identifying conserved regions...")
        jsd_scores = conservator.calculate_jsd(aligned_file)
        conserved_regions = [(i, i+30) for i, score in enumerate(jsd_scores) if score > 0.8]

        if not conserved_regions:
            raise ValueError("No conserved regions found!")
        print(f"✅ Found {len(conserved_regions)} conserved regions")

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
        print(f"\n❌ Critical Error: {e}")
        print("💡 Immediate Steps:")
        print("1. Check internet connection")
        print("2. Run: mafft --version")
        print("3. Inspect alignments/MAFFT_IN.fasta")
        sys.exit(1)

if __name__ == "__main__":
    main()
