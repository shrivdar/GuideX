from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq  # Added import
from pathlib import Path
from guidex.genome_fetcher import GenomeFetcher
from guidex.conservation import ConservationAnalyzer
from guidex.grna_designer import GuideXGrnaDesigner

def main():
    try:
        # 1. Fetch genomes with Seq conversion
        print("🕵️ Fetching genomes from NCBI...")
        fetcher = GenomeFetcher(email="darsh.shri123@gmail.com")
        raw_genomes = fetcher.fetch_ncbi("Influenza A virus[Organism]", limit=5)
        
        # Convert to proper SeqRecords
        genomes = []
        for g in raw_genomes:
            try:
                seq = Seq(str(g.seq).upper())  # Force uppercase
                seq = seq.ungap()  # Remove gaps
                if len(seq) < 100:  # Minimum length check
                    continue
                genomes.append(SeqRecord(seq, id=g.id, description=""))
            except Exception as e:
                print(f"⚠️ Invalid genome {g.id}: {str(e)}")
                continue

        if len(genomes) < 2:
            raise ValueError("Need at least 2 valid genomes for alignment")
        
        if not genomes:
            raise ValueError("No genomes found! Check your search term.")
        print(f"✅ Retrieved {len(genomes)} genomes")

        # 2. Find conserved regions (now uses fixed align_genomes)
        print("\n🧬 Identifying conserved regions...")
        conservator = ConservationAnalyzer(window_size=30)
        align_dir = Path("alignments")
        aligned_file = conservator.align_genomes(genomes, output_dir=align_dir)

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
        print(f"\n❌ Error: {e}")
        print("💡 Debugging tips:")
        print("- Verify NCBI API access with 'darsh.shri123@gmail.com'")
        print("- Check MAFFT installation: 'mafft --version'")
        print("- Inspect alignment files in 'alignments/' directory")

if __name__ == "__main__":
    main()
