from Bio.SeqRecord import SeqRecord
from pathlib import Path
from guidex.genome_fetcher import GenomeFetcher
from guidex.conservation import ConservationAnalyzer
from guidex.grna_designer import GuideXGrnaDesigner

def main():
    try:
        # 1. Fetch genomes
        print("🕵️ Fetching genomes from NCBI...")
        fetcher = GenomeFetcher(email="darsh.shri123@gmail.com")
        genomes = fetcher.fetch_ncbi("Influenza A virus[Organism]", limit=5)
        
        if not genomes:
            raise ValueError("No genomes found! Check your search term.")
        print(f"✅ Retrieved {len(genomes)} genomes")

        # 2. Find conserved regions
        print("\n🧬 Identifying conserved regions...")
        conservator = ConservationAnalyzer(window_size=30)
        
        # Create alignment directory
        align_dir = Path("alignments")
        align_dir.mkdir(exist_ok=True)
        
        aligned_file = conservator.align_genomes(genomes, output_dir=align_dir)
        jsd_scores = conservator.calculate_jsd(aligned_file)
        
        # Find regions with JSD > 0.8 (80% conservation)
        conserved_regions = [(i, i+30) for i, score in enumerate(jsd_scores) if score > 0.8]
        
        if not conserved_regions:
            raise ValueError("No conserved regions found! Try adjusting window_size or threshold.")
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
        print(f"\n❌ Error: {e}")
        print("💡 Debugging tips:")
        print("- Verify NCBI API access with 'darsh.shri123@gmail.com'")
        print("- Check MAFFT installation: 'mafft --version'")
        print("- Inspect alignment files in 'alignments/' directory")

if __name__ == "__main__":
    main()
