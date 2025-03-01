from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
from guidex.conservation import ConservationAnalyzer
from guidex.grna_designer import GuideXGrnaDesigner
import sys

# Predefined influenza genome segments (HA and NA genes)
TEST_GENOMES = [
    SeqRecord(
        Seq("ATGCGATAGCATCGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA" * 100 +
            "ATGCGATAGCATCGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA" * 100),
        id="HA_Genome_1",
        description="Hemagglutinin segment | Length: 1800bp"
    ),
    SeqRecord(
        Seq("ATGCGATAGCATGGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA" * 100 +
            "ATGCGATAGCATCGACTAGCATGACGTACGTACGTACGTACGTACGTACGTACGTA" * 100),
        id="NA_Genome_2",
        description="Neuraminidase segment | Length: 1800bp"
    )
]

def main():
    try:
        print("🔬 Running in LOCAL TEST MODE")
        
        # 1. Use predefined genomes
        genomes = TEST_GENOMES
        print(f"✅ Loaded {len(genomes)} test genomes")
        print(f"📏 Genome lengths: {[len(g.seq) for g in genomes]} bp")

        # 2. Alignment and Conservation Analysis
        print("\n🧬 Aligning genomes...")
        align_dir = Path("alignments").absolute()
        align_dir.mkdir(exist_ok=True)
        
        conservator = ConservationAnalyzer(window_size=30)
        aligned_file = conservator.align_genomes(genomes, align_dir)
        print(f"🔍 Alignment saved to: {aligned_file}")

        # Verify alignment file exists
        if not aligned_file.exists():
            raise FileNotFoundError(f"Alignment failed: {aligned_file}")

        # 3. Conservation analysis
        print("\n🔎 Identifying conserved regions...")
        jsd_scores = conservator.calculate_jsd(aligned_file)
        conserved_regions = [(i, i+30) for i, score in enumerate(jsd_scores) if score > 0.8]
        
        if not conserved_regions:
            raise ValueError("No conserved regions found!")
        print(f"✅ Found {len(conserved_regions)} conserved regions")

        # 4. gRNA Design
        print("\n🔬 Designing Cas13 gRNAs...")
        designer = GuideXGrnaDesigner(subtype="LwaCas13a")
        grnas = designer.design(str(genomes[0].seq), conserved_regions)
        
        print(f"\n🎉 Success! Designed {len(grnas)} gRNAs:")
        for i, grna in enumerate(grnas[:5], 1):
            print(f"{i}. {grna['spacer']} (pos {grna['start']}-{grna['end']})")

        # 5. Save results
        output_file = Path("results/grnas.txt")
        output_file.parent.mkdir(exist_ok=True)
        with open(output_file, "w") as f:
            f.write("\n".join(g['spacer'] for g in grnas))
        print(f"\n📁 Results saved to {output_file}")

    except Exception as e:
        print(f"\n❌ Error: {e}")
        print("💡 Debug Tips:")
        print("1. Check alignment input: cat alignments/MAFFT_IN.fasta")
        print("2. Verify MAFFT installation: mafft --version")
        sys.exit(1)

if __name__ == "__main__":
    main()
