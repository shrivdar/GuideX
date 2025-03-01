from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
from guidex.genome_fetcher import GenomeFetcher
from guidex.conservation import ConservationAnalyzer
from guidex.grna_designer import GuideXGrnaDesigner
import sys

def main():
    try:
        # 1. Genome fetching with validation
        print("🕵️ Fetching genomes from NCBI...")
        fetcher = GenomeFetcher(email="darsh.shri123@gmail.com")
        raw_genomes = fetcher.fetch_ncbi("Influenza A virus[Organism]", limit=5)
        
        # Convert to valid SeqRecords
        genomes = []
        for g in raw_genomes:
            try:
                seq = Seq(str(g.seq).upper().ungap()
                if len(seq) < 1000:  # Adjusted minimum for viruses
                    continue
                genomes.append(SeqRecord(seq, id=g.id, description=""))
            except Exception as e:
                print(f"⚠️ Invalid genome {g.id}: {str(e)}", file=sys.stderr)

        # Validate genome count
        if len(genomes) < 2:
            raise ValueError(f"Only {len(genomes)} valid genomes (need ≥2)")
        print(f"✅ Retrieved {len(genomes)} genomes (first 50bp: {genomes[0].seq[:50]}...)")

        # 2. Alignment and Conservation Analysis
        print("\n🧬 Aligning genomes...")
        align_dir = Path("alignments").absolute()
        conservator = ConservationAnalyzer(window_size=30)
        aligned_file = conservator.align_genomes(genomes, align_dir)
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
        print(f"\n❌ Critical Error: {e}", file=sys.stderr)
        print("💡 Debug Checklist:", file=sys.stderr)
        print("1. Verify 'alignments/MAFFT_IN.fasta' exists", file=sys.stderr)
        print("2. Check file contents: cat alignments/MAFFT_IN.fasta", file=sys.stderr)
        print("3. Test manual alignment: mafft --auto alignments/MAFFT_IN.fasta", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
