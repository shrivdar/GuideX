from guidex.genome_fetcher import GenomeFetcher
from guidex.conservation import ConservationAnalyzer
from guidex.grna_designer import GuideXGrnaDesigner

# Fetch genomes (example: Influenza A virus)
fetcher = GenomeFetcher(email="darsh.shri123@gmail.com")  # Use your NCBI email!
genomes = fetcher.fetch_ncbi("Influenza A virus[Organism]", limit=5)

# Align genomes and find conserved regions
conservator = ConservationAnalyzer(window_size=30)
aligned_file = conservator.align_genomes(genomes)
jsd_scores = conservator.calculate_jsd(aligned_file)
conserved_regions = [(i, i+30) for i, score in enumerate(jsd_scores) if score > 0.8]

# Design gRNAs
designer = GuideXGrnaDesigner(subtype="LwaCas13a")
grnas = designer.design(genomes[0].seq, conserved_regions)

print(f"Designed {len(grnas)} gRNAs!")
