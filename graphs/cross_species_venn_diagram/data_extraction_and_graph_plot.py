# guidex/viz/comparative.py
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
from typing import Dict, Set

def calculate_species_overlaps(alignment_path: str, 
                             conservation_threshold: float = 80) -> Dict[str, Set[int]]:
    """
    Identify conserved regions per species
    Returns dictionary mapping species names to sets of conserved positions
    """
    alignment = AlignIO.read(alignment_path, "fasta")
    species_regions = {}
    
    for record in alignment:
        species = record.id.split("|")[0]  # Assuming ID format: "Species|Accession"
        seq = str(record.seq)
        
        # Find conserved regions for this species
        conserved = set()
        in_region = False
        for i, base in enumerate(seq):
            if base != '-' and not in_region:
                start = i
                in_region = True
            elif base == '-' and in_region:
                end = i-1
                if (end - start) >= 20:  # Minimum region length
                    conserved.update(range(start, end+1))
                in_region = False
        
        species_regions[species] = conserved
    
    return species_regions

def plot_venn_diagram(species_sets: Dict[str, Set[int]], 
                     primary_species: List[str] = None):
    """
    Generate Venn diagram for conserved region overlap
    """
    if not primary_species:
        primary_species = list(species_sets.keys())[:3]  # Default to first 3 species
        
    if len(primary_species) == 2:
        venn = venn2([species_sets[s] for s in primary_species],
                    set_labels=primary_species)
    elif len(primary_species) == 3:
        venn = venn3([species_sets[s] for s in primary_species],
                    set_labels=primary_species)
    else:
        raise ValueError("Venn diagram supports 2-3 species comparison")
    
    plt.title("Cross-Species Conservation Overlap", fontsize=14)
    for text in venn.subset_labels:
        if text: text.set_fontsize(12)
    return plt.gcf()
