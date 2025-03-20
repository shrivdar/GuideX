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
    with plt.rc_context(POSTER_STYLE):
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Color scheme matching conservation plot
        colors = ['#2c7bb6', '#fdae61', '#d7191c']
        alpha = 0.4
        
        if len(primary_species) == 2:
            venn = venn2([species_sets[s] for s in primary_species],
                        set_labels=primary_species,
                        set_colors=colors[:2], alpha=alpha)
        elif len(primary_species) == 3:
            venn = venn3([species_sets[s] for s in primary_species],
                        set_labels=primary_species,
                        set_colors=colors, alpha=alpha)
            
        # Style adjustments
        for text in venn.subset_labels:
            if text: 
                text.set_fontsize(14)
                text.set_fontweight('bold')
                
        for label in venn.set_labels:
            if label: 
                label.set_fontsize(16)
                label.set_fontweight('semibold')
                
        # Panel label
        ax.text(-0.1, 0.95, 'B', transform=ax.transAxes,
               fontsize=24, weight='bold', va='top')
        
        # Annotation box
        total_regions = sum(len(s) for s in species_sets.values())
        annotation_text = (
            f"Total Conserved Regions: {total_regions:,}\n"
            f"Species: {', '.join(primary_species)}"
        )
        ax.text(0.95, 0.05, annotation_text,
               transform=ax.transAxes, ha='right', va='bottom',
               bbox=dict(facecolor='white', alpha=0.8))
        
        plt.title("Cross-Species Conservation Overlap", pad=20)
        return fig
