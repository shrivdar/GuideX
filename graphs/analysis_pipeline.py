"""
COMPLETE VISUALIZATION PIPELINE FOR GUIDEX
Generates 8 essential analysis graphs with cross-verification
"""

import sys
import logging
from pathlib import Path
import pandas as pd
import numpy as np

# Configure local imports
sys.path.append(str(Path(__file__).parent))
from guidex.viz import (
    conservation, comparative, selection, 
    features, plant_health, suppression
)
from guidex.analysis import conservation as cons_analysis

# Configuration
RESULTS_DIR = Path("results")
FIGURES = {
    # Conservation Analysis
    'conservation_profile': '1_conservation_profile.png',
    'venn_diagram': '2_species_venn.png',
    'heatmap': '3_conservation_heatmap.png',
    'violin_plot': '4_conservation_distribution.png',
    
    # gRNA Design
    'selection_landscape': '5_grna_selection.png',
    'feature_importance': '6_nn_features.png',
    
    # Experimental Results
    'plant_health': '7_plant_health.png',
    'viral_suppression': '8_viral_timeline.png'
}

REQUIRED_INPUTS = {
    'conserved_regions': RESULTS_DIR/"conserved_regions.csv",
    'grna_candidates': RESULTS_DIR/"grna_candidates.tsv",
    'alignment': RESULTS_DIR/"alignment.fasta",
    'plant_data': RESULTS_DIR/"plant_metrics.csv",
    'viral_data': RESULTS_DIR/"viral_loads.csv",
    'nn_model': RESULTS_DIR/"models/nn_model.h5"
}

def validate_inputs():
    """Ensure all required files exist"""
    missing = [str(p) for p in REQUIRED_INPUTS.values() if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required inputs: {', '.join(missing)}")
    return {k: pd.read_csv(v) if v.suffix == '.csv' else v 
            for k, v in REQUIRED_INPUTS.items()}

def generate_conservation_figures(data):
    """Generate all 4 conservation-related figures"""
    logging.info("Creating conservation visualizations...")
    
    # 1. Conservation Profile
    cons_scores = data['conserved_regions']['score'].values
    regions = [tuple(map(int, r[1:-1].split(','))) 
              for r in data['conserved_regions']['regions']]
    invalid = [(9000, 9030), (9030, 9060)]  # From pipeline warnings
    
    fig1 = conservation.plot_conservation(cons_scores, regions, invalid)
    fig1.savefig(RESULTS_DIR/FIGURES['conservation_profile'])
    
    # 2. Venn Diagram
    species_data = cons_analysis.calculate_species_overlaps(REQUIRED_INPUTS['alignment'])
    fig2 = comparative.plot_venn_diagram(species_data)
    fig2.savefig(RESULTS_DIR/FIGURES['venn_diagram'])
    
    # 3/4. Cross-isolate Analysis
    genomes = cons_analysis.fetch_genome_sequences(cons_analysis.CBSV_ISOLATES)
    blast_results = cons_analysis.analyze_conservation(
        {r['id']: r for r in data['conserved_regions'].to_dict('records')},
        genomes
    )
    matrix = cons_analysis.generate_conservation_matrix(blast_results)
    
    fig3 = comparative.plot_conservation_heatmap(matrix)
    fig3.savefig(RESULTS_DIR/FIGURES['heatmap'])
    
    fig4 = comparative.plot_conservation_distribution(matrix)
    fig4.savefig(RESULTS_DIR/FIGURES['violin_plot'])

def generate_grna_figures(data):
    """Generate gRNA design figures"""
    logging.info("Creating gRNA visualizations...")
    
    # 5. Selection Landscape
    grna_df = selection.extract_grna_metrics(data['grna_candidates'])
    fig5 = selection.plot_selection_landscape(grna_df)
    fig5.savefig(RESULTS_DIR/FIGURES['selection_landscape'])
    
    # 6. Feature Importance
    model = features.load_model(str(REQUIRED_INPUTS['nn_model']))
    importance = features.calculate_shap_importance(model, grna_df)
    fig6 = features.plot_feature_importance(importance)
    fig6.savefig(RESULTS_DIR/FIGURES['feature_importance'])

def generate_experimental_figures(data):
    """Generate experimental result figures"""
    logging.info("Creating experimental visualizations...")
    
    # 7. Plant Health Radar
    health_data = plant_health.prepare_health_data(
        data['plant_data']['control'].iloc[0],
        data['plant_data']['treated'].iloc[0]
    )
    fig7 = plant_health.plot_plant_health_radar(*health_data)
    fig7.savefig(RESULTS_DIR/FIGURES['plant_health'])
    
    # 8. Viral Suppression
    timeline = suppression.process_timeline_data(data['viral_data'])
    fig8 = suppression.plot_suppression_timeline(timeline)
    fig8.savefig(RESULTS_DIR/FIGURES['viral_suppression'])

def main():
    """Guaranteed visualization pipeline"""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    
    try:
        # 1. Input Validation
        inputs = validate_inputs()
        
        # 2. Generate All Figures
        RESULTS_DIR.mkdir(exist_ok=True)
        
        generate_conservation_figures(inputs)
        generate_grna_figures(inputs)
        generate_experimental_figures(inputs)
        
        # 3. Verification
        existing = [f.name for f in RESULTS_DIR.glob("*.png")]
        missing = [f for f in FIGURES.values() if f not in existing]
        
        if missing:
            raise RuntimeError(f"Failed to generate: {', '.join(missing)}")
            
        logging.info(f"Success! Generated {len(FIGURES)} figures in {RESULTS_DIR}")

    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
