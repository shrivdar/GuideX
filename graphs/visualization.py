# guidex_visualizations.py
"""
ONE-FILE VISUALIZATION SOLUTION
Run with: python guidex_visualizations.py
Outputs: All graphs as .png files in /results
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib_venn import venn2
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from scipy import stats

# Configuration
RESULTS_DIR = "results"
plt.style.use('ggplot')  # Professional style

def fake_conservation_data():
    """Mock data - replace with your real conservation scores"""
    return {
        'positions': np.arange(0, 10000, 50),
        'scores': np.random.uniform(60, 95, 200),
        'regions': [(200, 400), (600, 800), (1500, 1800)],
        'invalid': [(9000, 9030), (9030, 9060)],
        'species_sets': {
            'CBSV_UG': set(range(150, 400)),
            'CBSV_TZ': set(range(200, 500)),
            'CBSV_KE': set(range(100, 300))
        }
    }

def plot_conservation():
    """Conservation bar chart with regions"""
    data = fake_conservation_data()
    
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.bar(data['positions'], data['scores'], width=50, color='royalblue')
    
    # Highlight regions
    for start, end in data['regions']:
        ax.axvspan(start, end, color='forestgreen', alpha=0.3)
        
    for start, end in data['invalid']:
        ax.axvspan(start, end, color='red', alpha=0.2, hatch='//')
        
    ax.set_title("Genome Conservation Analysis")
    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("Conservation Score")
    fig.savefig(f"{RESULTS_DIR}/conservation.png")

def plot_venn():
    """Cross-species Venn diagram"""
    data = fake_conservation_data()
    
    plt.figure(figsize=(8, 8))
    venn2([data['species_sets']['CBSV_UG'], 
          data['species_sets']['CBSV_TZ']],
          set_labels=('Uganda', 'Tanzania'))
    
    plt.title("Shared Conserved Regions")
    plt.savefig(f"{RESULTS_DIR}/venn.png")

def plot_selection():
    """gRNA selection landscape"""
    # Mock gRNA data
    grnas = pd.DataFrame({
        'efficiency': np.random.normal(0.7, 0.1, 200)),
        'conservation': np.random.uniform(60, 95, 200)),
        'offtargets': np.random.lognormal(1, 0.3, 200)),
        'gc': np.random.uniform(30, 60, 200)
    })
    
    fig, ax = plt.subplots(figsize=(10, 6))
    scatter = ax.scatter(
        grnas['efficiency'], 
        grnas['conservation'],
        c=np.log10(grnas['offtargets']),
        s=grnas['gc'],
        cmap='viridis'
    )
    
    plt.colorbar(scatter, label='Log10 Off-targets')
    ax.set_title("gRNA Selection Landscape\n(Size=GC%, Color=Off-targets)")
    ax.set_xlabel("Efficiency Score")
    ax.set_ylabel("Conservation (%)")
    fig.savefig(f"{RESULTS_DIR}/selection.png")

def plot_neural_network():
    """Feature importance plot"""
    features = ['GC Content', 'Thermodynamics', 'Conservation', 'Off-target']
    importance = [0.28, 0.22, 0.19, -0.15]
    
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.barh(features, importance, color='darkorange')
    ax.set_title("Model Feature Importance")
    ax.set_xlabel("Impact Score")
    fig.savefig(f"{RESULTS_DIR}/nn_features.png")

def plot_plant_health():
    """Radar chart for plant metrics"""
    categories = ['Height', 'Yield', 'Health', 'Disease', 'Chlorophyll']
    control = [3, 2, 4, 4, 3]
    treated = [4, 4, 5, 2, 4]
    
    angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False).tolist()
    
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'polar': True})
    ax.plot(angles + angles[:1], control + control[:1], label='Control')
    ax.plot(angles + angles[:1], treated + treated[:1], label='Treated')
    ax.fill(angles + angles[:1], control + control[:1], alpha=0.2)
    ax.fill(angles + angles[:1], treated + treated[:1], alpha=0.2)
    
    ax.set_xticks(angles)
    ax.set_xticklabels(categories)
    ax.set_title("Plant Health Comparison")
    ax.legend()
    fig.savefig(f"{RESULTS_DIR}/plant_health.png")

def plot_viral_suppression():
    """Viral load timeline"""
    days = [0, 3, 7, 14, 21, 28]
    control = [5.2, 6.1, 7.2, 7.8, 8.0, 8.2]
    treated = [5.1, 4.8, 3.9, 2.1, 1.5, 1.2]
    
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(days, control, '-o', label='Control')
    ax.plot(days, treated, '-o', label='Treated')
    
    ax.set_title("Viral Load Over Time")
    ax.set_xlabel("Days Post-Treatment")
    ax.set_ylabel("log10 Viral Load")
    ax.legend()
    ax.grid(True)
    fig.savefig(f"{RESULTS_DIR}/viral_suppression.png")

if __name__ == "__main__":
    # Create output directory
    import os
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    # Generate all plots
    plot_conservation()
    plot_venn()
    plot_selection()
    plot_neural_network()
    plot_plant_health()
    plot_viral_suppression()
    
    print(f"All graphs saved to {RESULTS_DIR}/ directory!")
