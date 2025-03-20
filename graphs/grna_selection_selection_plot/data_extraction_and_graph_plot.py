
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def extract_grna_metrics(grna_data_path: str) -> pd.DataFrame:
    """
    Process gRNA candidate metrics from selection pipeline
    Returns DataFrame with columns:
    ['sequence', 'efficiency', 'conservation', 'offtargets', 'gc_content', 'selected']
    """
    df = pd.read_csv(grna_data_path, sep='\t')
    
    # Calculate composite metrics
    df['efficiency'] = df['cutting_frequency'] * df['accessibility']
    df['offtargets'] = np.log10(df['exact_matches'] + 0.1)
    df['gc_content'] = df['sequence'].apply(lambda x: (x.count('G')+x.count('C'))/len(x)*100
    
    # Normalize conservation scores 0-100
    df['conservation'] = (df['conservation_score'] - df['conservation_score'].min()) / \
                        (df['conservation_score'].max() - df['conservation_score'].min()) * 100
    
    return df[['sequence', 'efficiency', 'conservation', 'offtargets', 'gc_content', 'selected']]

def plot_selection_landscape(df: pd.DataFrame):
    with plt.rc_context(POSTER_STYLE):
        fig, ax = plt.subplots(figsize=(16, 10))
        
        # Hexbin background
        hexbin = ax.hexbin(df['efficiency'], df['conservation'], 
                          gridsize=40, cmap='Blues', norm=LogNorm(),
                          mincnt=1, alpha=0.8, edgecolors='face')
        
        # Selected candidates
        selected = df[df['selected']]
        scatter = ax.scatter(selected['efficiency'], selected['conservation'],
                            c=selected['offtargets'], s=selected['gc_content']*15,
                            cmap='Reds_r', ec='k', lw=0.8, alpha=0.9,
                            vmin=0, vmax=3, zorder=10)
        
        # Optimal zone
        ax.axvspan(0.7, 1.0, ymin=0.8, ymax=1.0, color='#2ca25f', alpha=0.2, lw=0)
        ax.plot([0.7, 1.0], [80, 80], color='#2ca25f', ls='--', lw=2, alpha=0.8)
        ax.plot([0.7, 0.7], [80, 100], color='#2ca25f', ls='--', lw=2, alpha=0.8)
        
        # Colorbars
        cb1 = fig.colorbar(hexbin, ax=ax, label='Candidate Density')
        cb2 = fig.colorbar(scatter, ax=ax, label='Off-target Score (log10)')
        
        # Legend for GC content
        gc_legend = [
            plt.Line2D([0], [0], marker='o', color='w', 
                      markersize=np.sqrt(40*15), markeredgecolor='k',
                      markerfacecolor='gray', label='40% GC'),
            plt.Line2D([0], [0], marker='o', color='w', 
                      markersize=np.sqrt(60*15), markeredgecolor='k',
                      markerfacecolor='gray', label='60% GC')
        ]
        ax.legend(handles=gc_legend, title='GC Content', loc='lower right')
        
        # Labels and titles
        ax.set_xlabel("Editing Efficiency Prediction", labelpad=10)
        ax.set_ylabel("Sequence Conservation (%)", labelpad=10)
        ax.set_title("gRNA Selection Landscape\n(Optimal Zone Highlighted)", pad=20)
        
        # Panel label
        ax.text(0.01, 0.97, 'C', transform=ax.transAxes,
               fontsize=24, weight='bold', va='top')
        
        return fig
