import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def plot_conservation_heatmap(matrix: pd.DataFrame):
    """Create annotated heatmap of conservation percentages"""
    with plt.rc_context(POSTER_STYLE):
        fig, ax = plt.subplots(figsize=(20, 15))
        
        # Custom colormap
        cmap = LinearSegmentedColormap.from_list(
            'conservation', ['#ffffff', '#2c7bb6', '#d7191c'], N=256
        )
        
        # Plot heatmap
        sns.heatmap(
            matrix,
            cmap=cmap,
            vmin=0,
            vmax=100,
            annot=True,
            fmt=".0f",
            annot_kws={"size": 8},
            linewidths=0.5,
            ax=ax
        )
        
        # Formatting
        ax.set_xlabel("CBSV Isolates", labelpad=15)
        ax.set_ylabel("Conserved Regions", labelpad=15)
        ax.set_title("Cross-Isolate Conservation Analysis", pad=25)
        
        # Rotate labels
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        
        # Panel label
        ax.text(0.01, 0.95, 'G', transform=ax.transAxes,
               fontsize=24, weight='bold', va='top')
        
        return fig

def plot_conservation_distribution(matrix: pd.DataFrame):
    """Violin plot showing conservation distribution"""
    with plt.rc_context(POSTER_STYLE):
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Melt matrix for distribution plotting
        melted = matrix.melt(var_name='Isolate', value_name='Identity')
        
        # Plot
        sns.violinplot(
            data=melted,
            x='Isolate',
            y='Identity',
            palette='Blues',
            inner='quartile',
            cut=0,
            ax=ax
        )
        
        # Formatting
        ax.set_ylim(0, 100)
        ax.set_xlabel("CBSV Isolates", labelpad=15)
        ax.set_ylabel("Conservation Identity (%)", labelpad=15)
        ax.set_title("Conservation Distribution Across Isolates", pad=25)
        plt.xticks(rotation=45, ha='right')
        
        # Panel label
        ax.text(0.01, 0.95, 'H', transform=ax.transAxes,
               fontsize=24, weight='bold', va='top')
        
        return fig
