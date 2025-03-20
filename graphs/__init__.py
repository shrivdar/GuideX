# guidex/viz/__init__.py
import matplotlib.pyplot as plt

# Global style configuration
plt.style.use('seaborn-whitegrid')
POSTER_STYLE = {
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.titlesize": 18,
    "font.family": "DejaVu Sans",
    "figure.dpi": 300,
    "savefig.bbox": "tight",
    "axes.prop_cycle": plt.cycler(color=plt.cm.Dark2.colors),
    "image.cmap": "viridis"
}
