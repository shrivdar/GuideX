from conservation_viz import plot_conservation
from grna_viz import plot_grna_selection
# ... other imports

def generate_all_figures():
    plot_conservation()
    plot_grna_selection()
    # ... other plots
    print("All figures saved to /results")

if __name__ == "__main__":
    generate_all_figures()
