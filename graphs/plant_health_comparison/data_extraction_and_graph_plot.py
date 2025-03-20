import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def prepare_health_data(control_data: dict, treated_data: dict) -> tuple:
    """
    Process plant health metrics for radar chart
    Returns:
        - categories: List of metric names
        - control_values: Normalized values [0-1]
        - treated_values: Normalized values [0-1]
        - actual_ranges: Dictionary of original value ranges
    """
    categories = list(control_data.keys())
    
    # Normalize data 0-1 within each category
    control_values = []
    treated_values = []
    actual_ranges = {}
    
    for cat in categories:
        c_val = control_data[cat]
        t_val = treated_data[cat]
        
        # Store original ranges
        actual_ranges[cat] = {
            'control': (min(c_val), max(c_val)),
            'treated': (min(t_val), max(t_val))
        }
        
        # Normalize
        all_vals = np.concatenate([c_val, t_val])
        c_norm = np.mean(c_val) / max(all_vals)
        t_norm = np.mean(t_val) / max(all_vals)
        
        control_values.append(c_norm)
        treated_values.append(t_norm)
    
    return categories, control_values, treated_values, actual_ranges

def plot_plant_health_radar(categories: list, 
                          control_values: list,
                          treated_values: list,
                          actual_ranges: dict):
    """
    Create radar chart comparing plant health metrics
    """
    with plt.rc_context(POSTER_STYLE):
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, polar=True)
        
        # Calculate angles
        N = len(categories)
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]
        
        # Complete the circular plot
        control = control_values + control_values[:1]
        treated = treated_values + treated_values[:1]
        
        # Plot
        ax.plot(angles, control, color='#d7191c', linewidth=2, linestyle='solid', label='Control')
        ax.fill(angles, control, color='#d7191c', alpha=0.25)
        ax.plot(angles, treated, color='#2c7bb6', linewidth=2, linestyle='solid', label='CRISPR-treated')
        ax.fill(angles, treated, color='#2c7bb6', alpha=0.25)
        
        # Axis settings
        ax.set_theta_offset(np.pi/2)
        ax.set_theta_direction(-1)
        plt.xticks(angles[:-1], categories)
        ax.tick_params(axis='x', pad=25)
        ax.set_rlabel_position(0)
        plt.yticks([0.25, 0.5, 0.75], ["25%", "50%", "75%"], color="grey", size=12)
        plt.ylim(0, 1)
        
        # Legend and title
        plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
        plt.title("Plant Health Metrics Comparison\n(Normalized Values)", pad=35)
        
        # Range annotations
        for i, (cat, vals) in enumerate(actual_ranges.items()):
            ax.text(angles[i], 1.08, 
                   f"C: {vals['control'][0]:.1f}-{vals['control'][1]:.1f}\n"
                   f"T: {vals['treated'][0]:.1f}-{vals['treated'][1]:.1f}",
                   ha='center', va='center', fontsize=10)
        
        # Panel label
        ax.text(0.05, 0.95, 'E', transform=ax.transAxes,
               fontsize=24, weight='bold', va='top')
        
        return fig
