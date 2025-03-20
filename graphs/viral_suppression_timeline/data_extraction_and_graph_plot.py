import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def process_timeline_data(time_data: dict) -> dict:
    """
    Process viral load time series data
    Returns:
        processed = {
            'days': list of days post-treatment,
            'control': {
                'mean': [],
                'ci': [],
            },
            'treated': {
                'mean': [],
                'ci': [],
            }
        }
    """
    processed = {'days': sorted(time_data.keys())}
    
    for group in ['control', 'treated']:
        means = []
        cis = []
        
        for day in processed['days']:
            values = time_data[day][group]
            mean = np.mean(values)
            sem = stats.sem(values)
            ci = sem * 1.96  # 95% CI
            
            means.append(mean)
            cis.append(ci)
            
        processed[group] = {'mean': means, 'ci': cis}
    
    return processed

def plot_suppression_timeline(processed: dict):
    """
    Create viral load timeline with confidence bands
    """
    with plt.rc_context(POSTER_STYLE):
        fig, ax = plt.subplots(figsize=(12, 6))
        days = processed['days']
        
        # Plot control
        ax.errorbar(days, processed['control']['mean'], 
                   yerr=processed['control']['ci'],
                   fmt='-o', color='#d7191c', capsize=5,
                   label='Control', linewidth=2)
        
        # Plot treated
        ax.errorbar(days, processed['treated']['mean'],
                   yerr=processed['treated']['ci'],
                   fmt='-o', color='#2c7bb6', capsize=5,
                   label='CRISPR-treated', linewidth=2)
        
        # Formatting
        ax.set_xlabel("Days Post-Treatment", labelpad=10)
        ax.set_ylabel("Viral Load (log10 copies/μg RNA)", labelpad=10)
        ax.set_title("Viral Suppression Timeline", pad=20)
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
        
        # Infection phase annotations
        ax.text(7, ax.get_ylim()[1]*0.9, 'Acute Phase', 
               color='#636363', fontstyle='italic')
        ax.text(21, ax.get_ylim()[1]*0.5, 'Suppression Phase',
               color='#636363', fontstyle='italic')
        
        # Panel label
        ax.text(0.01, 0.95, 'F', transform=ax.transAxes,
               fontsize=24, weight='bold', va='top')
        
        # Significance markers
        for i, day in enumerate(days):
            if processed['control']['mean'][i] - processed['treated']['mean'][i] > 1:
                ax.text(day, processed['control']['mean'][i]+0.5, '***',
                       ha='center', color='#2ca25f', fontsize=14)
        
        plt.tight_layout()
        return fig
