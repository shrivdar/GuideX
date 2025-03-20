import matplotlib.pyplot as plt
import numpy as np

# Sample data structure
regions = {
    'RdRp III': (5200, 5262, 97.8),
    'CP N-term': (8100, 8123, 96.2),
    'RdRp V': (5500, 5542, 95.5),
    'CP Central': (8340, 8362, 96.7)
}

# Plotting
fig, ax = plt.subplots(figsize=(12, 6))
x_pos = np.arange(len(regions))
colors = ['#4e79a7' if v[2] > 95 else '#a0a0a0' for v in regions.values()]

bars = ax.bar(x_pos, [v[2] for v in regions.values()], 
              color=colors, edgecolor='black', width=0.7)

ax.set_xticks(x_pos)
ax.set_xticklabels([k for k in regions.keys()], rotation=45, ha='right')
ax.set_ylabel('Conservation (%)', fontsize=12)
ax.set_title('Highly Conserved Regions in CBSV Genome', fontsize=14)
ax.spines[['top', 'right']].set_visible(False)

# Add genomic positions below labels
for i, (name, (start, end, _)) in enumerate(regions.items()):
    plt.text(i, -5, f"{start}-{end}", ha='center', fontsize=9)

plt.tight_layout()
plt.savefig('conservation.png', dpi=300, bbox_inches='tight')
