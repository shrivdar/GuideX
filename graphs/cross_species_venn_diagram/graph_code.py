from matplotlib_venn import venn2

plt.figure(figsize=(8, 8))
venn2(subsets=(len(cbsv_targets), len(ucbsv_seqs), len(matches)),
      set_labels=('CBSV Targets', 'UCBSV Targets'),
      set_colors=('#4e79a7', '#f28e2b'),
      alpha=0.7)

plt.title('Cross-Species Target Compatibility', fontsize=14)
plt.savefig('venn.png', dpi=300)
