plt.figure(figsize=(10, 8))
sns.set(style="whitegrid", palette="muted")

# Plot all gRNAs
scatter = plt.scatter(df['Predicted_Efficiency'], 
                     df['Validated_Efficiency'],
                     c=df['Predicted_Efficiency'], 
                     cmap='viridis',
                     alpha=0.6,
                     edgecolors='w',
                     s=80)

# Highlight selected
plt.scatter(selected_grnas['Predicted_Efficiency'],
           selected_grnas['Validated_Efficiency'],
           edgecolor='red',
           facecolor='none',
           s=200,
           linewidth=2,
           label='Selected gRNAs')

# Formatting
plt.plot([0, 100], [0, 100], 'k--', alpha=0.3)
plt.colorbar(scatter).set_label('Predicted Efficiency (%)')
plt.xlabel('Computational Prediction (%)', fontsize=12)
plt.ylabel('Experimental Validation (%)', fontsize=12)
plt.title('gRNA Efficiency Validation', fontsize=14)
plt.legend()

plt.savefig('grna_scatter.png', dpi=300)
