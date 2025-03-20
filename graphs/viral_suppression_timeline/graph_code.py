plt.figure(figsize=(10, 6))
plt.plot(time_points, viral_load['mean'], 
        marker='o', markersize=10,
        color='#59a14f', linewidth=3,
        label='Predicted Suppression')

plt.fill_between(time_points, viral_load['ci_low'], viral_load['ci_high'],
                color='#59a14f', alpha=0.2)

plt.xticks(time_points, ['24h', '48h', '72h'], fontsize=12)
plt.yticks(np.arange(0, 101, 20), fontsize=12)
plt.xlabel('Time Post-Treatment', fontsize=12)
plt.ylabel('Viral Load Reduction (%)', fontsize=12)
plt.title('Time Course of Predicted Viral Suppression', fontsize=14)
plt.grid(alpha=0.2)
sns.despine()

plt.savefig('viral_suppression.png', dpi=300)
