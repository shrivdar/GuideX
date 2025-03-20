fig, ax1 = plt.subplots(figsize=(10, 6))

# Biomass (Left Axis)
ax1.bar([0, 1], health_data['Treated'], width=0.4,
       label='Treated', color='#4e79a7',
       yerr=health_data['SEM'][0], capsize=10)

ax1.bar([0.4, 1.4], health_data['Untreated'], width=0.4,
       label='Untreated', color='#f28e2b',
       yerr=health_data['SEM'][0], capsize=10)

ax1.set_ylabel('Biomass (%)', color='#4e79a7')
ax1.tick_params(axis='y', labelcolor='#4e79a7')

# Symptom Score (Right Axis)
ax2 = ax1.twinx()
ax2.plot([0.2, 1.2], health_data['Treated'][1:], 'D--',
        color='#e15759', markersize=10, linewidth=2)
ax2.plot([0.6, 1.6], health_data['Untreated'][1:], 's--',
        color='#76b7b2', markersize=10, linewidth=2)
ax2.set_ylabel('Symptom Severity (0-5)', color='#e15759')
ax2.tick_params(axis='y', labelcolor='#e15759')

# Formatting
ax1.set_xticks([0.2, 1.2])
ax1.set_xticklabels(health_data['Metric'])
ax1.legend(loc='upper left')
plt.title('Predicted Plant Health Outcomes', fontsize=14)

plt.savefig('plant_health.png', dpi=300)
