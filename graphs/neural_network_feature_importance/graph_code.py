plt.figure(figsize=(10, 6))
sns.barplot(x='Importance', y='Feature', data=feature_importance,
           palette='rocket', orient='h')

plt.xlabel('Mean |SHAP Value|', fontsize=12)
plt.ylabel('')
plt.title('Feature Importance for gRNA Efficiency Prediction', fontsize=14)
plt.gca().invert_yaxis()
sns.despine()

plt.savefig('feature_importance.png', dpi=300)
